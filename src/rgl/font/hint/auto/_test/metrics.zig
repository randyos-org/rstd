//! Tests: stem width quantizing, and (under `-Dtest-fonts`) the style
//! coverage and Latin metrics of the reference fonts.

const std = @import("std");
const widths = @import("../widths.zig");
const coverage = @import("../coverage.zig");
const Source = @import("../Source.zig");
const LatinMetrics = @import("../latin/Metrics.zig");
const styles = @import("../data/styles.zig").styles;
const scripts = @import("../data/scripts.zig").scripts;
const Sfnt = @import("../../../sfnt/Sfnt.zig");
const test_fonts = @import("rstd_test_fonts");

const testing = std.testing;

test "autohint widths: clusters merge to a mean, and no widths keep one" {
    var t: [16]widths.Width = @splat(.{});
    for ([_]i32{ 60, 21, 58, 109, 20 }, 0..) |v, i| t[i].org = v;
    // 20 and 21 merge; as in FreeType, the next cluster starts past 58,
    // and its sum, 60, is divided by its end index, 4
    const n = widths.sortAndQuantize(&t, 5, 10);
    try testing.expectEqual(@as(usize, 4), n);
    for ([_]i32{ 20, 58, 15, 109 }, t[0..4]) |e, w| try testing.expectEqual(e, w.org);
    // FreeType's count of 0 becomes 1, its first entry kept
    var z: [16]widths.Width = @splat(.{});
    try testing.expectEqual(@as(usize, 1), widths.sortAndQuantize(&z, 0, 10));
}

test "autohint data: styles name existing scripts and stringsets" {
    for (styles) |s| try testing.expect(s.script < scripts.len);
    try testing.expectEqualStrings("hani_dflt", styles[coverage.fallback_style].name);
    try testing.expectEqualStrings("none_dflt", styles[coverage.none_style].name);
}

fn styleIndex(name: []const u8) u16 {
    for (styles, 0..) |s, i| if (std.mem.eql(u8, s.name, name)) return @intCast(i);
    unreachable;
}

const fonts = .{ "arimo", "jetbrains_mono", "ubuntu", "ibm_plex_serif", "commit_mono", "atkynson_mono" };
const monospaced = .{ "jetbrains_mono", "commit_mono", "atkynson_mono" };

test "autohint metrics: Latin coverage, widths and blue zones of the reference fonts (-Dtest-fonts)" {
    const al = testing.allocator;
    var any = false;
    inline for (fonts) |name| if (@field(test_fonts, name)) |bytes| {
        any = true;
        var arena: std.heap.ArenaAllocator = .init(al);
        defer arena.deinit();
        const sfnt = try Sfnt.init(arena.allocator(), bytes, .{});
        var src = try Source.init(al, &sfnt);
        defer src.deinit();
        const gstyles = try coverage.compute(al, &sfnt);
        defer al.free(gstyles);
        const latn = styleIndex("latn_dflt");
        // letters are Latin, digits marked
        const a = @backingInt(sfnt.getGlyphIndex('a').?);
        try testing.expectEqual(latn, gstyles[a] & coverage.style_mask);
        try testing.expect(gstyles[a] & coverage.has_cmap_entry != 0);
        try testing.expect(gstyles[@backingInt(sfnt.getGlyphIndex('7').?)] & coverage.digit != 0);
        var m = try LatinMetrics.init(al, &src, latn, gstyles);
        for (m.axis) |axis| {
            try testing.expect(axis.width_count > 0 and axis.standard_width > 0);
            try testing.expectEqual(@divTrunc(axis.standard_width, 5), axis.edge_distance_threshold);
        }
        // the Latin zones: capital and small tops and bottoms, descender
        const zones = m.axis[1].blueSlice();
        try testing.expect(zones.len >= 5);
        var x_height: usize = 0;
        for (zones) |z| {
            if (z.flags & LatinMetrics.blue_flag.adjustment != 0) x_height += 1;
            if (z.flags & LatinMetrics.blue_flag.top != 0) try testing.expect(z.shoot.org >= z.ref.org) else try testing.expect(z.shoot.org <= z.ref.org);
        }
        try testing.expectEqual(@as(usize, 1), x_height);
        inline for (monospaced) |mono| if (comptime std.mem.eql(u8, mono, name)) try testing.expect(m.digits_have_same_width);
        // at 16 ppem the x height lands on the grid
        const s = @import("../../fixed.zig").div(16 << 6, src.unitsPerEm());
        m.scale(s, s, 0, 0);
        for (m.axis[1].blueSlice()) |z| if (z.flags & LatinMetrics.blue_flag.active != 0) try testing.expectEqual(@as(i32, 0), z.ref.fit & 63);
    };
    if (!any) return error.SkipZigTest;
}
