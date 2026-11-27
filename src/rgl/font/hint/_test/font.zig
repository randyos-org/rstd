//! Tests: the hinting front end's default choice and loads against
//! FreeType's (`refs.zig`), and hinted bitmaps against unhinted ones and
//! FreeType's outlines, under `-Dtest-fonts`.

const std = @import("std");
const refs = @import("refs.zig");
const Sfnt = @import("../../sfnt/Sfnt.zig");
const Font = @import("../Font.zig");
const Glyph = @import("../Glyph.zig");
const hinted = @import("../../render/hinted.zig");
const test_fonts = @import("rstd_test_fonts");

const testing = std.testing;

const fonts = .{
    .{ "arimo", @embedFile("refs/arimo.bin.gz") },
    .{ "jetbrains_mono", @embedFile("refs/jetbrains_mono.bin.gz") },
    .{ "ubuntu", @embedFile("refs/ubuntu.bin.gz") },
    .{ "ibm_plex_serif", @embedFile("refs/ibm_plex_serif.bin.gz") },
    .{ "cica", @embedFile("refs/cica.bin.gz") },
    .{ "noto_sans_limbu", @embedFile("refs/noto_sans_limbu.bin.gz") },
    .{ "commit_mono", @embedFile("refs/commit_mono.bin.gz") },
    .{ "atkynson_mono", @embedFile("refs/atkynson_mono.bin.gz") },
};

/// The refs' mode a front-end mode and target load as.
fn refMode(m: Font.Mode, target: Font.Target) refs.Mode {
    return switch (m) {
        .none => .unhinted,
        .designer => if (target == .light) .native_light else .native,
        .auto => if (target == .light) .autohint_light else .autohint,
        .default => unreachable,
    };
}

fn same(g: Glyph, o: refs.Outline) bool {
    if (o.err != 0 or g.advance != o.advance or g.lsb_delta != o.lsb_delta or g.rsb_delta != o.rsb_delta) return false;
    if (!std.mem.eql(u16, g.contours, o.ends) or g.points.len != o.points.len) return false;
    for (g.points, o.points, g.tags, o.tags) |p, q, t, u| {
        if (p[0] != q[0] or p[1] != q[1] or t != @backingInt(u)) return false;
    }
    return true;
}

test "hinting front end: FreeType's default choice, and its default loads (-Dtest-fonts)" {
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    var any = false;
    inline for (fonts) |f| if (@field(test_fonts, f[0])) |bytes| {
        any = true;
        const r = try refs.parse(arena.allocator(), f[1]);
        const sfnt = try Sfnt.init(arena.allocator(), bytes, .{});
        var font: Font = .init(testing.allocator, &sfnt);
        defer font.deinit();
        try testing.expect(!font.tricky);
        try testing.expectEqual(r.default_mode, refMode(font.defaultMode(.normal), .normal));
        try testing.expectEqual(r.default_light_mode, refMode(font.defaultMode(.light), .light));
        // default loads, both targets, through the cached instances
        for ([_]Font.Target{ .normal, .light }) |target| {
            const mode = refMode(font.defaultMode(target), target);
            for (r.sizes, 0..) |ppem, s| for (r.glyphs) |g| {
                const o = r.outline(g, s, mode) orelse continue;
                var got = try font.load(testing.allocator, g.gid, ppem, .{ .target = target });
                defer got.deinit(testing.allocator);
                try testing.expect(same(got, o));
            };
        }
    };
    if (!any) return error.SkipZigTest;
}

/// How gray a bitmap is: each pixel's distance from fully off or on.
fn grayness(pixels: []const u8) usize {
    var n: usize = 0;
    for (pixels) |p| n += @min(p, 255 - p);
    return n;
}

test "hinted bitmaps: FreeType's box and advance, crisper than unhinted (-Dtest-fonts)" {
    const bytes = test_fonts.arimo orelse return error.SkipZigTest;
    const al = testing.allocator;
    var arena: std.heap.ArenaAllocator = .init(al);
    defer arena.deinit();
    const r = try refs.parse(arena.allocator(), fonts[0][1]);
    const sfnt = try Sfnt.init(arena.allocator(), bytes, .{});
    var font: Font = .init(al, &sfnt);
    defer font.deinit();
    for (r.glyphs) |g| {
        if (g.codepoint != 'H') continue;
        for (r.sizes, 0..) |ppem, s| {
            var h_px: std.ArrayList(u8) = .empty;
            defer h_px.deinit(al);
            var u_px: std.ArrayList(u8) = .empty;
            defer u_px.deinit(al);
            const hb = try hinted.glyphBitmapHinted(&font, al, &h_px, g.gid, ppem, .{ .mode = .designer });
            _ = try hinted.glyphBitmapHinted(&font, al, &u_px, g.gid, ppem, .{ .mode = .none });
            // FreeType's box: the control box of its hinted outline, to whole pixels
            const o = r.outline(g, s, .native).?;
            var lo = o.points[0];
            var hi = o.points[0];
            for (o.points) |p| for (0..2) |k| {
                lo[k] = @min(lo[k], p[k]);
                hi[k] = @max(hi[k], p[k]);
            };
            try testing.expectEqual(@as(i16, @intCast(lo[0] >> 6)), hb.bitmap.off_x);
            try testing.expectEqual(@as(i16, @intCast(-((hi[1] + 63) >> 6))), hb.bitmap.off_y);
            try testing.expectEqual(@as(u16, @intCast(((hi[0] + 63) >> 6) - (lo[0] >> 6))), hb.bitmap.width);
            try testing.expectEqual(@as(u16, @intCast(((hi[1] + 63) >> 6) - (lo[1] >> 6))), hb.bitmap.height);
            try testing.expectEqual(o.advance, hb.advance);
            // edges on the grid: less gray than unhinted, autohinted (both axes) the least
            var a_px: std.ArrayList(u8) = .empty;
            defer a_px.deinit(al);
            _ = try hinted.glyphBitmapHinted(&font, al, &a_px, g.gid, ppem, .{ .mode = .auto });
            try testing.expect(grayness(h_px.items) < grayness(u_px.items));
            try testing.expect(grayness(a_px.items) < grayness(u_px.items));
        }
    }
}
