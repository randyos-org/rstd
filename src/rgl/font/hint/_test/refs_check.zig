//! Tests: the FreeType references (`refs.zig`) parse and agree with
//! themselves, and with `Sfnt` on the fonts under `-Dtest-fonts`.

const std = @import("std");
const refs = @import("refs.zig");
const Sfnt = @import("../../sfnt/Sfnt.zig");
const test_fonts = @import("rstd_test_fonts");

const testing = std.testing;

/// ASCII from the space, then `hintrefs.c`'s `extra`.
const glyph_count = 0x7F - 0x20 + 46;

fn expectOutline(r: refs.Refs, mode: refs.Mode, o: refs.Outline, unhinted: refs.Outline) !void {
    try testing.expectEqual(@as(i32, 0), o.err);
    try testing.expectEqual(o.tags.len, o.points.len);
    if (o.ends.len == 0) {
        try testing.expectEqual(@as(usize, 0), o.points.len);
    } else try testing.expectEqual(o.points.len - 1, o.ends[o.ends.len - 1]);
    for (o.tags) |t| try testing.expect(t == .on or t == (if (r.cff) refs.Tag.cubic else .conic));
    // hinting moves points and keeps every contour and tag, except that the
    // CFF engine keeps a closing line it drops unhinted, when hinting moves
    // the contour's last point off its first (Atkynson Mono's `a`, say)
    if (r.cff and (mode == .native or mode == .native_light)) return;
    try testing.expectEqualSlices(u16, unhinted.ends, o.ends);
    try testing.expectEqualSlices(refs.Tag, unhinted.tags, o.tags);
}

test "hinting refs: each file parses, every glyph loads, hinting keeps the contours" {
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    inline for (refs.files) |f| {
        const r = try refs.parse(arena.allocator(), f.gz);
        try testing.expectEqual(@as(u32, 2), r.freetype[0]);
        try testing.expectEqualSlices(u16, &.{ 9, 12, 16, 24 }, r.sizes);
        try testing.expectEqual(@as(usize, glyph_count), r.glyphs.len);
        try testing.expect(r.units_per_em > 0);
        // each mode hints: some glyph moves from its unhinted outline
        var moved: [6]bool = @splat(false);
        for (r.glyphs) |g| for (0..r.sizes.len) |s| {
            const unhinted = r.outline(g, s, .unhinted) orelse continue;
            for (r.modes, 0..) |m, i| {
                const o = r.outline(g, s, m).?;
                try expectOutline(r, m, o, unhinted);
                if (!std.mem.eql([2]i32, o.points, unhinted.points)) moved[i] = true;
            }
        };
        // native hinting moves something in some mode when the face has
        // hints; not in each (Cica's programs change nothing at v40)
        var native_moved = false;
        for (r.modes, 0..) |m, i| switch (m) {
            .unhinted => try testing.expect(!moved[i]),
            .native, .native_v35, .native_light => native_moved = native_moved or moved[i],
            // the autohinter leaves a synthetic font's glyphs, of no script, alone
            .autohint, .autohint_light => if (!f.synthetic) try testing.expect(moved[i]),
        };
        try testing.expectEqual(r.cff or r.fpgm_len + r.prep_len > 0, native_moved);
    }
}

test "hinting refs: glyph indexes, units per em and outline format match Sfnt (-Dtest-fonts)" {
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    var any = false;
    inline for (refs.files) |f| if (@field(test_fonts, f.font)) |bytes| {
        any = true;
        const r = try refs.parse(arena.allocator(), f.gz);
        const sfnt = try Sfnt.init(arena.allocator(), bytes, .{});
        try testing.expectEqual(r.units_per_em, sfnt.head.?.unitsPerEm());
        try testing.expectEqual(r.cff, sfnt.cff != null);
        for (r.glyphs) |g| {
            const gid = if (sfnt.getGlyphIndex(g.codepoint)) |i| @backingInt(i) else 0;
            try testing.expectEqual(g.gid, gid);
        }
    };
    if (!any) return error.SkipZigTest;
}
