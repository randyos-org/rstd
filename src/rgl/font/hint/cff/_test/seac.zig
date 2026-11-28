//! Tests: SEAC accents (`load.zig`) on a hand-made CFF table: the accent
//! drawn at its offset with its stems counted from it, then the base, as
//! if one charstring drew the accent there; and (under `-Dtest-fonts`) the
//! charset lookup of a real font.

const std = @import("std");
const fixed = @import("../../fixed.zig");
const load = @import("../load.zig").load;
const Instance = @import("../Instance.zig");
const Sfnt = @import("../../../sfnt/Sfnt.zig");
const CffTable = @import("../../../sfnt/tables/cff/CffTable.zig");
const Index = @import("../../../sfnt/tables/cff/Index.zig");
const charset = @import("../../../sfnt/tables/cff/charset.zig");
const standard_encoding = @import("../../../sfnt/tables/cff/charset_data.zig").standard_encoding;
const test_fonts = @import("rstd_test_fonts");

const testing = std.testing;

/// `v` as a charstring integer.
fn num(comptime v: i32) []const u8 {
    return comptime if (v >= -107 and v <= 107)
        &[_]u8{@intCast(v + 139)}
    else if (v >= 108 and v <= 1131)
        &[_]u8{ @intCast(247 + @divTrunc(v - 108, 256)), @intCast(@mod(v - 108, 256)) }
    else
        &[_]u8{ @intCast(251 + @divTrunc(-v - 108, 256)), @intCast(@mod(-v - 108, 256)) };
}

const hstem = [_]u8{1};
const rmoveto = [_]u8{21};
const hlineto = [_]u8{6};
const vlineto = [_]u8{7};
const endchar = [_]u8{14};

/// A box with a bottom stem.
const base = num(0) ++ num(50) ++ hstem ++ num(100) ++ num(0) ++ rmoveto ++ num(300) ++ hlineto ++ num(600) ++ vlineto ++ num(-300) ++ hlineto ++ endchar;
/// A bar with a stem, and the same bar 50 right and 20 up.
const accent = num(650) ++ num(40) ++ hstem ++ num(150) ++ num(650) ++ rmoveto ++ num(100) ++ hlineto ++ num(40) ++ vlineto ++ num(-100) ++ hlineto ++ endchar;
const accent_moved = num(670) ++ num(40) ++ hstem ++ num(200) ++ num(670) ++ rmoveto ++ num(100) ++ hlineto ++ num(40) ++ vlineto ++ num(-100) ++ hlineto ++ endchar;

// the Standard Encoding codes of the base, the accent and the composite
const base_code = 65;
const accent_code = 194;
const seac_code = 197;

const glyphs = [_][]const u8{
    &endchar,
    base,
    accent,
    num(50) ++ num(20) ++ num(base_code) ++ num(accent_code) ++ endchar,
    accent_moved,
    // a SEAC of a SEAC, and of a code past 255
    num(0) ++ num(0) ++ num(base_code) ++ num(seac_code) ++ endchar,
    num(0) ++ num(0) ++ num(base_code) ++ num(300) ++ endchar,
};

/// The glyphs' INDEX, with one-byte offsets.
fn charstrings() [3 + glyphs.len + 1 + total()]u8 {
    var out: [3 + glyphs.len + 1 + total()]u8 = undefined;
    std.mem.writeInt(u16, out[0..2], glyphs.len, .big);
    out[2] = 1;
    var off: usize = 1;
    var at: usize = 3 + glyphs.len + 1;
    for (glyphs, 0..) |g, i| {
        out[3 + i] = @intCast(off);
        @memcpy(out[at..][0..g.len], g);
        at += g.len;
        off += g.len;
    }
    out[3 + glyphs.len] = @intCast(off);
    return out;
}

fn total() usize {
    var n: usize = 0;
    for (glyphs) |g| n += g.len;
    return n;
}

/// The table: four bytes, then a format 0 charset at offset 4.
fn table() [5 + 2 * (glyphs.len - 1)]u8 {
    const sids = [_]u16{ standard_encoding[base_code], standard_encoding[accent_code], standard_encoding[seac_code], 391, 392, 393 };
    var out: [5 + 2 * (glyphs.len - 1)]u8 = @splat(0);
    for (sids, 0..) |s, i| std.mem.writeInt(u16, out[5 + 2 * i ..][0..2], s, .big);
    return out;
}

test "cff SEAC: the accent at its offset with its stems moved too, then the base" {
    const al = testing.allocator;
    const cs = comptime charstrings();
    const data = comptime table();
    const empty: Index = try .init(&.{ 0, 0 }, true);
    const top = [_]u8{ num(4)[0], 15 };
    const cff: CffTable = .{
        .table_data = &data,
        .top_dict_data = &top,
        .charstrings_index = try .init(&cs, true),
        .global_subrs_index = empty,
        .is_cid = false,
        .local_subrs_index = null,
        .fd_array_index = null,
        .fdselect_offset = null,
    };
    try testing.expectEqual(@as(?u16, 1), charset.glyphForStdCode(cff, base_code));
    try testing.expectEqual(@as(?u16, 2), charset.glyphForStdCode(cff, accent_code));
    const sfnt: Sfnt = .{ .cff = cff };
    const inst: Instance = .{ .sfnt = &sfnt, .cff = cff, .scale = fixed.div(16 << 6, 1000) };
    for ([_]bool{ false, true }) |hinted| {
        var seac = try load(al, &inst, @fromBackingInt(3), hinted);
        defer seac.deinit(al);
        var a = try load(al, &inst, @fromBackingInt(4), hinted);
        defer a.deinit(al);
        var b = try load(al, &inst, @fromBackingInt(1), hinted);
        defer b.deinit(al);
        const points = try std.mem.concat(al, [2]i32, &.{ a.points, b.points });
        defer al.free(points);
        try testing.expectEqualSlices([2]i32, points, seac.points);
        const tags = try std.mem.concat(al, u8, &.{ a.tags, b.tags });
        defer al.free(tags);
        try testing.expectEqualSlices(u8, tags, seac.tags);
        try testing.expectEqual(a.contours.len + b.contours.len, seac.contours.len);
        for (b.contours, seac.contours[a.contours.len..]) |e, got| try testing.expectEqual(e + @as(u16, @intCast(a.points.len)), got);
    }
    try testing.expectError(error.InvalidSeac, load(al, &inst, @fromBackingInt(5), true));
    try testing.expectError(error.InvalidSeac, load(al, &inst, @fromBackingInt(6), true));
}

test "cff charset: Standard Encoding codes name a real font's glyphs (-Dtest-fonts)" {
    const bytes = test_fonts.commit_mono orelse return error.SkipZigTest;
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    const sfnt = try Sfnt.init(arena.allocator(), bytes, .{});
    // `A`, `a` and the acute accent, by their glyph names
    for ([_][2]u21{ .{ 65, 'A' }, .{ 97, 'a' }, .{ 194, 0xB4 } }) |c| {
        const want = @backingInt(sfnt.getGlyphIndex(c[1]).?);
        try testing.expectEqual(@as(?u16, want), charset.glyphForStdCode(sfnt.cff.?, c[0]));
    }
}
