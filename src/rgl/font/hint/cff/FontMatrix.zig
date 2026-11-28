//! A CFF font's FontMatrix as FreeType 2.14.3 reads it
//! (`src/cff/cffparse.c` `cff_parse_font_matrix`), normalizes it so `yy`
//! is 1 with the scale in `upm`, and concatenates a CID subfont's onto the
//! top font's (`src/cff/cffobjs.c` `cff_face_init`); the loader applies
//! it (`cffgload.c` `cff_slot_load`).

const std = @import("std");
const fixed = @import("../fixed.zig");
const dict = @import("dict.zig");
const CffTable = @import("../../sfnt/tables/cff/CffTable.zig");
const GlyphIndex = @import("../../sfnt/types.zig").GlyphIndex;

const FontMatrix = @This();

pub const identity = [4]i32{ 0x10000, 0, 0, 0x10000 };

/// The DICT has a FontMatrix.
has: bool = false,
/// xx, yx, xy, yy, 16.16.
matrix: [4]i32 = identity,
/// In font units once normalized.
offset: [2]i32 = .{ 0, 0 },
upm: i32 = 1,

const power_tens = [_]i32{ 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000 };

/// The FontMatrix operator (12 7) of a DICT, as parsed; the default when
/// absent, odd or degenerate.
pub fn parse(data: []const u8) dict.Error!FontMatrix {
    var p: dict.Parser = .{ .data = data };
    var m: FontMatrix = .{};
    while (try p.next()) |item| {
        if (item.op != 0x0C07 or item.args.len < 6) continue;
        m = fromOperands(item.args[0..6]);
    }
    return m;
}

/// `cff_parse_font_matrix`: the six numbers scaled alike by the largest's power of ten.
fn fromOperands(args: []const dict.Operand) FontMatrix {
    const unlikely: FontMatrix = .{ .has = true };
    var values: [6]i32 = undefined;
    var scalings: [6]i32 = undefined;
    var max: i32 = std.math.minInt(i32);
    var min: i32 = std.math.maxInt(i32);
    for (args, &values, &scalings) |a, *v, *s| {
        v.*, s.* = a.dynamic();
        if (v.* != 0) {
            max = @max(max, s.*);
            min = @min(min, s.*);
        }
    }
    if (max < -9 or max > 0 or max - min < 0 or max - min > 9) return unlikely;
    for (&values, scalings) |*v, s| {
        if (v.* == 0) continue;
        const divisor = power_tens[@intCast(max - s)];
        const half = divisor >> 1;
        v.* = if (v.* < 0) @divTrunc(v.* - half, divisor) else @divTrunc(v.* + half, divisor);
    }
    const m: FontMatrix = .{ .has = true, .matrix = values[0..4].*, .offset = values[4..6].*, .upm = power_tens[@intCast(-max)] };
    if (!check(m.matrix)) return unlikely;
    return m;
}

/// `FT_Matrix_Check`: nonzero, and not nearly singular.
fn check(m: [4]i32) bool {
    var xx: i64 = m[0];
    var yx: i64 = m[1];
    var xy: i64 = m[2];
    var yy: i64 = m[3];
    const val: u64 = @abs(xx) | @abs(xy) | @abs(yx) | @abs(yy);
    if (val == 0 or val > 0x7FFFFFFF) return false;
    const shift: i32 = 63 - @as(i32, @clz(val)) - 12;
    if (shift > 0) {
        const s: u6 = @intCast(shift);
        xx >>= s;
        xy >>= s;
        yx >>= s;
        yy >>= s;
    }
    const temp1: u64 = 32 * @abs(xx * yy - xy * yx);
    const temp2: u64 = @bitCast(xx * xx + xy * xy + yx * yx + yy * yy);
    return temp1 > temp2;
}

/// `cff_face_init`'s normalization: `yy` (or `yx`) made 1, its scale moved
/// into `upm`, the offset to whole units.
fn normalize(m: *FontMatrix) void {
    const temp: i32 = @intCast(@abs(if (m.matrix[3] != 0) m.matrix[3] else m.matrix[1]));
    if (temp != 0x10000) {
        m.upm = fixed.div(m.upm, temp);
        for (&m.matrix) |*v| v.* = fixed.div(v.*, temp);
        for (&m.offset) |*v| v.* = fixed.div(v.*, temp);
    }
    for (&m.offset) |*v| v.* >>= 16;
}

/// What a glyph loads with: its subfont's matrix and offset, and the top
/// font's and the subfont's units per em.
pub const ForGlyph = struct { matrix: [4]i32, offset: [2]i32, top_upm: i32, sub_upm: i32 };

/// The matrices `cff_face_init` leaves for `glyph`'s font; `face_upem` is
/// the `head` table's, a top font without a matrix's.
pub fn forGlyph(cff: CffTable, face_upem: u16, glyph: GlyphIndex) dict.Error!ForGlyph {
    var top = try parse(cff.top_dict_data);
    if (!top.has) top.upm = face_upem;
    normalize(&top);
    const fd_array = cff.fd_array_index orelse return .{ .matrix = top.matrix, .offset = top.offset, .top_upm = top.upm, .sub_upm = top.upm };
    // an FD past the array is its last (`cff_slot_load`)
    const fd = @min(@as(usize, cff.getFdIndexForGlyph(glyph)), fd_array.count() -| 1);
    var sub = try parse(fd_array.get(fd));
    if (sub.has) {
        if (top.has) {
            const scaling: i32 = if (top.upm > 1 and sub.upm > 1) @min(top.upm, sub.upm) else 1;
            sub.matrix = multiplyScaled(top.matrix, sub.matrix, scaling);
            sub.offset = transformScaled(sub.offset, top.matrix, scaling);
            sub.upm = fixed.mulDiv(sub.upm, top.upm, scaling);
        }
    } else {
        sub.matrix = top.matrix;
        sub.offset = top.offset;
        sub.upm = top.upm;
    }
    // a copied top matrix is normalized again, its offset shifted again
    normalize(&sub);
    return .{ .matrix = sub.matrix, .offset = sub.offset, .top_upm = top.upm, .sub_upm = sub.upm };
}

/// `FT_Matrix_Multiply_Scaled`: `b` = `a` * `b`, over 16.16 `scaling`.
fn multiplyScaled(a: [4]i32, b: [4]i32, scaling: i32) [4]i32 {
    const val = 0x10000 *% scaling;
    const axx, const ayx, const axy, const ayy = a;
    const bxx, const byx, const bxy, const byy = b;
    return .{
        fixed.mulDiv(axx, bxx, val) +% fixed.mulDiv(axy, byx, val),
        fixed.mulDiv(ayx, bxx, val) +% fixed.mulDiv(ayy, byx, val),
        fixed.mulDiv(axx, bxy, val) +% fixed.mulDiv(axy, byy, val),
        fixed.mulDiv(ayx, bxy, val) +% fixed.mulDiv(ayy, byy, val),
    };
}

/// `FT_Vector_Transform_Scaled`.
fn transformScaled(v: [2]i32, m: [4]i32, scaling: i32) [2]i32 {
    const val = 0x10000 *% scaling;
    return .{
        fixed.mulDiv(v[0], m[0], val) +% fixed.mulDiv(v[1], m[2], val),
        fixed.mulDiv(v[0], m[1], val) +% fixed.mulDiv(v[1], m[3], val),
    };
}

/// `FT_Vector_Transform`, as `FT_Outline_Transform` moves each point.
pub fn transform(p: [2]i32, m: [4]i32) [2]i32 {
    return .{ fixed.mul(p[0], m[0]) +% fixed.mul(p[1], m[2]), fixed.mul(p[0], m[1]) +% fixed.mul(p[1], m[3]) };
}
