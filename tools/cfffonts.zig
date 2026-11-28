//! `cfffonts <cid|seac> <out.otf>`: writes a synthetic CFF test font for the hinting tests and
//! `hint-refs` (issue #832), as no real font at hand has these features:
//! - `cid`: CID-keyed; FD 0's FontMatrix the identity, FD 1's half scale (2000 units per em), sheared
//!   and offset, under a top FontMatrix of 1/1000; boxes and curves with stems in both FDs.
//! - `seac`: name-keyed, sheared and offset by its top FontMatrix; a base, an accent, and two SEAC
//!   composites of them, one with a width.
//! Both have `A` to `D` or `A`, '`', `À` and `Å` in `cmap`, 1000 units per em and 600-unit advances.

const std = @import("std");

const Allocator = std.mem.Allocator;
const Buf = std.ArrayList(u8);

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    const arena = init.arena.allocator();
    const args = try init.minimal.args.toSlice(arena);
    if (args.len != 3) {
        std.debug.print("usage: cfffonts <cid|seac> <out.otf>\n", .{});
        std.process.exit(2);
    }
    const font = if (std.mem.eql(u8, args[1], "cid")) try cidFont(arena) else if (std.mem.eql(u8, args[1], "seac")) try seacFont(arena) else {
        std.debug.print("cfffonts: unknown font {s}\n", .{args[1]});
        std.process.exit(2);
    };
    try std.Io.Dir.cwd().writeFile(io, .{ .sub_path = args[2], .data = font });
}

// charstring operators
const hstem = 1;
const vstem = 3;
const rlineto = 5;
const hlineto = 6;
const vlineto = 7;
const rrcurveto = 8;
const endchar = 14;
const rmoveto = 21;

/// A Type 2 charstring: numbers, and operators as `.{op}`.
fn charstring(a: Allocator, items: []const []const i32) ![]u8 {
    var b: Buf = .empty;
    for (items) |it| {
        for (it[0 .. it.len - 1]) |v| try csNumber(a, &b, v);
        try b.append(a, @intCast(it[it.len - 1]));
    }
    return b.items;
}

fn csNumber(a: Allocator, b: *Buf, v: i32) !void {
    if (v >= -107 and v <= 107) return b.append(a, @intCast(v + 139));
    if (v >= 108 and v <= 1131) return b.appendSlice(a, &.{ @intCast(247 + @divTrunc(v - 108, 256)), @intCast(@mod(v - 108, 256)) });
    if (v >= -1131 and v <= -108) return b.appendSlice(a, &.{ @intCast(251 + @divTrunc(-v - 108, 256)), @intCast(@mod(-v - 108, 256)) });
    try b.append(a, 28);
    try b.appendSlice(a, &std.mem.toBytes(std.mem.nativeToBig(i16, @intCast(v))));
}

/// A DICT integer, shortest form.
fn dictInt(a: Allocator, b: *Buf, v: i32) !void {
    if (v >= -1131 and v <= 1131) return csNumber(a, b, v);
    try dictInt32(a, b, v);
}

/// A DICT integer in 5 bytes, for offsets patched later.
fn dictInt32(a: Allocator, b: *Buf, v: i32) !void {
    try b.append(a, 29);
    try b.appendSlice(a, &std.mem.toBytes(std.mem.nativeToBig(i32, v)));
}

/// A DICT real from its decimal text (digits, `.`, `-`, `E`, `E-`).
fn dictReal(a: Allocator, b: *Buf, text: []const u8) !void {
    var nibs: std.ArrayList(u4) = .empty;
    var i: usize = 0;
    while (i < text.len) : (i += 1) {
        const c = text[i];
        try nibs.append(a, switch (c) {
            '0'...'9' => @intCast(c - '0'),
            '.' => 0xA,
            '-' => 0xE,
            'E' => if (i + 1 < text.len and text[i + 1] == '-') blk: {
                i += 1;
                break :blk 0xC;
            } else 0xB,
            else => unreachable,
        });
    }
    try nibs.append(a, 0xF);
    if (nibs.items.len % 2 != 0) try nibs.append(a, 0xF);
    try b.append(a, 30);
    var k: usize = 0;
    while (k < nibs.items.len) : (k += 2) try b.append(a, (@as(u8, nibs.items[k]) << 4) | nibs.items[k + 1]);
}

fn dictOp(a: Allocator, b: *Buf, op: u16) !void {
    if (op >= 0x0C00) try b.append(a, 12);
    try b.append(a, @truncate(op));
}

/// An INDEX of `items`.
fn index(a: Allocator, b: *Buf, items: []const []const u8) !void {
    try b.appendSlice(a, &std.mem.toBytes(std.mem.nativeToBig(u16, @intCast(items.len))));
    if (items.len == 0) return;
    var total: usize = 1;
    for (items) |it| total += it.len;
    const off_size: u8 = if (total <= 0xFF) 1 else if (total <= 0xFFFF) 2 else 4;
    try b.append(a, off_size);
    var off: usize = 1;
    for (0..items.len + 1) |k| {
        const v: u32 = @intCast(off);
        try b.appendSlice(a, (&std.mem.toBytes(std.mem.nativeToBig(u32, v)))[4 - off_size ..]);
        if (k < items.len) off += items[k].len;
    }
    for (items) |it| try b.appendSlice(a, it);
}

fn patch32(b: []u8, at: usize, v: usize) void {
    std.mem.writeInt(i32, b[at + 1 ..][0..4], @intCast(v), .big);
}

fn u16be(a: Allocator, b: *Buf, v: u16) !void {
    try b.appendSlice(a, &std.mem.toBytes(std.mem.nativeToBig(u16, v)));
}

/// A Private DICT: blue zones of the given unit, stem widths.
fn privateDict(a: Allocator, unit: i32) ![]u8 {
    var b: Buf = .empty;
    for ([_]i32{ -15, 15, 700, 15 }) |v| try dictInt(a, &b, v * unit);
    try dictOp(a, &b, 6);
    try dictInt(a, &b, 50 * unit);
    try dictOp(a, &b, 10);
    try dictInt(a, &b, 60 * unit);
    try dictOp(a, &b, 11);
    return b.items;
}

/// The box and the curve, in `unit`s of FD 0's size.
fn box(a: Allocator, unit: i32) ![]u8 {
    return charstring(a, &.{
        &.{ 0, 50 * unit, hstem },    &.{ 100 * unit, 60 * unit, vstem },
        &.{ 100 * unit, 0, rmoveto }, &.{ 400 * unit, hlineto },
        &.{ 700 * unit, vlineto },    &.{ -400 * unit, hlineto },
        &.{endchar},
    });
}

fn curve(a: Allocator, unit: i32) ![]u8 {
    return charstring(a, &.{
        &.{ 0, 40 * unit, hstem },     &.{ 100 * unit, 0, rmoveto },
        &.{ 300 * unit, 0, rlineto },  &.{ 100 * unit, 100 * unit, 0, 300 * unit, -100 * unit, 100 * unit, rrcurveto },
        &.{ -300 * unit, 0, rlineto }, &.{endchar},
    });
}

fn cidFont(a: Allocator) ![]u8 {
    const glyphs = [_][]const u8{ try charstring(a, &.{&.{endchar}}), try box(a, 1), try box(a, 2), try curve(a, 2), try curve(a, 1) };
    // FDSelect: glyphs 0 and 1 in FD 0, 2 and 3 in FD 1, 4 in FD 0
    const fd_ranges = [_][2]u16{ .{ 0, 0 }, .{ 2, 1 }, .{ 4, 0 } };
    var cff: Buf = .empty;
    try cff.appendSlice(a, &.{ 1, 0, 4, 4 });
    try index(a, &cff, &.{"RstdCid"});
    // the Top DICT, its offsets patched once known
    var top: Buf = .empty;
    for ([_]i32{ 391, 392, 0 }) |v| try dictInt(a, &top, v);
    try dictOp(a, &top, 0x0C1E);
    for ([_][]const u8{ "0.001", "0", "0", "0.001", "0", "0" }) |t| try dictReal(a, &top, t);
    try dictOp(a, &top, 0x0C07);
    try dictInt(a, &top, glyphs.len);
    try dictOp(a, &top, 0x0C22);
    const at_charset = top.items.len;
    try dictInt32(a, &top, 0);
    try dictOp(a, &top, 15);
    const at_fdselect = top.items.len;
    try dictInt32(a, &top, 0);
    try dictOp(a, &top, 0x0C25);
    const at_charstrings = top.items.len;
    try dictInt32(a, &top, 0);
    try dictOp(a, &top, 17);
    const at_fdarray = top.items.len;
    try dictInt32(a, &top, 0);
    try dictOp(a, &top, 0x0C24);
    const top_at = cff.items.len + 2 + 1 + 1 + 1; // the INDEX header with one-byte offsets
    try index(a, &cff, &.{top.items});
    try index(a, &cff, &.{ "Adobe", "Identity", "RstdCid-FD0", "RstdCid-FD1" });
    try index(a, &cff, &.{});
    // charset: CIDs 1 to 4
    const charset_at = cff.items.len;
    try cff.append(a, 2);
    try u16be(a, &cff, 1);
    try u16be(a, &cff, glyphs.len - 2);
    const fdselect_at = cff.items.len;
    try cff.append(a, 3);
    try u16be(a, &cff, fd_ranges.len);
    for (fd_ranges) |r| {
        try u16be(a, &cff, r[0]);
        try cff.append(a, @intCast(r[1]));
    }
    try u16be(a, &cff, glyphs.len);
    const charstrings_at = cff.items.len;
    try index(a, &cff, &glyphs);
    // the FDs' Private DICTs, then the FDArray naming them
    const privates = [_][]u8{ try privateDict(a, 1), try privateDict(a, 2) };
    var private_at: [2]usize = undefined;
    for (privates, &private_at) |p, *at| {
        at.* = cff.items.len;
        try cff.appendSlice(a, p);
    }
    var fds: [2][]u8 = undefined;
    const matrices = [2][6][]const u8{
        .{ "1", "0", "0", "1", "0", "0" },
        .{ "0.5", "0", "0.1", "0.5", "0.02", "0" },
    };
    for (&fds, 0..) |*fd, k| {
        var d: Buf = .empty;
        try dictInt(a, &d, @intCast(393 + k));
        try dictOp(a, &d, 0x0C26);
        for (matrices[k]) |t| try dictReal(a, &d, t);
        try dictOp(a, &d, 0x0C07);
        try dictInt(a, &d, @intCast(privates[k].len));
        try dictInt32(a, &d, @intCast(private_at[k]));
        try dictOp(a, &d, 18);
        fd.* = d.items;
    }
    const fdarray_at = cff.items.len;
    try index(a, &cff, &fds);
    patch32(cff.items, top_at + at_charset, charset_at);
    patch32(cff.items, top_at + at_fdselect, fdselect_at);
    patch32(cff.items, top_at + at_charstrings, charstrings_at);
    patch32(cff.items, top_at + at_fdarray, fdarray_at);
    return sfnt(a, cff.items, glyphs.len, &.{ .{ 'A', 1 }, .{ 'B', 2 }, .{ 'C', 3 }, .{ 'D', 4 } });
}

fn seacFont(a: Allocator) ![]u8 {
    const grave = try charstring(a, &.{
        &.{ 760, 40, hstem }, &.{ 200, 760, rmoveto }, &.{ 150, hlineto },
        &.{ 40, vlineto },    &.{ -150, hlineto },     &.{endchar},
    });
    const glyphs = [_][]const u8{
        try charstring(a, &.{&.{endchar}}),
        try box(a, 1),
        grave,
        // Agrave: the grave, raised and moved right, over A; Aring the same with a width
        try charstring(a, &.{&.{ 50, 60, 65, 193, endchar }}),
        try charstring(a, &.{&.{ 600, 100, 80, 65, 193, endchar }}),
    };
    // A, grave (by their Standard Encoding names), then two names of no code
    const sids = [_]u16{ 34, 124, 174, 175 };
    var cff: Buf = .empty;
    try cff.appendSlice(a, &.{ 1, 0, 4, 4 });
    try index(a, &cff, &.{"RstdSeac"});
    var top: Buf = .empty;
    for ([_][]const u8{ "0.001", "0", "0.00025", "0.001", "0.01", "0" }) |t| try dictReal(a, &top, t);
    try dictOp(a, &top, 0x0C07);
    const at_charset = top.items.len;
    try dictInt32(a, &top, 0);
    try dictOp(a, &top, 15);
    const at_charstrings = top.items.len;
    try dictInt32(a, &top, 0);
    try dictOp(a, &top, 17);
    const priv = try privateDict(a, 1);
    try dictInt(a, &top, @intCast(priv.len));
    const at_private = top.items.len;
    try dictInt32(a, &top, 0);
    try dictOp(a, &top, 18);
    const top_at = cff.items.len + 2 + 1 + 1 + 1;
    try index(a, &cff, &.{top.items});
    try index(a, &cff, &.{});
    try index(a, &cff, &.{});
    const charset_at = cff.items.len;
    try cff.append(a, 0);
    for (sids) |s| try u16be(a, &cff, s);
    const charstrings_at = cff.items.len;
    try index(a, &cff, &glyphs);
    const private_at = cff.items.len;
    try cff.appendSlice(a, priv);
    patch32(cff.items, top_at + at_charset, charset_at);
    patch32(cff.items, top_at + at_charstrings, charstrings_at);
    patch32(cff.items, top_at + at_private, private_at);
    return sfnt(a, cff.items, glyphs.len, &.{ .{ 'A', 1 }, .{ '`', 2 }, .{ 0xC0, 3 }, .{ 0xC5, 4 } });
}

/// An OpenType font of `cff` and the tables FreeType needs: `head`, `hhea`,
/// `hmtx`, `maxp`, and a format 4 `cmap` of `map` (codepoint, glyph).
fn sfnt(a: Allocator, cff: []const u8, n_glyphs: usize, map: []const [2]u16) ![]u8 {
    var head: Buf = .empty;
    for ([_]u16{ 1, 0, 1, 0, 0, 0, 0x5F0F, 0x3CF5, 3, 1000 }) |v| try u16be(a, &head, v);
    try head.appendSlice(a, &@as([16]u8, @splat(0)));
    for ([_]u16{ 0, @bitCast(@as(i16, -200)), 1000, 900, 0, 8, 2, 0, 0 }) |v| try u16be(a, &head, v);
    var hhea: Buf = .empty;
    for ([_]u16{ 1, 0, 800, @bitCast(@as(i16, -200)), 0, 600, 0, 0, 1000, 1, 0, 0, 0, 0, 0, 0, 0 }) |v| try u16be(a, &hhea, v);
    try u16be(a, &hhea, @intCast(n_glyphs));
    var hmtx: Buf = .empty;
    for (0..n_glyphs) |_| {
        try u16be(a, &hmtx, 600);
        try u16be(a, &hmtx, 0);
    }
    var maxp: Buf = .empty;
    try maxp.appendSlice(a, &.{ 0, 0, 0x50, 0 });
    try u16be(a, &maxp, @intCast(n_glyphs));
    // cmap: one segment per codepoint, then 0xFFFF
    var cmap: Buf = .empty;
    const seg = map.len + 1;
    for ([_]u16{ 0, 1, 3, 1 }) |v| try u16be(a, &cmap, v);
    try cmap.appendSlice(a, &.{ 0, 0, 0, 12 });
    // the binary search fields: the largest power of two segments
    const pow: u16 = @as(u16, 1) << @intCast(std.math.log2_int(usize, seg));
    for ([_]u16{ 4, @intCast(16 + 8 * seg), 0, @intCast(2 * seg), 2 * pow, std.math.log2_int(u16, pow), @intCast(2 * seg - 2 * pow) }) |v| try u16be(a, &cmap, v);
    for (map) |m| try u16be(a, &cmap, m[0]);
    try u16be(a, &cmap, 0xFFFF);
    try u16be(a, &cmap, 0);
    for (map) |m| try u16be(a, &cmap, m[0]);
    try u16be(a, &cmap, 0xFFFF);
    for (map) |m| try u16be(a, &cmap, m[1] -% m[0]);
    try u16be(a, &cmap, 1);
    for (0..seg) |_| try u16be(a, &cmap, 0);
    const tables = [_]struct { tag: *const [4]u8, data: []const u8 }{
        .{ .tag = "CFF ", .data = cff },
        .{ .tag = "cmap", .data = cmap.items },
        .{ .tag = "head", .data = head.items },
        .{ .tag = "hhea", .data = hhea.items },
        .{ .tag = "hmtx", .data = hmtx.items },
        .{ .tag = "maxp", .data = maxp.items },
    };
    var out: Buf = .empty;
    try out.appendSlice(a, "OTTO");
    for ([_]u16{ tables.len, 64, 2, 32 }) |v| try u16be(a, &out, v);
    var at: usize = 12 + 16 * tables.len;
    for (tables) |t| {
        try out.appendSlice(a, t.tag);
        try out.appendSlice(a, &std.mem.toBytes(std.mem.nativeToBig(u32, checksum(t.data))));
        try out.appendSlice(a, &std.mem.toBytes(std.mem.nativeToBig(u32, @intCast(at))));
        try out.appendSlice(a, &std.mem.toBytes(std.mem.nativeToBig(u32, @intCast(t.data.len))));
        at += (t.data.len + 3) & ~@as(usize, 3);
    }
    for (tables) |t| {
        try out.appendSlice(a, t.data);
        while (out.items.len % 4 != 0) try out.append(a, 0);
    }
    return out.items;
}

fn checksum(data: []const u8) u32 {
    var sum: u32 = 0;
    var i: usize = 0;
    while (i < data.len) : (i += 4) {
        var w: [4]u8 = @splat(0);
        const n = @min(4, data.len - i);
        @memcpy(w[0..n], data[i..][0..n]);
        sum +%= std.mem.readInt(u32, &w, .big);
    }
    return sum;
}
