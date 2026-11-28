//! Reading a glyph's `glyf` record as FreeType 2.14.3 does
//! (`src/truetype/ttgload.c` `TT_Load_Glyph_Header`, `TT_Load_Simple_Glyph`,
//! `TT_Load_Composite_Glyph`), with its checks: points in font units, tags
//! kept to the on-curve bit, components with their 16.16 transforms.

const std = @import("std");

const Allocator = std.mem.Allocator;

pub const Error = error{ InvalidOutline, InvalidComposite, TooManyHints } || Allocator.Error;

pub const on_curve: u8 = 0x01;
const x_short: u8 = 0x02;
const y_short: u8 = 0x04;
const repeat: u8 = 0x08;
const x_same: u8 = 0x10;
const y_same: u8 = 0x20;

/// Component flags.
pub const args_are_words: u16 = 0x0001;
pub const args_are_xy_values: u16 = 0x0002;
pub const round_xy_to_grid: u16 = 0x0004;
pub const we_have_a_scale: u16 = 0x0008;
pub const more_components: u16 = 0x0020;
pub const we_have_an_xy_scale: u16 = 0x0040;
pub const we_have_a_2x2: u16 = 0x0080;
pub const we_have_instructions: u16 = 0x0100;
pub const use_my_metrics: u16 = 0x0200;
pub const scaled_component_offset: u16 = 0x0800;

pub const Header = struct { n_contours: i16, x_min: i16, y_min: i16, x_max: i16, y_max: i16 };

fn int(comptime T: type, b: []const u8) T {
    return std.mem.readInt(T, b[0..@sizeOf(T)], .big);
}

pub fn header(data: []const u8) Error!Header {
    if (data.len < 10) return error.InvalidOutline;
    return .{ .n_contours = int(i16, data[0..]), .x_min = int(i16, data[2..]), .y_min = int(i16, data[4..]), .x_max = int(i16, data[6..]), .y_max = int(i16, data[8..]) };
}

/// A simple glyph: contour ends, its program, on-curve tags and points in
/// font units, owned by the caller's allocator.
pub const Simple = struct {
    contours: []u16,
    instructions: []const u8,
    tags: []u8,
    points: [][2]i32,

    pub fn deinit(s: *Simple, allocator: Allocator) void {
        allocator.free(s.points);
        allocator.free(s.tags);
        allocator.free(s.contours);
    }
};

/// The simple glyph after its header; `instructions` borrows `data`.
pub fn simple(allocator: Allocator, data: []const u8, n_contours: usize) Error!Simple {
    var p: usize = 10;
    if (n_contours >= 0xFFF or p + 2 * n_contours + 2 > data.len) return error.InvalidOutline;
    const contours = try allocator.alloc(u16, n_contours);
    errdefer allocator.free(contours);
    var last: i32 = -1;
    for (contours) |*c| {
        c.* = int(u16, data[p..]);
        p += 2;
        if (c.* <= last) return error.InvalidOutline;
        last = c.*;
    }
    const n_points: usize = @intCast(last + 1);
    const n_ins = int(u16, data[p..]);
    p += 2;
    if (p + n_ins > data.len) return error.TooManyHints;
    const instructions = data[p..][0..n_ins];
    p += n_ins;
    const tags = try allocator.alloc(u8, n_points);
    errdefer allocator.free(tags);
    var f: usize = 0;
    while (f < n_points) {
        if (p + 1 > data.len) return error.InvalidOutline;
        const c = data[p];
        p += 1;
        tags[f] = c;
        f += 1;
        if (c & repeat != 0) {
            if (p + 1 > data.len) return error.InvalidOutline;
            const count = data[p];
            p += 1;
            if (f + count > n_points) return error.InvalidOutline;
            @memset(tags[f..][0..count], c);
            f += count;
        }
    }
    const points = try allocator.alloc([2]i32, n_points);
    errdefer allocator.free(points);
    inline for (0..2) |axis| {
        const short = if (axis == 0) x_short else y_short;
        const same = if (axis == 0) x_same else y_same;
        var v: i32 = 0;
        for (points, tags) |*pt, t| {
            var d: i32 = 0;
            if (t & short != 0) {
                if (p + 1 > data.len) return error.InvalidOutline;
                d = data[p];
                p += 1;
                if (t & same == 0) d = -d;
            } else if (t & same == 0) {
                if (p + 2 > data.len) return error.InvalidOutline;
                d = int(i16, data[p..]);
                p += 2;
            }
            v +%= d;
            pt[axis] = v;
        }
    }
    for (tags) |*t| t.* &= on_curve;
    return .{ .contours = contours, .instructions = instructions, .tags = tags, .points = points };
}

pub const Component = struct {
    flags: u16,
    glyph: u16,
    arg1: i32,
    arg2: i32,
    /// 16.16: x' = xx x + xy y, y' = yx x + yy y.
    xx: i32 = 0x10000,
    xy: i32 = 0,
    yx: i32 = 0,
    yy: i32 = 0x10000,
};

/// A composite glyph's components, and where its instruction count sits.
pub const Composite = struct {
    components: []Component,
    instructions_at: usize,
};

/// The components after the header; each glyph index must be below
/// `num_glyphs`.
pub fn composite(allocator: Allocator, data: []const u8, num_glyphs: u16) Error!Composite {
    var list: std.ArrayList(Component) = .empty;
    errdefer list.deinit(allocator);
    var p: usize = 10;
    while (true) {
        if (p + 4 > data.len) return error.InvalidComposite;
        var c: Component = .{ .flags = int(u16, data[p..]), .glyph = int(u16, data[p + 2 ..]), .arg1 = 0, .arg2 = 0 };
        p += 4;
        if (c.glyph >= num_glyphs) return error.InvalidComposite;
        var count: usize = 2;
        if (c.flags & args_are_words != 0) count += 2;
        if (c.flags & we_have_a_scale != 0) {
            count += 2;
        } else if (c.flags & we_have_an_xy_scale != 0) {
            count += 4;
        } else if (c.flags & we_have_a_2x2 != 0) count += 8;
        if (p + count > data.len) return error.InvalidComposite;
        const words = c.flags & args_are_words != 0;
        if (c.flags & args_are_xy_values != 0) {
            c.arg1 = if (words) int(i16, data[p..]) else @as(i8, @bitCast(data[p]));
            c.arg2 = if (words) int(i16, data[p + 2 ..]) else @as(i8, @bitCast(data[p + 1]));
        } else {
            c.arg1 = if (words) int(u16, data[p..]) else data[p];
            c.arg2 = if (words) int(u16, data[p + 2 ..]) else data[p + 1];
        }
        p += if (words) 4 else 2;
        // F2Dot14 to 16.16
        if (c.flags & we_have_a_scale != 0) {
            c.xx = @as(i32, int(i16, data[p..])) * 4;
            c.yy = c.xx;
            p += 2;
        } else if (c.flags & we_have_an_xy_scale != 0) {
            c.xx = @as(i32, int(i16, data[p..])) * 4;
            c.yy = @as(i32, int(i16, data[p + 2 ..])) * 4;
            p += 4;
        } else if (c.flags & we_have_a_2x2 != 0) {
            c.xx = @as(i32, int(i16, data[p..])) * 4;
            c.yx = @as(i32, int(i16, data[p + 2 ..])) * 4;
            c.xy = @as(i32, int(i16, data[p + 4 ..])) * 4;
            c.yy = @as(i32, int(i16, data[p + 6 ..])) * 4;
            p += 8;
        }
        try list.append(allocator, c);
        if (c.flags & more_components == 0) break;
    }
    return .{ .components = try list.toOwnedSlice(allocator), .instructions_at = p };
}
