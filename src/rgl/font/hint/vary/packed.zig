//! A tuple variation store's packed point numbers and deltas, read as
//! FreeType 2.14.3 reads them (`src/truetype/ttgxvar.c`
//! `ft_var_readpackedpoints`, `ft_var_readpackeddeltas`): a read that runs
//! out of data is null, and the caller skips that tuple.

const std = @import("std");

const Allocator = std.mem.Allocator;

/// FreeType's stream in a frame: reads past the end give 0 (`FT_GET_*`),
/// a seek past it lands on it (`FT_Stream_SeekSet`).
pub const Cursor = struct {
    data: []const u8,
    at: usize = 0,

    pub fn seek(c: *Cursor, at: usize) void {
        c.at = @min(at, c.data.len);
    }

    pub fn left(c: *const Cursor) usize {
        return c.data.len - c.at;
    }

    pub fn byte(c: *Cursor) u8 {
        if (c.left() < 1) return 0;
        defer c.at += 1;
        return c.data[c.at];
    }

    pub fn ushort(c: *Cursor) u16 {
        if (c.left() < 2) return 0;
        defer c.at += 2;
        return std.mem.readInt(u16, c.data[c.at..][0..2], .big);
    }

    pub fn short(c: *Cursor) i16 {
        return @bitCast(c.ushort());
    }

    pub fn ulong(c: *Cursor) u32 {
        if (c.left() < 4) return 0;
        defer c.at += 4;
        return std.mem.readInt(u32, c.data[c.at..][0..4], .big);
    }
};

/// The points a tuple's deltas are for.
pub const Points = union(enum) {
    /// Every point, in order (a count of 0).
    all,
    list: []u16,

    pub fn deinit(p: Points, allocator: Allocator) void {
        switch (p) {
            .all => {},
            .list => |l| allocator.free(l),
        }
    }
};

/// Point numbers: a count, then runs of byte or word increments.
pub fn readPoints(allocator: Allocator, c: *Cursor) Allocator.Error!?Points {
    var n: usize = c.byte();
    if (n == 0) return .all;
    if (n & 0x80 != 0) n = ((n & 0x7f) << 8) | c.byte();
    const points = try allocator.alloc(u16, n);
    var first: u16 = 0;
    var i: usize = 0;
    while (i < n) {
        if (c.left() == 0) {
            allocator.free(points);
            return null;
        }
        const run = c.data[c.at];
        c.at += 1;
        // the run count leaves out its first point
        const cnt = @min(@as(usize, run & 0x7f) + 1, n - i);
        const words = run & 0x80 != 0;
        if ((if (words) 2 * cnt else cnt) > c.left()) {
            allocator.free(points);
            return null;
        }
        for (0..cnt) |_| {
            first +%= if (words) c.ushort() else c.byte();
            points[i] = first;
            i += 1;
        }
    }
    return .{ .list = points };
}

/// `count` deltas, as 16.16: runs of zeros, bytes or words.
pub fn readDeltas(allocator: Allocator, c: *Cursor, count: usize) Allocator.Error!?[]i32 {
    const deltas = try allocator.alloc(i32, count);
    var i: usize = 0;
    while (i < count) {
        if (c.left() == 0) {
            allocator.free(deltas);
            return null;
        }
        const run = c.data[c.at];
        c.at += 1;
        const cnt = @min(@as(usize, run & 0x3f) + 1, count - i);
        if (run & 0x80 != 0) {
            @memset(deltas[i..][0..cnt], 0);
            i += cnt;
            continue;
        }
        const words = run & 0x40 != 0;
        if ((if (words) 2 * cnt else cnt) > c.left()) {
            allocator.free(deltas);
            return null;
        }
        for (0..cnt) |_| {
            const v: i32 = if (words) c.short() else @as(i8, @bitCast(c.byte()));
            deltas[i] = v << 16;
            i += 1;
        }
    }
    return deltas;
}
