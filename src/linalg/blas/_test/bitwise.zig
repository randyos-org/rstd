//! Bit-for-bit comparison: value bits without an `f80`'s padding, their bytes
//! and digests for reproducibility checks, and the all-ones poison.

const std = @import("std");
const sysinfo = @import("builtin");
const zigstd = @import("../../../zigstd.zig");
const math = @import("../../../math.zig");
const types = @import("../types.zig");
const check = @import("check.zig");

const isComplex = types.isComplex;
const FloatOf = types.Float;

/// Whether two slices hold the same bytes.
pub fn sameBits(comptime T: type, a: []const T, b: []const T) bool {
    return std.mem.eql(u8, std.mem.sliceAsBytes(a), std.mem.sliceAsBytes(b));
}

/// `sameBits`, but for `f64` on x87 within a few ulps of the largest entry:
/// there the compiler picks, per code path, where excess precision is
/// rounded off (issue #264).
pub fn sameBitsOrX87(comptime T: type, a: []const T, b: []const T) bool {
    if (comptime !(math.x87_f64 and FloatOf(T) == f64)) return sameBits(T, a, b);
    var big: f64 = 0;
    for (a) |x| big = @max(big, mag(x));
    const tol = 8 * std.math.floatEps(f64) * big;
    for (a, b) |x, y| {
        const d = if (comptime isComplex(T)) @abs(x.re - y.re) + @abs(x.im - y.im) else @abs(x - y);
        if (!(d <= tol)) return false;
    }
    return true;
}

fn mag(x: anytype) f64 {
    return if (comptime isComplex(@TypeOf(x))) @abs(x.re) + @abs(x.im) else @abs(x);
}

/// Fails unless `want` and `got`, a scalar or a `@Vector`, have the same bits.
pub fn expectSameBits(want: anytype, got: @TypeOf(want)) !void {
    const T = @TypeOf(want);
    const U = switch (@typeInfo(T)) {
        .vector => |info| @Vector(info.len, @Int(.unsigned, @bitSizeOf(info.child))),
        else => @Int(.unsigned, @bitSizeOf(T)),
    };
    try zigstd.testExpectEqual(@as(U, @bitCast(want)), @as(U, @bitCast(got)));
}

fn Bits(comptime F: type) type {
    return @Int(.unsigned, @bitSizeOf(F));
}

/// The value bits of `a` and `b` agree (not an `f80`'s padding); with
/// `signed_zero` false, zeros of either sign agree too.
pub fn same(comptime T: type, a: T, b: T, signed_zero: bool) bool {
    if (comptime isComplex(T)) return same(FloatOf(T), a.re, b.re, signed_zero) and same(FloatOf(T), a.im, b.im, signed_zero);
    if (!signed_zero and a == 0 and b == 0) return true;
    return @as(Bits(T), @bitCast(a)) == @as(Bits(T), @bitCast(b));
}

/// Appends the value bits of `v` (a scalar, a complex, or a slice of them)
/// to `out`, not an `f80`'s padding, which Debug leaves undefined.
pub fn appendBits(al: std.mem.Allocator, out: *std.ArrayList(u8), v: anytype) std.mem.Allocator.Error!void {
    const V = @TypeOf(v);
    switch (@typeInfo(V)) {
        .pointer => |p| {
            if (comptime p.size != .slice) return appendBits(al, out, v.*);
            if (comptime !padded(p.child)) return out.appendSlice(al, std.mem.sliceAsBytes(v));
            for (v) |e| try appendBits(al, out, e);
        },
        .float => {
            var b: [bitsLen(V)]u8 = undefined;
            std.mem.writeInt(Bits(V), &b, @bitCast(v), sysinfo.target.cpu.arch.endian());
            try out.appendSlice(al, &b);
        },
        .@"struct" => if (comptime isComplex(V)) {
            try appendBits(al, out, v.re);
            try appendBits(al, out, v.im);
        } else try out.appendSlice(al, std.mem.asBytes(&v)),
        else => try out.appendSlice(al, std.mem.asBytes(&v)),
    }
}

/// The bytes `appendBits` gives an element of `T`.
pub fn bitsLen(comptime T: type) usize {
    if (comptime !padded(T)) return @sizeOf(T);
    const n = @divExact(@bitSizeOf(FloatOf(T)), 8);
    return if (isComplex(T)) 2 * n else n;
}

/// An element of `T` back from `appendBits`' bytes `b` (`bitsLen(T)` of them).
pub fn fromBits(comptime T: type, b: []const u8) T {
    if (comptime !padded(T)) return std.mem.bytesToValue(T, b[0..@sizeOf(T)]);
    if (comptime isComplex(T)) {
        const n = bitsLen(FloatOf(T));
        return .init(fromBits(FloatOf(T), b[0..n]), fromBits(FloatOf(T), b[n..][0..n]));
    }
    return @bitCast(std.mem.readInt(Bits(T), b[0..comptime bitsLen(T)], sysinfo.target.cpu.arch.endian()));
}

/// A float, or a complex of one, whose storage has bytes past its value.
fn padded(comptime T: type) bool {
    const F = if (isComplex(T)) FloatOf(T) else T;
    return @typeInfo(F) == .float and @bitSizeOf(F) != 8 * @sizeOf(F);
}

/// Appends an 8-byte digest of `appendBits`' bytes of `v`: for outputs too
/// large for a rig's heap to keep twice.
pub fn appendDigest(al: std.mem.Allocator, out: *std.ArrayList(u8), v: anytype) std.mem.Allocator.Error!void {
    var bits: std.ArrayList(u8) = .empty;
    defer bits.deinit(al);
    try appendBits(al, &bits, v);
    const d = std.hash.Wyhash.hash(0, bits.items);
    try out.appendSlice(al, std.mem.asBytes(&d));
}

/// All ones: what `Vector` and `Matrix` fill unread storage with.
pub fn poisoned(comptime T: type, a: T) bool {
    if (comptime isComplex(T)) return poisoned(FloatOf(T), a.re) and poisoned(FloatOf(T), a.im);
    return @as(Bits(T), @bitCast(a)) == std.math.maxInt(Bits(T));
}

/// `got` has exactly `want`'s bits.
pub fn expectBits(comptime F: type, want: F, got: F, what: []const u8, n: usize, inc: usize) !void {
    const U = @Int(.unsigned, @bitSizeOf(F));
    if (@as(U, @bitCast(want)) == @as(U, @bitCast(got))) return;
    check.say("{s} {s} n={d} inc={d}: got {e}, want {e}\n", .{ what, @typeName(F), n, inc, got, want });
    return error.TestExpectedEqual;
}

/// `std.mem.eql` on bytes a word at a time: its vector path is slow in debug
/// builds, and the race stress tests compare megabytes per round.
pub fn sameBytes(a: []const u8, b: []const u8) bool {
    if (a.len != b.len) return false;
    const words = a.len / 8 * 8;
    const x = std.mem.bytesAsSlice(u64, a[0..words]);
    const y = std.mem.bytesAsSlice(u64, b[0..words]);
    for (x, y) |p, q| if (p != q) return false;
    for (a[words..], b[words..]) |p, q| if (p != q) return false;
    return true;
}

test "appendBits leaves out an f80's padding, and fromBits reads it back" {
    const al = std.testing.allocator;
    inline for (.{ f64, f80, std.math.Complex(f80) }) |T| {
        var out: std.ArrayList(u8) = .empty;
        defer out.deinit(al);
        var v: [2]T = undefined;
        // padding bytes of 0xAA, as Debug leaves them, must not reach `out`
        @memset(std.mem.asBytes(&v), 0xaa);
        for (&v, [_]f80{ 1.5, -2.25 }) |*e, x| e.* = if (comptime isComplex(T)) .init(x, -x) else @floatCast(x);
        try appendBits(al, &out, @as([]const T, &v));
        const z = bitsLen(T);
        try zigstd.testExpectEqual(2 * z, out.items.len);
        try zigstd.testExpectEqual(@as(usize, if (T == f64) 8 else if (T == f80) 10 else 20), z);
        for (v, 0..) |e, i| try zigstd.testExpectTrue(same(T, e, fromBits(T, out.items[i * z ..][0..z]), true));
    }
}
