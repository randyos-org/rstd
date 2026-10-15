//! bfloat16 conversions to and from `f32`/`f64`, and the portable bfloat16
//! kernels of `bscal` and `sbdot` (all in `kernel.zig`). The vtable wrappers
//! and their tests.
//!
//! Argument names are copied from the OpenBLAS kernel sources (see the
//! `CNAME(...)` definitions under `kernel/<arch>/*.c`). Several are named
//! `dummy*` in the reference: kernel-ABI padding slots this shape does not use.

const blas = @import("../../__root__.zig");
const types = @import("../../types.zig");

const BlasLong = types.BlasLong;
const BlasInt = types.BlasInt;
const HF16 = types.HF16;
const BF16 = types.BF16;
const XDouble = types.XDouble;
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;
const ComplexXD = types.ComplexXD;

pub fn sbstobf16_k(n: BlasLong, in: [*]f32, inc_in: BlasLong, out: [*]BF16, inc_out: BlasLong) void {
    blas.kernel(.sbstobf16_k)(n, in, inc_in, out, inc_out);
}
pub fn sbdtobf16_k(n: BlasLong, in: [*]f64, inc_in: BlasLong, out: [*]BF16, inc_out: BlasLong) void {
    blas.kernel(.sbdtobf16_k)(n, in, inc_in, out, inc_out);
}
pub fn sbf16tos_k(n: BlasLong, in: [*]BF16, inc_in: BlasLong, out: [*]f32, inc_out: BlasLong) void {
    blas.kernel(.sbf16tos_k)(n, in, inc_in, out, inc_out);
}
pub fn dbf16tod_k(n: BlasLong, in: [*]BF16, inc_in: BlasLong, out: [*]f64, inc_out: BlasLong) void {
    blas.kernel(.dbf16tod_k)(n, in, inc_in, out, inc_out);
}

// ---- portable kernels ------------------------------------------------------

const std = @import("std");
const kernel = @import("kernel.zig");
const Vtable = @import("../../Vtable.zig");

pub const toBf16 = kernel.toBf16;
pub const fromBf16 = kernel.fromBf16;
pub const sbstobf16_generic = kernel.sbstobf16;
pub const sbdtobf16_generic = kernel.sbdtobf16;
pub const sbf16tos_generic = kernel.sbf16tos;
pub const dbf16tod_generic = kernel.dbf16tod;
pub const bscal_generic = kernel.bscal_k;
pub const sbdot_generic = kernel.sbdot_portable;

/// Points the conversion, `bscal` and `sbdot` slots at the portable kernels.
pub const installGeneric = kernel.installGeneric;

/// A bfloat16 scalar as `f32`, through the vtable's conversion, as the
/// interfaces' `SBF16TOS_K(1, &alpha, 1, &alpha_float, 1)` does.
pub fn scalarToF32(v: BF16) f32 {
    var in = v;
    var out: f32 = undefined;
    blas.kernel(.sbf16tos_k)(1, @ptrCast(&in), 1, @ptrCast(&out), 1);
    return out;
}

// ---- tests -----------------------------------------------------------------

fn bitsOf(v: f32) u32 {
    return @bitCast(v);
}

test "toBf16 rounds to nearest even, flushes denormals, quiets NaNs" {
    const cases = [_][2]u32{
        .{ 0x3f800000, 0x3f80 }, // 1
        .{ 0x3f808000, 0x3f80 }, // tie, even below
        .{ 0x3f818000, 0x3f82 }, // tie, even above
        .{ 0x3f808001, 0x3f81 }, // just over half
        .{ 0x00000001, 0x0000 }, // denormal
        .{ 0x80000001, 0x8000 }, // negative denormal
        .{ 0x7f800000, 0x7f80 }, // infinity
        .{ 0xff800000, 0xff80 },
        .{ 0x7f800001, 0x7fc0 }, // signaling NaN
        .{ 0x7f7fffff, 0x7f80 }, // the largest float rounds up to infinity
    };
    for (cases) |c| {
        try std.testing.expectEqual(c[1], @backingInt(toBf16(@bitCast(c[0]))));
    }
}

test "fromBf16 widens exactly, flushes denormals, quiets NaNs" {
    const cases = [_][2]u32{
        .{ 0x4049, 0x40490000 },
        .{ 0x0001, 0x00000000 },
        .{ 0x8001, 0x80000000 },
        .{ 0xff80, 0xff800000 },
        .{ 0x7f81, 0x7fc10000 },
    };
    for (cases) |c| {
        try std.testing.expectEqual(c[1], bitsOf(fromBf16(@fromBackingInt(@intCast(c[0])))));
    }
}

test "portable conversions round-trip through strides" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    var prng = std.Random.DefaultPrng.init(1);
    const r = prng.random();
    var f: [64]f32 = undefined;
    var d: [64]f64 = undefined;
    for (&f, &d) |*a, *b| {
        a.* = (r.float(f32) - 0.5) * 1000;
        b.* = a.*;
    }
    var h: [96]BF16 = @splat(@fromBackingInt(@intCast(0)));
    var back: [64]f32 = @splat(0);
    var back_d: [64]f64 = @splat(0);
    vt.sbstobf16_k.?(32, &f, 2, &h, 3);
    vt.sbf16tos_k.?(32, &h, 3, &back, 1);
    for (0..32) |i| try std.testing.expectEqual(fromBf16(toBf16(f[2 * i])), back[i]);
    vt.sbdtobf16_k.?(32, &d, 2, &h, 3);
    vt.dbf16tod_k.?(32, &h, 3, &back_d, 2);
    for (0..32) |i| try std.testing.expectEqual(@as(f64, fromBf16(toBf16(f[2 * i]))), back_d[2 * i]);
}

const check = @import("../../blas1/check.zig");

/// Every bfloat16 and half-precision slot a reachable interface calls.
const slots = blk: {
    var out: []const []const u8 = &.{ "sbstobf16_k", "sbdtobf16_k", "sbf16tos_k", "dbf16tod_k", "bscal_k", "sbdot_k" };
    for (kernel.aliases) |a| out = out ++ [_][]const u8{a[0]};
    for (kernel.families) |f| out = out ++ [_][]const u8{ f[0] ++ "_n", f[0] ++ "_t" };
    break :blk out;
};

fn checkAll(vt: *const Vtable) !void {
    inline for (slots) |s| try std.testing.expect(@field(vt, s) != null);
    // the aliases follow the vtable's own `f32` kernels
    inline for (kernel.aliases) |a| try std.testing.expect(@field(vt, a[0]) == @field(vt, a[1]));
    inline for (kernel.aliases_l23) |a| try std.testing.expect(@field(vt, a[0]) == @field(vt, a[1]));
    // sbdot sums the widened products
    for (check.sizes) |n| for (check.strides) |inc| {
        var f: [check.capacity]f32 = undefined;
        check.fill(f32, &f, n + 3);
        var x: [check.capacity]BF16 = undefined;
        var y: [check.capacity]BF16 = undefined;
        for (0..check.capacity) |i| {
            x[i] = toBf16(f[i]);
            y[i] = toBf16(f[check.capacity - 1 - i]);
        }
        var want: f64 = 0;
        for (0..n) |i| want += @as(f64, fromBf16(x[i * inc[0]])) * fromBf16(y[i * inc[1]]);
        const got = vt.sbdot_k.?(@intCast(n), &x, @intCast(inc[0]), &y, @intCast(inc[1]));
        check.expectClose(f32, @floatCast(want), got, check.tolerance(f32, 4 * n)) catch |e| return check.case(e, "sbdot", BF16, n, inc);
    };
}

test "every bfloat16 and half-precision slot is set, on the portable kernels" {
    var vt: Vtable = .{};
    @import("../../portable.zig").install(&vt);
    try checkAll(&vt);
}

test "every bfloat16 and half-precision slot is set on every linked core, where this machine runs it" {
    try check.eachCore(checkAll);
}

test "sbdot sums a long, cancelling vector in f64" {
    // products of 2^28 against 1s: +2^28 over the first half, -2^28 over the
    // second; an f32 partial sum loses the 1s
    const n = 4096;
    var x: [n]BF16 = undefined;
    var y: [n]BF16 = undefined;
    for (0..n) |i| {
        const big = i % 2 == 0;
        x[i] = toBf16(if (big) 16384 else 1);
        y[i] = toBf16(if (!big) 1 else if (i < n / 2) 16384 else -16384);
    }
    try std.testing.expectEqual(@as(f32, n / 2), kernel.sbdot_portable(n, &x, 1, &y, 1));
    try std.testing.expectEqual(@as(f32, 0), kernel.sbdot_portable(0, &x, 1, &y, 1));
}

/// `generic/scal.c` as written, for one element.
fn bscalModel(flag: i64, da: u16, x: u16) u16 {
    const da_f = fromBf16(@fromBackingInt(@intCast(da)));
    var v = fromBf16(@fromBackingInt(@intCast(x)));
    if (flag == 0) {
        v = if (da_f == 0) 0 else da_f * v;
    } else if (da == 0) {
        v = if (!std.math.isNan(v) and !std.math.isInf(v)) 0 else std.math.nan(f32);
    } else {
        v = da_f * v;
    }
    return @backingInt(toBf16(v));
}

test "bscal follows generic/scal.c, with the flag, signed zeros and denormals" {
    const das = [_]u16{ 0x0000, 0x8000, 0x0001, 0x3fc0, 0xbf00, 0x7f80, 0x7fc0 };
    const xs = [_]u16{ 0x0000, 0x8000, 0x0001, 0x3f80, 0xc040, 0x7f80, 0xff80, 0x7fc1, 0x4049 };
    for ([_]i64{ 0, 1 }) |flag| for (das) |da| {
        var pad: [xs.len * 2]BF16 = undefined;
        for (xs, 0..) |s, i| {
            pad[2 * i] = @fromBackingInt(@intCast(s));
            pad[2 * i + 1] = @fromBackingInt(@intCast(0x1234));
        }
        _ = kernel.bscal_k(xs.len, 0, 0, @fromBackingInt(@intCast(da)), &pad, 2, &pad, 0, &pad, @intCast(flag));
        for (xs, 0..) |s, i| {
            try std.testing.expectEqual(bscalModel(flag, da, s), @backingInt(pad[2 * i]));
            try std.testing.expectEqual(@as(u16, 0x1234), @backingInt(pad[2 * i + 1]));
        }
    };
    // nothing for n <= 0 or a non-positive stride
    var one = [_]BF16{@fromBackingInt(@intCast(0x3f80))};
    _ = kernel.bscal_k(0, 0, 0, @fromBackingInt(@intCast(0)), &one, 1, &one, 0, &one, 1);
    _ = kernel.bscal_k(1, 0, 0, @fromBackingInt(@intCast(0)), &one, 0, &one, 0, &one, 1);
    try std.testing.expectEqual(@as(u16, 0x3f80), @backingInt(one[0]));
}
