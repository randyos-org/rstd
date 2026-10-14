//! bfloat16 conversions to and from `f32`/`f64`. The vtable wrappers, the
//! portable kernels (`*_generic`) and their tests.
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
const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");

/// `kernel/x86_64/tobf16.c`: denormals flush to a signed zero, NaNs come out
/// quiet, and the rest round to nearest even (overflowing to infinity).
pub fn toBf16(v: f32) BF16 {
    var bits: u32 = @bitCast(v);
    const top: u16 = switch (bits & 0xff800000) {
        0x00000000 => 0x0000,
        0x80000000 => 0x8000,
        0x7f800000, 0xff800000 => @as(u16, @truncate(bits >> 16)) | @as(u16, if (bits & 0x007fffff != 0) 0x0040 else 0),
        else => blk: {
            bits +%= ((bits >> 16) & 1) + 0x7fff;
            break :blk @truncate(bits >> 16);
        },
    };
    return @fromBackingInt(@intCast(top));
}

/// `kernel/x86_64/bf16to.c`: denormals flush to a signed zero and NaNs come
/// out quiet; the rest widen exactly.
pub fn fromBf16(b: BF16) f32 {
    const h: u16 = @backingInt(b);
    const top: u16 = switch (h & 0xff80) {
        0x0000 => 0x0000,
        0x8000 => 0x8000,
        0x7f80, 0xff80 => h | @as(u16, if (h & 0x007f != 0) 0x0040 else 0),
        else => h,
    };
    return @bitCast(@as(u32, top) << 16);
}

/// To bfloat16 from `F`; a double rounds to float first, as upstream does.
fn To(comptime F: type) type {
    return struct {
        fn k(n: BlasLong, in: [*]F, inc_in: BlasLong, out: [*]BF16, inc_out: BlasLong) callconv(.c) void {
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                util.at(BF16, out, i * inc_out).* = toBf16(@floatCast(util.at(F, in, i * inc_in).*));
            }
        }
    };
}

/// From bfloat16 to `F`.
fn From(comptime F: type) type {
    return struct {
        fn k(n: BlasLong, in: [*]BF16, inc_in: BlasLong, out: [*]F, inc_out: BlasLong) callconv(.c) void {
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                util.at(F, out, i * inc_out).* = fromBf16(util.at(BF16, in, i * inc_in).*);
            }
        }
    };
}

pub const sbstobf16_generic = &To(f32).k;
pub const sbdtobf16_generic = &To(f64).k;
pub const sbf16tos_generic = &From(f32).k;
pub const dbf16tod_generic = &From(f64).k;

/// Points every conversion slot at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.sbstobf16_k = sbstobf16_generic;
    vt.sbdtobf16_k = sbdtobf16_generic;
    vt.sbf16tos_k = sbf16tos_generic;
    vt.dbf16tod_k = dbf16tod_generic;
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
