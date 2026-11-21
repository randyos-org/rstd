//! The elementwise drivers: `y[i] = f(x[i])` and `z[i] = f(x[i], y[i])`
//! over strided views, `@Vector` blocks where every view has stride 1, a
//! scalar loop otherwise. Threaded by chunks that start on the BLAS's
//! quantum, so each element takes the path, and gets the bits, one call
//! would give it (#614).

const std = @import("std");
const types = @import("../blas/types.zig");
const threading = @import("../blas/threading.zig");
const strided = @import("../blas/strided.zig");
const shared = @import("../blas/vector/util.zig");
const core = @import("../blas/blas1/interface.zig");
const scalar = @import("../../scalar.zig");

const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const ConstStrided = strided.ConstStrided;
const Strided = strided.Strided;

pub const Error = strided.LengthError || Io.Cancelable;

/// The one-input functions.
pub const Unary = enum { sqrt, exp, exp2, log, log2, log10, sin, cos, tan, abs, floor, ceil, trunc, round };

/// The two-input functions.
pub const Binary = enum { add, sub, mul, div };

/// Whether `T` goes through `@Vector` blocks: the IEEE types the targets'
/// vectors hold, not `f64` on the x87.
fn vectorized(comptime T: type) bool {
    return switch (T) {
        f16, f32 => true,
        f64 => !scalar.x87_f64,
        else => false,
    };
}

/// `op` of a real `v`, a scalar or a vector.
inline fn unaryReal(comptime op: Unary, v: anytype) @TypeOf(v) {
    return switch (op) {
        .sqrt => @sqrt(v),
        .exp => @exp(v),
        .exp2 => @exp2(v),
        .log => @log(v),
        .log2 => @log2(v),
        .log10 => @log10(v),
        .sin => @sin(v),
        .cos => @cos(v),
        .tan => @tan(v),
        .abs => @abs(v),
        // exact on the x87 too (issue #264)
        .floor => if (@typeInfo(@TypeOf(v)) == .vector) @floor(v) else scalar.floor(v),
        .ceil => if (@typeInfo(@TypeOf(v)) == .vector) @ceil(v) else scalar.ceil(v),
        .trunc => @trunc(v),
        .round => if (@typeInfo(@TypeOf(v)) == .vector) @round(v) else scalar.round(v),
    };
}

/// `op` of one element: complex `abs` by `std.math.hypot`, no overflow.
inline fn unaryOne(comptime T: type, comptime op: Unary, v: T) Out(T, op) {
    if (comptime isComplex(T)) return std.math.hypot(v.re, v.im);
    return unaryReal(op, v);
}

inline fn binaryOne(comptime T: type, comptime op: Binary, a: T, b: T) T {
    if (comptime isComplex(T)) return switch (op) {
        .add => a.add(b),
        .sub => a.sub(b),
        .mul => a.mul(b),
        .div => a.div(b),
    };
    return binaryReal(op, a, b);
}

inline fn binaryReal(comptime op: Binary, a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    return switch (op) {
        .add => a + b,
        .sub => a - b,
        .mul => a * b,
        .div => a / b,
    };
}

/// `op`'s output element: `Float(T)` for `abs`, else `T`.
pub fn Out(comptime T: type, comptime op: Unary) type {
    return if (op == .abs) Float(T) else T;
}

/// Element `i` of a view, through `index` for any stride's sign.
inline fn ref(v: anytype, i: usize) @TypeOf(&v.data[0]) {
    return &v.data[v.index(i)];
}

fn unit(v: anytype) bool {
    return v.n <= 1 or v.stride == 1;
}

/// `y[i] = op(x[i])`; `y` may be `x`.
pub fn unary(io: Io, comptime T: type, comptime op: Unary, x: ConstStrided(T), y: Strided(Out(T, op))) Error!void {
    try x.sameLength(y);
    const Ctx = struct {
        x: ConstStrided(T),
        y: Strided(Out(T, op)),

        fn run(c: @This(), ch: threading.Chunk) void {
            var i = ch.start;
            if (comptime vectorized(T) and !isComplex(T)) if (unit(c.x) and unit(c.y)) {
                const L = shared.lanes(T);
                while (i + L <= ch.end) : (i += L) {
                    const v: @Vector(L, T) = c.x.data[i..][0..L].*;
                    c.y.data[i..][0..L].* = unaryReal(op, v);
                }
            };
            while (i < ch.end) : (i += 1) ref(c.y, i).* = unaryOne(T, op, c.x.get(i));
        }
    };
    return run(io, T, x.n, Ctx{ .x = x, .y = y }, Ctx.run);
}

/// `z[i] = op(x[i], y[i])`; `z` may be `x` or `y`.
pub fn binary(io: Io, comptime T: type, comptime op: Binary, x: ConstStrided(T), y: ConstStrided(T), z: Strided(T)) Error!void {
    try x.sameLength(y);
    try x.sameLength(z);
    const Ctx = struct {
        x: ConstStrided(T),
        y: ConstStrided(T),
        z: Strided(T),

        fn run(c: @This(), ch: threading.Chunk) void {
            var i = ch.start;
            if (comptime vectorized(T) and !isComplex(T)) if (unit(c.x) and unit(c.y) and unit(c.z)) {
                const L = shared.lanes(T);
                while (i + L <= ch.end) : (i += L) {
                    const a: @Vector(L, T) = c.x.data[i..][0..L].*;
                    const b: @Vector(L, T) = c.y.data[i..][0..L].*;
                    c.z.data[i..][0..L].* = binaryReal(op, a, b);
                }
            };
            while (i < ch.end) : (i += 1) ref(c.z, i).* = binaryOne(T, op, c.x.get(i), c.y.get(i));
        }
    };
    return run(io, T, x.n, Ctx{ .x = x, .y = y, .z = z }, Ctx.run);
}

/// `f` over `n` elements' chunks, on `axpy`'s thresholds.
fn run(io: Io, comptime T: type, n: usize, ctx: anytype, comptime f: fn (@TypeOf(ctx), threading.Chunk) void) Io.Cancelable!void {
    if (n == 0) return;
    if (threading.serial(.axpy, T, n)) return f(ctx, .{ .index = 0, .start = 0, .end = n });
    const p = threading.plan(.axpy, T, n, n, .even, core.quantum);
    return threading.forEach(io, &p, ctx, f);
}
