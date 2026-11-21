//! Accelerate's level-1 routines (`apple_accel.zig`).

const std = @import("std");
const types = @import("../types.zig");
const threading = @import("../threading.zig");
const portable = @import("../portable.zig");
const util = @import("../util.zig");
const accel = @import("../apple_accel.zig");
const binding = @import("binding.zig");
const parts = @import("split.zig");
const quirks = @import("quirks.zig");

const Io = std.Io;
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;
const isComplex = types.isComplex;
const Float = types.Float;
const isZero = util.isZero;
const x86_64 = binding.x86_64;
const Blas = binding.Blas;
const covers = accel.covers;
const rosetta = accel.rosetta;
const use = accel.use;
const call = accel.call;
const span = accel.span;
const vector_quantum = parts.vector_quantum;
const split = parts.split;
const fold = parts.fold;
const Sum = parts.Sum;
const rounded = parts.rounded;
const fits = parts.fits;
const mag = parts.mag;
const whole = parts.whole;
const sub = parts.sub;
const units = parts.units;
const scalar = parts.scalar;
const sent = quirks.sent;
const zeroPart = quirks.zeroPart;

/// The element-wise level-1 routines `split` takes.
const Elem = enum { copy, swap, scal, rot, axpy };

/// One element-wise call: `x` and, but for `scal`, `y`, with its scalars.
fn Elementwise(comptime T: type, comptime r: Elem) type {
    return struct {
        const Self = @This();
        const name = "cblas_" ++ switch (r) {
            .rot => switch (T) {
                Complex32 => "csrot",
                Complex64 => "zdrot",
                else => portable.letter(T).? ++ "rot",
            },
            else => portable.letter(T).? ++ @tagName(r),
        };

        b: *const Blas,
        n: isize,
        x: [*]T,
        incx: isize,
        y: [*]T = undefined,
        incy: isize = 1,
        alpha: T = undefined,
        c: Float(T) = 0,
        s: Float(T) = 0,

        fn part(e: *const Self, ch: threading.Chunk) void {
            const len = units(ch);
            const x = sub(e.x, e.n, e.incx, ch);
            const y = if (r == .scal) e.y else sub(e.y, e.n, e.incy, ch);
            _ = switch (r) {
                .copy, .swap => call(name, e.b, &.{}, .{ len, x, e.incx, y, e.incy }),
                .scal => call(name, e.b, &.{}, .{ len, scalar(T, &e.alpha), x, e.incx }),
                .rot => call(name, e.b, &.{}, .{ len, x, e.incx, y, e.incy, e.c, e.s }),
                .axpy => call(name, e.b, &.{}, .{ len, scalar(T, &e.alpha), x, e.incx, y, e.incy }),
            } orelse unreachable;
        }

        /// The call over its parts; false where rstd runs it.
        fn run(e: *const Self, io: Io) Io.Cancelable!bool {
            if (!fits(name, e.b, &.{ e.n, mag(e.incx), mag(e.incy), span(e.n, e.incx), span(e.n, e.incy) })) return false;
            if (e.incy == 0 or (e.incx == 0 and (r == .swap or r == .scal or r == .rot))) return false;
            // a zero stride's elements depend on each other, as for rstd's own
            if (e.incx == 0) {
                e.part(whole(e.n));
                return true;
            }
            const len: usize = @intCast(e.n);
            try split(io, @field(threading.Family, @tagName(r)), T, len, len, vector_quantum, e, part);
            return true;
        }
    };
}

pub inline fn copy(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) Io.Cancelable!bool {
    if (comptime !covers(T)) return false else {
        if (!sent(.copy, T, incx == 1 and incy == 1)) return false;
        const b = use(T) orelse return false;
        const e: Elementwise(T, .copy) = .{ .b = b, .n = n, .x = x, .incx = incx, .y = y, .incy = incy };
        return e.run(io);
    }
}

pub inline fn swap(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) Io.Cancelable!bool {
    if (comptime !covers(T)) return false else {
        if (!sent(.swap, T, incx == 1 and incy == 1)) return false;
        const b = use(T) orelse return false;
        const e: Elementwise(T, .swap) = .{ .b = b, .n = n, .x = x, .incx = incx, .y = y, .incy = incy };
        return e.run(io);
    }
}

pub inline fn scal(io: Io, comptime T: type, n: isize, alpha: T, x: [*]T, incx: isize) Io.Cancelable!bool {
    if (comptime !covers(T)) return false else {
        if (!sent(.scal, T, incx == 1)) return false;
        const b = use(T) orelse return false;
        // arm64 skips a zero part's products; x86_64 zeroes `x` for a zero `alpha`
        if (zeroPart(T, alpha) and (!x86_64 or isZero(T, alpha))) return false;
        const e: Elementwise(T, .scal) = .{ .b = b, .n = n, .x = x, .incx = incx, .alpha = alpha };
        return e.run(io);
    }
}

/// The reductions `fold` takes: `dot`, `dotc`, `dsdot` (`f32` only) and
/// `asum` (`x` only, `incx > 0`).
const Red = enum { dot, dotc, dsdot, asum };

/// One reduction, its parts summed in `Sum(T)`.
fn Reduction(comptime T: type, comptime r: Red) type {
    return struct {
        const Self = @This();
        /// `sdot` is `dsdot`: rstd's contract sums in `f64` (readme, One
        /// contract on every core).
        const name = if (r == .asum) switch (T) {
            Complex32 => "cblas_scasum",
            Complex64 => "cblas_dzasum",
            else => "cblas_" ++ portable.letter(T).? ++ "asum",
        } else switch (T) {
            f32 => "cblas_dsdot",
            f64 => "cblas_ddot",
            else => "cblas_" ++ portable.letter(T).? ++ (if (r == .dotc) "dotc" else "dotu") ++ "_sub",
        };
        const W = if (r == .asum) Sum(Float(T)) else Sum(T);

        b: *const Blas,
        n: isize,
        x: [*]T,
        incx: isize,
        y: [*]T = undefined,
        incy: isize = 1,

        fn part(e: *const Self, ch: threading.Chunk) W {
            const len = units(ch);
            const x = sub(e.x, e.n, e.incx, ch);
            if (r == .asum) return call(name, e.b, &.{}, .{ len, x, e.incx }).?;
            const y = sub(e.y, e.n, e.incy, ch);
            if (comptime !isComplex(T)) return call(name, e.b, &.{}, .{ len, x, e.incx, y, e.incy }).?;
            var v: T = undefined;
            _ = call(name, e.b, &.{}, .{ len, x, e.incx, y, e.incy, &v }).?;
            return .init(v.re, v.im);
        }

        /// The sum over its parts; null where rstd runs it.
        fn run(e: *const Self, io: Io) Io.Cancelable!?W {
            if (!fits(name, e.b, &.{ e.n, mag(e.incx), mag(e.incy), span(e.n, e.incx), span(e.n, e.incy) })) return null;
            // a zero stride runs serially, as rstd's own
            if (e.incx == 0 or e.incy == 0) return e.part(whole(e.n));
            return try fold(io, if (r == .asum) .asum else .dot, T, W, @intCast(e.n), e, part);
        }
    };
}

/// `x . y`, or `conj(x) . y`.
pub inline fn dot(io: Io, comptime T: type, comptime conj: bool, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) Io.Cancelable!?T {
    if (comptime !covers(T)) return null else {
        if (!sent(if (conj) .dotc else .dot, T, incx == 1 and incy == 1)) return null;
        const b = use(T) orelse return null;
        if (T == f32 and rosetta()) return null;
        const e: Reduction(T, if (conj) .dotc else .dot) = .{ .b = b, .n = n, .x = x, .incx = incx, .y = y, .incy = incy };
        return rounded(T, try e.run(io) orelse return null);
    }
}

pub inline fn dsdot(io: Io, n: isize, x: [*]f32, incx: isize, y: [*]f32, incy: isize) Io.Cancelable!?f64 {
    if (!sent(.dsdot, f32, incx == 1 and incy == 1)) return null;
    const b = use(f32) orelse return null;
    if (rosetta()) return null;
    const e: Reduction(f32, .dsdot) = .{ .b = b, .n = n, .x = x, .incx = incx, .y = y, .incy = incy };
    return e.run(io);
}

/// `s` and `c` as rstd's contract has them, the squares summed in `f64` by
/// `dsdot` and the root rounded once; `d` and `z` for `incx > 0`, the
/// reference's only stride.
pub inline fn nrm2(comptime T: type, n: isize, x: [*]T, incx: isize) ?Float(T) {
    if (comptime !covers(T)) return null else {
        if (!sent(.nrm2, T, incx == 1)) return null;
        const b = use(T) orelse return null;
        if (incx <= 0) return null;
        const f: [*]Float(T) = @ptrCast(x);
        if (comptime Float(T) == f32) {
            if (rosetta()) return null;
            const sq = if (comptime !isComplex(T))
                call("cblas_dsdot", b, &.{span(n, incx)}, .{ n, f, incx, f, incx }) orelse return null
            else if (incx == 1)
                call("cblas_dsdot", b, &.{2 * n}, .{ 2 * n, f, @as(isize, 1), f, @as(isize, 1) }) orelse return null
            else blk: {
                const re = call("cblas_dsdot", b, &.{span(n, 2 * incx)}, .{ n, f, 2 * incx, f, 2 * incx }) orelse return null;
                const im = call("cblas_dsdot", b, &.{span(n, 2 * incx)}, .{ n, f + 1, 2 * incx, f + 1, 2 * incx }) orelse return null;
                break :blk re + im;
            };
            return @floatCast(@sqrt(sq));
        }
        return call(if (comptime isComplex(T)) "cblas_dznrm2" else "cblas_dnrm2", b, &.{span(n, incx)}, .{ n, x, incx });
    }
}

/// For `incx > 0`; rstd's kernels give 0 otherwise.
pub inline fn asum(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize) Io.Cancelable!?Float(T) {
    if (comptime !covers(T)) return null else {
        if (!sent(.asum, T, incx == 1)) return null;
        const b = use(T) orelse return null;
        if (incx <= 0) return null;
        const e: Reduction(T, .asum) = .{ .b = b, .n = n, .x = x, .incx = incx };
        return rounded(Float(T), try e.run(io) orelse return null);
    }
}

pub inline fn rot(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize, c: Float(T), s: Float(T)) Io.Cancelable!bool {
    if (comptime !covers(T)) return false else {
        if (!sent(.rot, T, incx == 1 and incy == 1)) return false;
        const b = use(T) orelse return false;
        const e: Elementwise(T, .rot) = .{ .b = b, .n = n, .x = x, .incx = incx, .y = y, .incy = incy, .c = c, .s = s };
        return e.run(io);
    }
}

pub inline fn rotm(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize, param: *const [5]T) bool {
    if (comptime !covers(T)) return false else {
        if (!sent(.rotm, T, incx == 1 and incy == 1)) return false;
        const b = use(T) orelse return false;
        if (incx == 0 or incy == 0) return false;
        return call("cblas_" ++ portable.letter(T).? ++ "rotm", b, &.{ span(n, incx), span(n, incy) }, .{ n, x, incx, y, incy, param }) != null;
    }
}

pub inline fn axpy(io: Io, comptime T: type, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize) Io.Cancelable!bool {
    if (comptime !covers(T)) return false else {
        if (!sent(.axpy, T, incx == 1 and incy == 1)) return false;
        const b = use(T) orelse return false;
        // arm64 only
        if (!x86_64 and zeroPart(T, alpha)) return false;
        const e: Elementwise(T, .axpy) = .{ .b = b, .n = n, .x = x, .incx = incx, .y = y, .incy = incy, .alpha = alpha };
        return e.run(io);
    }
}

pub inline fn axpby(comptime T: type, n: isize, alpha: T, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) bool {
    if (comptime !covers(T)) return false else {
        if (!sent(.axpby, T, incx == 1 and incy == 1)) return false;
        const b = use(T) orelse return false;
        if (isZero(T, alpha)) return false;
        // arm64 only
        if (!x86_64 and (zeroPart(T, alpha) or zeroPart(T, beta))) return false;
        if (incy == 0) return false;
        return call("catlas_" ++ portable.letter(T).? ++ "axpby", b, &.{ span(n, incx), span(n, incy) }, .{ n, scalar(T, &alpha), x, incx, scalar(T, &beta), y, incy }) != null;
    }
}
