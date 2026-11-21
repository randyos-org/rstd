//! The level-1 interface's extremes (`interface.zig`): `index` and `extreme`.

const std = @import("std");
const blas = @import("../__root__.zig");
const types = @import("../types.zig");
const util = @import("../util.zig");
const portable = @import("../portable.zig");
const interface = @import("interface.zig");

const Io = std.Io;
const BlasLong = types.BlasLong;
const isComplex = types.isComplex;
const Float = types.Float;
const raw = util.raw;
const Error = interface.Error;
const spanX = interface.spanX;
const chunked = interface.chunked;
const fold = interface.fold;
const shift = interface.shift;
const Pair = interface.Pair;

/// Which extreme: by magnitude (`|re| + |im|` for complex) or by value.
pub const Extreme = enum {
    amax,
    amin,
    max,
    min,

    /// The slot of the kernel giving the extreme's value.
    pub fn op(e: Extreme) portable.Op {
        return switch (e) {
            .amax => .amax_k,
            .amin => .amin_k,
            .max => .max_k,
            .min => .min_k,
        };
    }

    /// The slot of the kernel giving its index.
    pub fn indexOp(e: Extreme) portable.Op {
        return switch (e) {
            .amax => .iamax_k,
            .amin => .iamin_k,
            .max => .imax_k,
            .min => .imin_k,
        };
    }
};

/// `I?AMAX`, `I?AMIN`, `I?MAX`, `I?MIN` (`interface/imax.c`): the 1-based
/// index of the first extreme element; 0 for `n <= 0` (and, from the
/// kernel, for `incx <= 0`). `max` and `min` are real only.
pub fn index(io: Io, comptime T: type, comptime which: Extreme, n: isize, x: [*]T, incx: isize) Error!isize {
    comptime realOnly(T, which);
    if (n <= 0) return 0;
    try spanX(T, n, incx);
    const B = Best(T, which);
    if (incx > 0) if (chunked(.index, T, n)) |p| {
        const c: Pair(T, void) = .{ .x = x, .incx = incx, .y = x, .incy = incx, .e = {} };
        const S = struct {
            fn part(q: @TypeOf(c), len: isize, at: isize) B.Pick {
                const ret = portable.fetch(T, comptime which.indexOp(), blas.vtable)(len, raw(T, shift(T, q.x, at, q.incx)), q.incx);
                const i = at + @as(isize, @intCast(ret)) - 1;
                return B.of(shift(T, q.x, i, q.incx)[0], i);
            }
        };
        return (try fold(io, &p, B.Pick, c, S.part, B.pick)).at + 1;
    };
    const ret: BlasLong = portable.fetch(T, comptime which.indexOp(), blas.vtable)(n, raw(T, x), incx);
    return @intCast(@min(ret, n));
}

/// `?AMAX`, `?AMIN`, `?MAX`, `?MIN` (`interface/max.c`, OpenBLAS extensions):
/// the extreme value itself, `|x|` for `amax`/`amin`. 0 for `n <= 0`; for
/// `incx == 0`, `x[0]` (by magnitude for `amax`/`amin`).
pub fn extreme(io: Io, comptime T: type, comptime which: Extreme, n: isize, x: [*]T, incx: isize) Error!Float(T) {
    comptime realOnly(T, which);
    if (n <= 0) return 0;
    try spanX(T, n, incx);
    if (incx == 0) {
        if (comptime isComplex(T)) return util.abs(x[0].re) + util.abs(x[0].im);
        return if (which == .amax or which == .amin) @abs(x[0]) else x[0];
    }
    const c: Pair(T, void) = .{ .x = x, .incx = incx, .y = x, .incy = incx, .e = {} };
    const S = struct {
        fn part(q: @TypeOf(c), len: isize, at: isize) Float(T) {
            return portable.fetch(T, comptime which.op(), blas.vtable)(len, raw(T, shift(T, q.x, at, q.incx)), q.incx);
        }
    };
    if (incx > 0) if (chunked(.index, T, n)) |p| return fold(io, &p, Float(T), c, S.part, Best(T, which).value);
    return S.part(c, n, 0);
}

/// The extremes' fold, exact: the first NaN, else the first extreme, as one
/// kernel call over all of `x` finds them (`extremum.zig`). Where complex
/// keys overflow, a chunk's kernel has picked by `|re|/4 + |im|/4`, and so
/// does the fold between two such chunks.
pub fn Best(comptime T: type, comptime which: Extreme) type {
    const F = Float(T);
    return struct {
        /// An element's index from 0, its key and its quarter key.
        pub const Pick = struct { at: isize, key: F, q: F };

        fn better(k: F, best: F) bool {
            return switch (which) {
                .amax, .max => k > best,
                .amin, .min => k < best,
            };
        }

        fn value(a: F, b: F) F {
            if (a != a) return a;
            if (b != b) return b;
            return if (better(b, a)) b else a;
        }

        pub fn of(v: T, at: isize) Pick {
            if (comptime isComplex(T)) return .{ .at = at, .key = util.abs(v.re) + util.abs(v.im), .q = util.abs(v.re) / 4 + util.abs(v.im) / 4 };
            const k = if (which == .amax or which == .amin) @abs(v) else v;
            return .{ .at = at, .key = k, .q = k };
        }

        pub fn pick(a: Pick, b: Pick) Pick {
            if (a.key != a.key) return a;
            if (b.key != b.key) return b;
            if (std.math.isInf(a.key) and std.math.isInf(b.key)) return if (better(b.q, a.q)) b else a;
            return if (better(b.key, a.key)) b else a;
        }
    };
}

fn realOnly(comptime T: type, comptime which: Extreme) void {
    if (isComplex(T) and (which == .max or which == .min)) @compileError(@tagName(which) ++ ": " ++ @typeName(T) ++ " has no ordering; use ." ++ (if (which == .max) "amax" else "amin") ++ ", which compares |re| + |im|");
}
