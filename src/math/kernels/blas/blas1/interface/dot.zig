//! The level-1 interface's dot products: `dot`, `dotc` and their `.fast` forms, `dsdot`, `sdsdot`.

const std = @import("std");
const blas = @import("../../__root__.zig");
const threading = @import("../../threading.zig");
const types = @import("../../types.zig");
const util = @import("../../util.zig");
const portable = @import("../../portable.zig");
const accel = @import("../../apple_accel.zig");
const iroot = @import("../interface.zig");

const Io = std.Io;
const isComplex = types.isComplex;
const Float = types.Float;
const raw = util.raw;
const Error = iroot.Error;
const Pair = iroot.Pair;
const Wide = iroot.Wide;
const chunked = iroot.chunked;
const dsdotPart = iroot.dsdotPart;
const first = iroot.first;
const fold = iroot.fold;
const plus = iroot.plus;
const rounded = iroot.rounded;
const shift = iroot.shift;
const spanXY = iroot.spanXY;
const sums = iroot.sums;
const zero = iroot.zero;

// ---- dot ---------------------------------------------------------------------

/// Which slots a single-precision `asum`, `sum` or complex `dot` takes
/// (`api.Summation`, issue #584): `.plain`'s `*_wide_k`, summed in `f64` on
/// every core, or `.fast`'s, the core's own (OpenBLAS's where it has one) in
/// `f32` lanes. Other precisions have one slot.
pub const Lanes = enum { plain, fast };

/// Whether `T` is single precision, whose `asum`, `sum` and complex `dot`
/// never go to Accelerate: `.plain` sums in rstd's `f64` slots, as
/// Accelerate's precision there is untested, and `.fast` is the core's own
/// kernel, which on the Intel Mac beats Accelerate's complex `dot` 4 to 1
/// (issue #584).
pub fn single(comptime T: type) bool {
    return Float(T) == f32;
}

/// Whether `T` takes `.plain`'s own slot under `lanes`.
pub fn wide(comptime T: type, comptime lanes: Lanes) bool {
    return lanes == .plain and single(T);
}

/// `?DOT`/`?DOTU` (`interface/dot.c`, `zdot.c`): `x . y`, 0 for `n <= 0`;
/// complex `f32` summed in `f64`.
pub inline fn dot(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) Error!T {
    return dotImpl(io, T, false, .plain, n, x, incx, y, incy);
}

/// `?DOTC`: `conj(x) . y`. Complex only.
pub inline fn dotc(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) Error!T {
    comptime util.requireComplex(T, "dotc");
    return dotImpl(io, T, true, .plain, n, x, incx, y, incy);
}

/// `dot` for `Summation.fast`: complex `f32` on the core's own kernel.
pub fn dotFast(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) Error!T {
    return dotImpl(io, T, false, .fast, n, x, incx, y, incy);
}

/// `dotc` for `Summation.fast`.
pub fn dotcFast(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) Error!T {
    comptime util.requireComplex(T, "dotcFast");
    return dotImpl(io, T, true, .fast, n, x, incx, y, incy);
}

inline fn dotImpl(io: Io, comptime T: type, comptime conj: bool, comptime lanes: Lanes, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) Error!T {
    if (n <= 0) return zero(T);
    try spanXY(T, T, n, incx, incy);
    const complex_op: portable.Op = if (comptime wide(T, lanes)) (if (conj) .dotc_wide_k else .dotu_wide_k) else (if (conj) .dotc_k else .dotu_k);
    if (comptime !(single(T) and isComplex(T))) if (try accel.dot(io, T, conj, n, x, incx, y, incy)) |r| return r;
    // `kernel/x86_64/ddot.c`, `zdot.c`: a zero stride runs serially; the
    // serial test first, as `summed`'s (issue #677)
    if (comptime sums(T)) if (incx != 0 and incy != 0) if (!threading.serial(.dot, T, @intCast(n))) if (try dotChunks(io, T, complex_op, n, x, incx, y, incy)) |r| return r;
    const px = raw(T, first(T, x, n, incx));
    const py = raw(T, first(T, y, n, incy));
    if (comptime isComplex(T)) {
        return types.complexOf(portable.fetch(T, complex_op, blas.vtable)(n, px, incx, py, incy));
    }
    return portable.fetch(T, .dot_k, blas.vtable)(n, px, incx, py, incy);
}

/// `dotImpl`'s threaded path, null for one chunk.
noinline fn dotChunks(io: Io, comptime T: type, comptime complex_op: portable.Op, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) Io.Cancelable!?T {
    const p = chunked(.dot, T, n) orelse return null;
    const c: Pair(T, void) = .{ .x = first(T, x, n, incx), .incx = incx, .y = first(T, y, n, incy), .incy = incy, .e = {} };
    const S = struct {
        fn part(q: @TypeOf(c), len: isize, at: isize) Wide(T) {
            // `sdot` sums in `f64` and rounds once (`readme.md`)
            if (comptime T == f32) return dsdotPart(q, len, at);
            const a = raw(T, shift(T, q.x, at, q.incx));
            const b = raw(T, shift(T, q.y, at, q.incy));
            if (comptime isComplex(T)) {
                const r = types.complexOf(portable.fetch(T, complex_op, blas.vtable)(len, a, q.incx, b, q.incy));
                return .init(r.re, r.im);
            }
            return portable.fetch(T, .dot_k, blas.vtable)(len, a, q.incx, b, q.incy);
        }
    };
    return rounded(T, try fold(io, &p, Wide(T), c, S.part, plus(Wide(T))));
}

/// `DSDOT` (`interface/dsdot.c`): the dot of two `f32` vectors, summed and
/// returned in `f64`.
pub fn dsdot(io: Io, n: isize, x: [*]f32, incx: isize, y: [*]f32, incy: isize) Error!f64 {
    if (n <= 0) return 0;
    try spanXY(f32, f32, n, incx, incy);
    if (try accel.dsdot(io, n, x, incx, y, incy)) |r| return r;
    const c: Pair(f32, void) = .{ .x = first(f32, x, n, incx), .incx = incx, .y = first(f32, y, n, incy), .incy = incy, .e = {} };
    if (incx != 0 and incy != 0) if (chunked(.dot, f32, n)) |p| return fold(io, &p, f64, c, dsdotPart, plus(f64));
    return dsdotPart(c, n, 0);
}

/// `SDSDOT` (`interface/sdsdot.c`): `alpha + x . y`, summed in `f64` and
/// returned as `f32`; `alpha` alone for `n <= 0`.
pub fn sdsdot(io: Io, n: isize, alpha: f32, x: [*]f32, incx: isize, y: [*]f32, incy: isize) Error!f32 {
    if (n <= 0) return alpha;
    return @floatCast(try dsdot(io, n, x, incx, y, incy) + alpha);
}
