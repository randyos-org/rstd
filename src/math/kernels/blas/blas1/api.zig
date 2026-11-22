//! The public level-1 routines, on strided views (`strided.zig`, issue #309):
//! the views' counts must agree, and a negative stride is the view's vector
//! read backward wherever the reference would do nothing (`scal`, `asum`,
//! `sum`, the extremes). Each delegates to `interface.zig`, the
//! reference-shaped core the LAPACK routines and the ported tests call; a view
//! always fits memory, so the core's stride errors cannot occur.

const std = @import("std");
const types = @import("../types.zig");
const strided = @import("../strided.zig");
const core = @import("interface.zig");
const compensated = @import("compensated.zig");
const fast = @import("fast.zig");
const util = @import("../util.zig");

const Io = std.Io;
const BFloat16 = types.BFloat16;
const Float = types.Float;
const Strided = strided.Strided;
const ConstStrided = strided.ConstStrided;
const LengthError = strided.LengthError;

pub const Extreme = core.Extreme;
/// How `sumWith`, `asumWith`, `dotWith` and `dotcWith` sum; each gives the
/// same bits at any worker count. `.plain`, the default: `f32` in `f64` on
/// every core. `.compensated`: Neumaier's, about one rounding of the result,
/// at portable speed (`compensated.zig`). `.fast`: speed over accuracy for
/// single precision in `f32` lanes: `asum`, `sum` and complex `dot` on the
/// core's own kernels, OpenBLAS's where it has one, so their bits may differ
/// between cores (issue #584); real `dot` in `fast.zig`; `.plain` for wider
/// types. The four are inline, so a mode known at the call compiles to its
/// path alone, without the others' frame (issue #677).
pub const Summation = enum { plain, compensated, fast };

/// What a routine of two vectors that may thread returns.
pub const PairError = LengthError || Io.Cancelable;

fn len(v: anytype) isize {
    return @intCast(v.n);
}

/// The core's stride: any, for at most one element, where it has no meaning.
fn inc(v: anytype) isize {
    return if (v.n <= 1) 1 else v.stride;
}

/// The stride's size, for the routines whose result does not depend on order.
fn step(v: anytype) isize {
    return if (v.n <= 1) 1 else @intCast(@abs(v.stride));
}

fn ptr(comptime T: type, v: anytype) [*]T {
    return @constCast(v.data.ptr);
}

/// The core's errors past a view: only cancellation.
fn canceled(e: core.Error) Io.Cancelable {
    return switch (e) {
        error.Canceled => error.Canceled,
        error.InvalidIncx, error.InvalidIncy => unreachable,
    };
}

/// For the core's routines that never cancel.
fn never(e: core.Error) noreturn {
    switch (e) {
        error.Canceled, error.InvalidIncx, error.InvalidIncy => unreachable,
    }
}

// ---- copy, swap, scal --------------------------------------------------------

/// `y := x`.
pub fn copy(io: Io, comptime T: type, x: ConstStrided(T), y: Strided(T)) PairError!void {
    try y.sameLength(x);
    core.copy(io, T, len(x), ptr(T, x), inc(x), ptr(T, y), inc(y)) catch |e| return canceled(e);
}

/// `x <-> y`.
pub fn swap(io: Io, comptime T: type, x: Strided(T), y: Strided(T)) PairError!void {
    try x.sameLength(y);
    core.swap(io, T, len(x), ptr(T, x), inc(x), ptr(T, y), inc(y)) catch |e| return canceled(e);
}

/// `x := alpha*x`; with `alpha` 0, NaN and infinity stay NaN, as `alpha*x`
/// gives them.
pub fn scal(io: Io, comptime T: type, alpha: T, x: Strided(T)) Io.Cancelable!void {
    core.scal(io, T, len(x), alpha, ptr(T, x), step(x)) catch |e| return canceled(e);
}

/// Complex `x` scaled by a real `alpha`.
pub fn scalReal(io: Io, comptime T: type, alpha: Float(T), x: Strided(T)) Io.Cancelable!void {
    core.scalReal(io, T, len(x), alpha, ptr(T, x), step(x)) catch |e| return canceled(e);
}

// ---- dot ---------------------------------------------------------------------

/// `x . y`, 0 for empty vectors.
pub fn dot(io: Io, comptime T: type, x: ConstStrided(T), y: ConstStrided(T)) PairError!T {
    try x.sameLength(y);
    return core.dot(io, T, len(x), ptr(T, x), inc(x), ptr(T, y), inc(y)) catch |e| return canceled(e);
}

/// `conj(x) . y`. Complex only.
pub fn dotc(io: Io, comptime T: type, x: ConstStrided(T), y: ConstStrided(T)) PairError!T {
    try x.sameLength(y);
    return core.dotc(io, T, len(x), ptr(T, x), inc(x), ptr(T, y), inc(y)) catch |e| return canceled(e);
}

/// `dot` summed as `s` says.
pub inline fn dotWith(io: Io, comptime T: type, x: ConstStrided(T), y: ConstStrided(T), s: Summation) PairError!T {
    return dotsWith(io, T, false, x, y, s);
}

/// `dotc` summed as `s` says.
pub inline fn dotcWith(io: Io, comptime T: type, x: ConstStrided(T), y: ConstStrided(T), s: Summation) PairError!T {
    return dotsWith(io, T, true, x, y, s);
}

inline fn dotsWith(io: Io, comptime T: type, comptime conj: bool, x: ConstStrided(T), y: ConstStrided(T), s: Summation) PairError!T {
    const plain = s == .plain or (s == .fast and !fast.serves(T));
    if (plain) return if (conj) dotc(io, T, x, y) else dot(io, T, x, y);
    try x.sameLength(y);
    // complex: the core's own kernel, OpenBLAS's where it has one (issue #584)
    if (comptime types.isComplex(T)) if (s == .fast) {
        const f = if (conj) core.dotcFast else core.dotFast;
        return f(io, T, len(x), ptr(T, x), inc(x), ptr(T, y), inc(y)) catch |e| return canceled(e);
    };
    const px = util.first(T, ptr(T, x), len(x), inc(x));
    const py = util.first(T, ptr(T, y), len(y), inc(y));
    if (comptime fast.serves(T) and !types.isComplex(T)) if (s == .fast) return fast.dot(io, x.n, px, inc(x), py, inc(y));
    return compensated.dot(io, T, conj, x.n, px, inc(x), py, inc(y));
}

/// The dot of two `f32` vectors, summed and returned in `f64`.
pub fn dsdot(io: Io, x: ConstStrided(f32), y: ConstStrided(f32)) PairError!f64 {
    try x.sameLength(y);
    return core.dsdot(io, len(x), ptr(f32, x), inc(x), ptr(f32, y), inc(y)) catch |e| return canceled(e);
}

/// `alpha + x . y`, summed in `f64`, returned as `f32`.
pub fn sdsdot(io: Io, alpha: f32, x: ConstStrided(f32), y: ConstStrided(f32)) PairError!f32 {
    try x.sameLength(y);
    return core.sdsdot(io, len(x), alpha, ptr(f32, x), inc(x), ptr(f32, y), inc(y)) catch |e| return canceled(e);
}

// ---- norms and sums ----------------------------------------------------------

/// The Euclidean norm, 0 for an empty vector.
pub fn nrm2(io: Io, comptime T: type, x: ConstStrided(T)) Float(T) {
    return core.nrm2(io, T, len(x), ptr(T, x), inc(x)) catch |e| never(e);
}

/// The sum of `|x|`, of `|re| + |im|` for complex; summed in memory order.
pub fn asum(io: Io, comptime T: type, x: ConstStrided(T)) Io.Cancelable!Float(T) {
    return core.asum(io, T, len(x), ptr(T, x), step(x)) catch |e| return canceled(e);
}

/// The sum of `x`, of `re + im` for complex; summed in memory order.
pub fn sum(io: Io, comptime T: type, x: ConstStrided(T)) Io.Cancelable!Float(T) {
    return core.sum(io, T, len(x), ptr(T, x), step(x)) catch |e| return canceled(e);
}

/// `asum` summed as `s` says.
pub inline fn asumWith(io: Io, comptime T: type, x: ConstStrided(T), s: Summation) Io.Cancelable!Float(T) {
    return switch (s) {
        .plain => asum(io, T, x),
        .compensated => compensated.asum(io, T, x.n, ptr(T, x), step(x)),
        // the core's own kernel, OpenBLAS's where it has one (issue #584)
        .fast => core.asumFast(io, T, len(x), ptr(T, x), step(x)) catch |e| return canceled(e),
    };
}

/// `sum` summed as `s` says.
pub inline fn sumWith(io: Io, comptime T: type, x: ConstStrided(T), s: Summation) Io.Cancelable!Float(T) {
    return switch (s) {
        .plain => sum(io, T, x),
        .compensated => compensated.sum(io, T, x.n, ptr(T, x), step(x)),
        .fast => core.sumFast(io, T, len(x), ptr(T, x), step(x)) catch |e| return canceled(e),
    };
}

// ---- extremes ----------------------------------------------------------------

/// The index of the first extreme element, or of the first NaN; null for an
/// empty vector. By magnitude (`|re| + |im|` for complex) or by value
/// (`max`, `min`: real only).
pub fn index(io: Io, comptime T: type, comptime which: Extreme, x: ConstStrided(T)) Io.Cancelable!?usize {
    if (x.n == 0) return null;
    if (x.n > 1 and x.stride < 0) {
        // first in the view's order, which the kernels, forward only, would
        // find last on ties
        const B = core.Best(T, which);
        var best = B.of(x.get(0), 0);
        for (1..x.n) |i| best = B.pick(best, B.of(x.get(i), @intCast(i)));
        return @intCast(best.at);
    }
    const r = core.index(io, T, which, len(x), ptr(T, x), inc(x)) catch |e| return canceled(e);
    return @intCast(r - 1);
}

/// The extreme value itself (`|x|` for `amax` and `amin`), or the first NaN;
/// null for an empty vector.
pub fn extreme(io: Io, comptime T: type, comptime which: Extreme, x: ConstStrided(T)) Io.Cancelable!?Float(T) {
    if (x.n == 0) return null;
    return core.extreme(io, T, which, len(x), ptr(T, x), step(x)) catch |e| return canceled(e);
}

// ---- rotations ---------------------------------------------------------------

/// The plane rotation `x := c*x + s*y`, `y := c*y - s*x`, real `c` and `s`.
pub fn rot(io: Io, comptime T: type, x: Strided(T), y: Strided(T), c: Float(T), s: Float(T)) PairError!void {
    try x.sameLength(y);
    core.rot(io, T, len(x), ptr(T, x), inc(x), ptr(T, y), inc(y), c, s) catch |e| return canceled(e);
}

/// The modified rotation `param` describes, `param[0]` its flag. Real only.
pub fn rotm(io: Io, comptime T: type, x: Strided(T), y: Strided(T), param: *const [5]T) LengthError!void {
    try x.sameLength(y);
    core.rotm(io, T, len(x), ptr(T, x), inc(x), ptr(T, y), inc(y), param) catch |e| never(e);
}

// ---- axpy, axpby -------------------------------------------------------------

/// `y := alpha*x + y`.
pub fn axpy(io: Io, comptime T: type, alpha: T, x: ConstStrided(T), y: Strided(T)) PairError!void {
    try y.sameLength(x);
    core.axpy(io, T, len(x), alpha, ptr(T, x), inc(x), ptr(T, y), inc(y)) catch |e| return canceled(e);
}

/// `y := alpha*conj(x) + y`. Complex only.
pub fn axpyc(io: Io, comptime T: type, alpha: T, x: ConstStrided(T), y: Strided(T)) PairError!void {
    try y.sameLength(x);
    core.axpyc(io, T, len(x), alpha, ptr(T, x), inc(x), ptr(T, y), inc(y)) catch |e| return canceled(e);
}

/// `y := alpha*x + beta*y`.
pub fn axpby(io: Io, comptime T: type, alpha: T, x: ConstStrided(T), beta: T, y: Strided(T)) LengthError!void {
    try y.sameLength(x);
    core.axpby(io, T, len(x), alpha, ptr(T, x), inc(x), beta, ptr(T, y), inc(y)) catch |e| never(e);
}

// ---- bfloat16 ----------------------------------------------------------------

/// The dot of two bfloat16 vectors in `f32`, 0 for empty ones.
pub fn sbdot(io: Io, x: ConstStrided(BFloat16), y: ConstStrided(BFloat16)) LengthError!f32 {
    try x.sameLength(y);
    return core.sbdot(io, len(x), ptr(BFloat16, x), inc(x), ptr(BFloat16, y), inc(y)) catch |e| never(e);
}

/// `x := alpha*x` over bfloat16.
pub fn bscal(io: Io, alpha: BFloat16, x: Strided(BFloat16)) void {
    core.bscal(io, len(x), alpha, ptr(BFloat16, x), step(x)) catch |e| never(e);
}

/// `out := bf16(in)`, rounding to nearest even; `f32` or `f64` in.
pub fn sbtobf16(io: Io, comptime F: type, in: ConstStrided(F), out: Strided(BFloat16)) LengthError!void {
    try out.sameLength(in);
    core.sbtobf16(io, F, len(in), ptr(F, in), inc(in), ptr(BFloat16, out), inc(out)) catch |e| never(e);
}

/// `out := in` widened to `f32` or `f64`.
pub fn bf16to(io: Io, comptime F: type, in: ConstStrided(BFloat16), out: Strided(F)) LengthError!void {
    try out.sameLength(in);
    core.bf16to(io, F, len(in), ptr(BFloat16, in), inc(in), ptr(F, out), inc(out)) catch |e| never(e);
}

pub const rotg = core.rotg;
pub const zrotg = core.zrotg;
pub const rotmg = core.rotmg;
pub const cabs1 = core.cabs1;

// ---- tests -------------------------------------------------------------------

const check = @import("../_test/check.zig");

/// `v`'s elements in its order.
fn gather(comptime T: type, v: anytype, out: []T) []T {
    for (out[0..v.n], 0..) |*o, i| o.* = v.get(i);
    return out[0..v.n];
}

test "api: a negative stride is the vector read backward" {
    inline for (.{ f32, f64, types.Complex64 }) |T| {
        var a: [9]T = undefined;
        var b: [9]T = undefined;
        check.fill(T, &a, 1);
        check.fill(T, &b, 2);
        const x = try ConstStrided(T).init(&a, 3, -4);
        const y = try Strided(T).init(&b, 3, 3);
        var buf: [3]T = undefined;
        var rx: [3]T = undefined;
        @memcpy(&rx, gather(T, x, &buf));
        // against contiguous copies in the views' order
        var cy: [3]T = undefined;
        @memcpy(&cy, gather(T, y.toConst(), &buf));
        const d = try dot(check.io, T, x, y.toConst());
        // strided and contiguous take different kernel paths: close, not the same bits
        try check.expectClose(T, try dot(check.io, T, .of(&rx), .of(&cy)), d, check.tolerance(Float(T), 2));
        try copy(check.io, T, x, y);
        for (0..3) |i| try check.expectClose(T, rx[i], y.get(i), 0);
        try check.expectClose(Float(T), nrm2(check.io, T, .of(&rx)), nrm2(check.io, T, x), check.tolerance(Float(T), 2));
    }
}

test "api: scal, asum, sum and the extremes take every element of a backward view" {
    var a = [_]f64{ 1, 9, -2, 9, 3, 0, -4 };
    const x = try Strided(f64).init(&a, 4, -2);
    try check.expectClose(f64, 10, try asum(check.io, f64, x.toConst()), 0);
    try check.expectClose(f64, -2, try sum(check.io, f64, x.toConst()), 0);
    try check.expectClose(f64, 4, (try extreme(check.io, f64, .largest_magnitude, x.toConst())).?, 0);
    try scal(check.io, f64, 2, x);
    try std.testing.expectEqualSlices(f64, &.{ 2, 9, -4, 9, 6, 0, -8 }, &a);
}

test "api: index is the view's first extreme, backward too" {
    const a = [_]f64{ 5, 0, -5, 0, 1 };
    // view order: 1, -5, 5; the first |5| in it is -5, at 1
    const back = try ConstStrided(f64).init(&a, 3, -2);
    try std.testing.expectEqual(@as(?usize, 1), try index(check.io, f64, .largest_magnitude, back));
    try std.testing.expectEqual(@as(?usize, 0), try index(check.io, f64, .largest_magnitude, try ConstStrided(f64).init(&a, 3, 2)));
    try std.testing.expectEqual(@as(?usize, 2), try index(check.io, f64, .largest, back));
    const nan = [_]f64{ std.math.nan(f64), 1, std.math.nan(f64) };
    try std.testing.expectEqual(@as(?usize, 0), try index(check.io, f64, .largest_magnitude, try ConstStrided(f64).init(&nan, 3, -1)));
}

test "api: empty views, and counts that disagree" {
    var a = [_]f64{ 1, 2, 3 };
    const e = try ConstStrided(f64).init(&a, 0, 1);
    try std.testing.expectEqual(@as(?usize, null), try index(check.io, f64, .largest_magnitude, e));
    try std.testing.expectEqual(@as(?f64, null), try extreme(check.io, f64, .largest, e));
    try check.expectClose(f64, 0, nrm2(check.io, f64, e), 0);
    try check.expectClose(f64, 0, try dot(check.io, f64, e, e), 0);
    const x = ConstStrided(f64).of(&a);
    const y = try Strided(f64).init(&a, 2, 1);
    try std.testing.expectError(error.LengthMismatch, dot(check.io, f64, x, y.toConst()));
    try std.testing.expectError(error.LengthMismatch, axpy(check.io, f64, 1, x, y));
    try std.testing.expectError(error.LengthMismatch, rotm(check.io, f64, try Strided(f64).init(&a, 3, 1), y, &.{ -2, 0, 0, 0, 0 }));
}
