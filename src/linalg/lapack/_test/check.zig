//! Test support shared by the LAPACK tests.

const std = @import("std");
const types = @import("../../blas/types.zig");
const portable = @import("../../blas/portable.zig");
const test_io = @import("../../blas/_test/io.zig");
const Vtable = @import("../../blas/Vtable.zig");
const check = @import("../../blas/_test/check.zig");
const guard = @import("../../blas/_test/guard.zig");
const reference = @import("../../blas/_test/reference.zig");
const ops = @import("../ops.zig");

const testing = std.testing;
const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const Uplo = types.Uplo;
const re = ops.re;

/// Every element type the tests take: the four precisions, `XDouble` where
/// it is its own, and `f80` where it has no slots.
pub const test_types = .{ f32, f64, types.Complex32, types.Complex64 } ++ check.x_real ++ check.x_complex ++ check.f80_unslotted;

pub fn W(comptime T: type) type {
    return reference.Wide(T);
}

pub fn fromF64(comptime T: type, v: f64) T {
    return if (comptime isComplex(T)) .init(@floatCast(v), 0) else @floatCast(v);
}

pub fn toT(comptime T: type, v: W(T)) T {
    return if (comptime isComplex(T)) .init(@floatCast(v.re), @floatCast(v.im)) else @floatCast(v);
}

/// NaN in every part.
pub fn poisoned(comptime T: type) T {
    const x = std.math.nan(Float(T));
    return if (comptime isComplex(T)) .init(x, x) else x;
}

/// Equal bits, part by part (`f80`'s padding aside).
pub fn sameBits(comptime T: type, x: T, y: T) bool {
    const F = Float(T);
    const U = @Int(.unsigned, @bitSizeOf(F));
    if (comptime isComplex(T)) return sameBits(F, x.re, y.re) and sameBits(F, x.im, y.im);
    return @as(U, @bitCast(x)) == @as(U, @bitCast(y));
}

pub fn isNanT(comptime T: type, v: T) bool {
    return if (comptime isComplex(T)) std.math.isNan(v.re) and std.math.isNan(v.im) else std.math.isNan(v);
}

/// The worker counts the threaded tests take past one.
pub const test_workers = [_]u32{ 2, 3, 8 };

/// `run(ctx, io, want, out)` at one worker on `Io.failing` against
/// `test_workers`, on `Io.failing` (tasks inline) and on a threaded `Io`,
/// bit for bit over the `len` elements it writes.
pub fn expectSameBits(comptime T: type, ctx: anytype, comptime run: fn (@TypeOf(ctx), Io, u32, []T) anyerror!void, len: usize, what: []const u8) !void {
    const al = testing.allocator;
    const want = try al.alloc(T, len);
    defer al.free(want);
    const got = try al.alloc(T, len);
    defer al.free(got);
    try run(ctx, Io.failing, 1, want);
    var th: if (test_io.threaded_tests) test_io.TestIo else void = if (comptime test_io.threaded_tests) .init(al, .{}) else {};
    defer if (comptime test_io.threaded_tests) th.deinit();
    for (test_workers) |w| for (0..if (test_io.threaded_tests) 2 else 1) |k| {
        const io = if (comptime test_io.threaded_tests) (if (k == 0) Io.failing else th.io()) else Io.failing;
        try run(ctx, io, w, got);
        for (want, got, 0..) |x, y, i| if (!sameBits(T, x, y)) {
            check.say("{s} {s}: {d} workers ({s}) differ from one at element {d}: {any} against {any}\n", .{ what, @typeName(T), w, if (k == 0) "inline" else "threaded", i, y, x });
            return error.TestUnexpectedResult;
        };
    };
}

/// `f(ctx, io, out)` on `Io.failing` (tasks inline) and on `Io.Threaded`,
/// bit for bit; the BLAS calls take their workers by size.
pub fn expectIoBits(comptime T: type, ctx: anytype, comptime f: fn (@TypeOf(ctx), Io, []T) anyerror!void, len: usize, what: []const u8) !void {
    if (comptime !test_io.hosted) return;
    const al = testing.allocator;
    const want = try al.alloc(T, len);
    defer al.free(want);
    const got = try al.alloc(T, len);
    defer al.free(got);
    try f(ctx, Io.failing, want);
    var th: Io.Threaded = .init(al, .{});
    defer th.deinit();
    try f(ctx, th.io(), got);
    for (want, got, 0..) |x, y, i| if (!sameBits(T, x, y)) {
        check.say("{s} {s}: threaded differs from inline at element {d}\n", .{ what, @typeName(T), i });
        return error.TestUnexpectedResult;
    };
}

/// The order of the inline-against-threaded cases.
pub const big_n = 300;

pub fn uploOf(upper: bool) Uplo {
    return if (upper) .u else .l;
}

/// `vt` with every precision's `q` and `dtb_entries` shrunk, so `potrf`'s
/// base case ends at 2 columns and its blocks are 3 wide past 12; static, as
/// a vtable is large for a rig's stack.
var shrunk: Vtable = undefined;

pub fn shrink(vt: *const Vtable) *const Vtable {
    shrunk = vt.*;
    shrunk.dtb_entries = 4;
    inline for ([_]portable.Gemm{ .s, .d, .q, .c, .z, .x }) |g| portable.paramPtr(&shrunk, g, .q).* = 3;
    return &shrunk;
}

/// `f(T, v, vi)` for each `T` of `Ts` on `vt` shrunk by `shrunk_by` (`vi`
/// 0), then on `vt` itself (1), each installed for the call.
pub fn eachBlocking(vt: *const Vtable, comptime shrunk_by: fn (*const Vtable) *const Vtable, comptime Ts: anytype, comptime f: anytype) anyerror!void {
    inline for (Ts) |T| for ([_]*const Vtable{ shrunk_by(vt), vt }, 0..) |v, vi| {
        guard.use(v);
        defer guard.restore();
        try f(T, v, vi);
    };
}

/// `got` (`n x nrhs` at `ldb`), solved from `want`'s `B` with a NaN put at
/// `at`: column `at[1]` all NaN, the others `want`'s bits.
pub fn expectNanColumn(comptime T: type, what: []const u8, variant: []const u8, n: usize, nrhs: usize, ldb: usize, at: [2]usize, got: []const T, want: []const T) !void {
    for (0..nrhs) |j| for (0..n) |i| {
        const x = got[i + j * ldb];
        const ok = if (j == at[1]) isNanT(T, x) else sameBits(T, x, want[i + j * ldb]);
        if (!ok) {
            check.say("{s} {s} {s} n={d}: NaN at B({d},{d}), X({d},{d}) is {any}\n", .{ what, @typeName(T), variant, n, at[0], at[1], i, j, x });
            return error.TestUnexpectedResult;
        }
    };
}

/// `THRESH` of the reference's linear-equation tests (`TESTING/dtest.in`).
pub const thresh = 30;

/// `T` from parts.
pub fn cplx(comptime T: type, r: f64, i: f64) T {
    return if (comptime isComplex(T)) .init(@floatCast(r), @floatCast(i)) else @floatCast(r);
}

/// `DGET06` (`TESTING/LIN/dget06.f`): how far an estimate `rcond` is from
/// the true `rcondc`, as a ratio less `1 - eps`.
pub fn ratio(comptime F: type, rcond: F, rcondc: F) F {
    const eps = std.math.floatEps(F);
    if (rcond > 0) return if (rcondc > 0) @max(rcond, rcondc) / @min(rcond, rcondc) - (1 - eps) else rcond / eps;
    return if (rcondc > 0) rcondc / eps else 0;
}

// ---- generators and checks -------------------------------------------------------

/// `a := X Y` of rank `r`, `X` (`a.rows x r`) and `Y` (`r x a.cols`)
/// random from `seed` and `seed + 1`, each sum in order.
pub fn rankProduct(comptime T: type, al: Allocator, a: ops.Mat(T), r: usize, seed: u64) Allocator.Error!void {
    const m = a.rows;
    const n = a.cols;
    const x = try al.alloc(T, m * r);
    defer al.free(x);
    const y = try al.alloc(T, r * n);
    defer al.free(y);
    check.fill(T, x, seed);
    check.fill(T, y, seed + 1);
    for (0..n) |j| for (0..m) |i| {
        var s = ops.zero(T);
        for (0..r) |p| s = ops.add(T, s, ops.mul(T, x[i + p * m], y[p + j * r]));
        a.at(i, j).* = s;
    };
}

/// A random Hermitian `A + Aᴴ`.
pub fn hermitian(comptime T: type, a: ops.Mat(T), seed: u64) void {
    const n = a.rows;
    check.fill(T, a.ptr[0 .. a.ld * n], seed);
    for (0..n) |j| for (0..j + 1) |i| {
        const s = ops.add(T, a.get(i, j), ops.conj(T, a.get(j, i)));
        a.at(i, j).* = s;
        a.at(j, i).* = ops.conj(T, s);
    };
}

/// `x` within `tol` of `y`, element by element (`|x - y|`).
pub fn expectNear(comptime T: type, what: []const u8, x: ops.Mat(T), y: ops.Mat(T), tol: Float(T)) !void {
    for (0..x.cols) |j| for (0..x.rows) |i| {
        const d = ops.abs(T, ops.sub(T, x.get(i, j), y.get(i, j)));
        if (!(d <= tol)) {
            check.say("{s} {s} {d} x {d}: ({d}, {d}) off by {e}, tolerance {e}\n", .{ what, @typeName(T), x.rows, x.cols, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(tol)) });
            return error.TestUnexpectedResult;
        }
    };
}
