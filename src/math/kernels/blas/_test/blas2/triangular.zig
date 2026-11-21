//! `blas2/triangular.zig`'s tests.

const std = @import("std");
const types = @import("../../types.zig");
const Vtable = @import("../../Vtable.zig");
const portable = @import("../../portable.zig");
const Num = @import("../../blas2/ops.zig").Num;
const triangular = @import("../../blas2/triangular.zig");

const Trans = types.Trans;
const Uplo = types.Uplo;
const Diag = types.Diag;
const isComplex = types.isComplex;
const Float = types.Float;
const canon = triangular.canon;
const Kind = triangular.Kind;
const call = triangular.call;

const check = @import("../check.zig");
const guard = @import("../guard.zig");
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;

const real = .{ f32, f64 } ++ check.x_real;
const complex = .{ Complex32, Complex64 } ++ check.x_complex;
const ops = [_]Trans{ .n, .t, .r, .c };
const sizes = [_]usize{ 1, 2, 3, 5, 8, 17, 70 };
const incs = [_]isize{ 1, 2, -1, -3 };
const max_n = 70;

const Layout = enum { full, band, packed_ };

fn layoutOf(comptime kind: Kind) Layout {
    return switch (kind) {
        .trmv, .trsv => .full,
        .tbmv, .tbsv => .band,
        .tpmv, .tpsv => .packed_,
    };
}

fn isSolve(comptime kind: Kind) bool {
    return kind == .trsv or kind == .tbsv or kind == .tpsv;
}

/// Where `A(i,j)` lives in the storage, for `i` and `j` inside the triangle.
fn index(layout: Layout, uplo: Uplo, n: usize, k: usize, lda: usize, i: usize, j: usize) usize {
    return switch (layout) {
        .full => i + j * lda,
        .band => (if (uplo == .u) k + i - j else i - j) + j * lda,
        .packed_ => if (uplo == .u) j * (j + 1) / 2 + i else j * (2 * n - j + 1) / 2 + i - j,
    };
}

fn inTriangle(layout: Layout, uplo: Uplo, k: usize, i: usize, j: usize) bool {
    if (uplo == .u) {
        if (i > j) return false;
        return layout != .band or j - i <= k;
    }
    if (i < j) return false;
    return layout != .band or i - j <= k;
}

/// `op(A)` as a dense row-major `n` by `n` matrix.
fn dense(comptime T: type, layout: Layout, uplo: Uplo, op: Trans, unit: bool, n: usize, k: usize, lda: usize, a: []const T, out: []T) void {
    const N = Num(T);
    const real_op = canon(T, op);
    for (0..n) |i| for (0..n) |j| {
        // (i,j) of op(A) is (r,c) of A
        const r = if (real_op.transposed()) j else i;
        const c = if (real_op.transposed()) i else j;
        var v: T = if (comptime isComplex(T)) .init(0, 0) else 0;
        if (inTriangle(layout, uplo, k, r, c)) {
            if (r == c and unit) {
                v = if (comptime isComplex(T)) .init(1, 0) else 1;
            } else {
                v = a[index(layout, uplo, n, k, lda, r, c)];
            }
            if (real_op.conjA()) v = N.conj(v, true);
        }
        out[i * n + j] = v;
    };
}

fn storageLen(layout: Layout, n: usize, lda: usize) usize {
    return switch (layout) {
        .full => lda * n,
        .band => lda * n,
        .packed_ => n * (n + 1) / 2,
    };
}

/// Physical index of logical element `i` of a vector of `n` at stride `inc`.
fn phys(i: usize, n: usize, inc: isize) usize {
    const a: usize = @intCast(@abs(inc));
    return if (inc > 0) i * a else (n - 1 - i) * a;
}

fn expectSame(comptime T: type, want: []const T, got: []const T, tol: Float(T)) !void {
    try check.expectCloseSlices(T, want, got, tol);
}

fn checkKind(comptime T: type, comptime kind: Kind, vt: *const Vtable) !void {
    const F = Float(T);
    const N = Num(T);
    const layout = comptime layoutOf(kind);
    const solve = comptime isSolve(kind);
    // scratch from the heap holds NaN, so a missing clear shows
    var poison: guard.Poison = .{ .child = std.testing.allocator };
    const gpa = poison.allocator();
    for (sizes) |n| for ([_]usize{ 0, 1, 3, 1000 }) |k_req| {
        if (layout != .band and k_req != 0) continue;
        const k = if (layout == .band) @min(k_req, n - 1) else 0;
        const lda = if (layout == .band) k + 1 + (n % 3) else n + (n % 3);
        for ([_]Uplo{ .u, .l }) |uplo| for (ops) |op| for (std.enums.values(Diag)) |diag| for (incs) |inc| {
            // static: at max_n a complex f80 frame would be ~330 KB
            const a = guard.operand(T, .a, max_n * (max_n + 3));
            check.fill(T, a[0..storageLen(layout, n, lda)], n * 31 + k + @backingInt(op));
            const alen = storageLen(layout, n, lda);
            // small off-diagonals and a dominant diagonal keep solves tame
            const sc: F = 1 / @sqrt(@as(F, @floatFromInt(n + 1)));
            for (a[0..alen]) |*v| v.* = if (comptime isComplex(T)) .init(v.re * sc, v.im * sc) else v.* * sc;
            for (0..n) |i| {
                const dg = &a[index(layout, uplo, n, k, lda, i, i)];
                dg.* = if (comptime isComplex(T)) .init(3 + @abs(dg.re), dg.im) else 3 + @abs(dg.*);
            }
            const d = guard.operand(T, .want, max_n * max_n);
            dense(T, layout, uplo, op, diag == .u, n, k, lda, a[0..alen], d[0 .. n * n]);

            const span = (n - 1) * @abs(inc) + 1;
            var x: [max_n * 3]T = undefined;
            check.fill(T, x[0..span], n + 17);
            var want = x;
            var xl: [max_n]T = undefined;
            for (0..n) |i| xl[i] = x[phys(i, n, inc)];
            var rl: [max_n]T = undefined;
            const zero: T = if (comptime isComplex(T)) .init(0, 0) else 0;
            if (!solve) {
                for (0..n) |i| {
                    var acc = zero;
                    for (0..n) |j| acc = N.add(acc, N.mul(d[i * n + j], xl[j], false));
                    rl[i] = acc;
                }
                for (0..n) |i| want[phys(i, n, inc)] = rl[i];
            }

            const nn: isize = @intCast(n);
            const kk: isize = @intCast(k);
            const ll: isize = @intCast(lda);
            const ii: isize = @intCast(inc);
            try call(gpa, check.io, T, kind, vt, uplo, op, diag, nn, kk, a.ptr, ll, &x, ii);
            const tol = check.tolerance(F, 4 * n + 4);
            if (!solve) {
                expectSame(T, want[0..span], x[0..span], tol) catch |e| return report(e, T, kind, uplo, op, diag, n, k, inc, vt.dtb_entries);
            } else {
                // op(A) * solution must give back the right-hand side
                for (0..n) |i| {
                    var acc = zero;
                    for (0..n) |j| acc = N.add(acc, N.mul(d[i * n + j], x[phys(j, n, inc)], false));
                    rl[i] = acc;
                }
                for (0..n) |i| {
                    check.expectClose(T, xl[i], rl[i], tol) catch |e| return report(e, T, kind, uplo, op, diag, n, k, inc, vt.dtb_entries);
                }
                // and the gaps between strided elements are untouched
                for (0..span) |p| {
                    var hit = false;
                    for (0..n) |i| hit = hit or phys(i, n, inc) == p;
                    if (!hit) try check.expectClose(T, want[p], x[p], 0);
                }
            }
        };
    };
}

fn report(e: anyerror, comptime T: type, comptime kind: Kind, uplo: Uplo, op: Trans, diag: Diag, n: usize, k: usize, inc: isize, dtb: i32) anyerror {
    check.say("{s} {s} {s} {s} {s} n={d} k={d} inc={d} dtb={d}: ", .{ @tagName(kind), @typeName(T), @tagName(uplo), @tagName(op), @tagName(diag), n, k, inc, dtb });
    return e;
}

const every_kind = std.enums.values(Kind);

/// `list`'s types through `kinds`, at the vtable's block size and a tiny one
/// that reaches the `gemv` paths.
fn Check(comptime list: anytype, comptime kinds: []const Kind) fn (*const Vtable) anyerror!void {
    return struct {
        fn run(vt: *const Vtable) anyerror!void {
            for ([_]i32{ vt.dtb_entries, 3 }) |dtb| {
                var v = vt.*;
                v.dtb_entries = dtb;
                inline for (list) |T| {
                    inline for (kinds) |kind| try checkKind(T, kind, &v);
                }
            }
        }
    }.run;
}

test "triangular routines match dense arithmetic with portable kernels" {
    var vt: Vtable = .{};
    portable.install(&vt);
    try Check(real ++ complex, every_kind)(&vt);
}

test "f80 without slots: triangular routines match dense arithmetic with portable kernels" {
    if (comptime check.f80_unslotted.len == 0) return error.SkipZigTest;
    var vt: Vtable = .{};
    portable.install(&vt);
    try Check(check.f80_unslotted, every_kind)(&vt);
}

// in three, so each fits a rig image (`build/qemu/machines.zig`'s `sparc64-debug`)
test "every linked core's kernels give the same triangular results, real" {
    try check.eachCore(Check(real, every_kind));
}

test "every linked core's kernels give the same triangular results, complex products" {
    try check.eachCore(Check(complex, &.{ .trmv, .tbmv, .tpmv }));
}

test "every linked core's kernels give the same triangular results, complex solves" {
    try check.eachCore(Check(complex, &.{ .trsv, .tbsv, .tpsv }));
}

test "triangular work space beyond the stack comes from the allocator" {
    var vt: Vtable = .{};
    portable.install(&vt);
    const n = 300;
    const gpa = std.testing.allocator;
    const a = try gpa.alloc(f64, n * n);
    defer gpa.free(a);
    const x = try gpa.alloc(f64, 2 * n);
    defer gpa.free(x);
    check.fill(f64, a, 1);
    for (a) |*v| v.* *= 0.05;
    for (0..n) |i| a[i + i * n] = 4;
    check.fill(f64, x, 2);
    const b = try gpa.dupe(f64, x);
    defer gpa.free(b);
    try call(gpa, check.io, f64, .trsv, &vt, .l, .n, .n, n, 0, a.ptr, n, x.ptr, 2);
    for (0..n) |i| {
        var acc: f64 = 0;
        for (0..i + 1) |j| acc += a[i + j * n] * x[2 * j];
        try check.expectClose(f64, b[2 * i], acc, 1e-9);
    }
    var failing = std.testing.FailingAllocator.init(gpa, .{ .fail_index = 0 });
    try std.testing.expectError(error.OutOfMemory, call(failing.allocator(), check.io, f64, .trmv, &vt, .u, .n, .u, n, 0, a.ptr, n, x.ptr, 2));
}

test "triangular argument errors, in the reference's order" {
    var vt: Vtable = .{};
    portable.install(&vt);
    const gpa = std.testing.allocator;
    var a: [16]f64 = @splat(1);
    var x: [4]f64 = @splat(1);
    const e = std.testing.expectError;
    try e(error.InvalidN, call(gpa, check.io, f64, .trmv, &vt, .u, .n, .u, -1, 0, &a, 0, &x, 0));
    try e(error.InvalidLda, call(gpa, check.io, f64, .trmv, &vt, .u, .n, .u, 4, 0, &a, 3, &x, 0));
    try e(error.InvalidIncx, call(gpa, check.io, f64, .trsv, &vt, .u, .n, .u, 4, 0, &a, 4, &x, 0));
    try e(error.InvalidK, call(gpa, check.io, f64, .tbmv, &vt, .u, .n, .u, 4, -1, &a, 4, &x, 1));
    try e(error.InvalidLda, call(gpa, check.io, f64, .tbsv, &vt, .u, .n, .u, 4, 2, &a, 2, &x, 1));
    try e(error.InvalidIncx, call(gpa, check.io, f64, .tpmv, &vt, .l, .t, .u, 4, 0, &a, 1, &x, 0));
    // n == 0 returns before touching anything, even a bad lda
    try call(gpa, check.io, f64, .trmv, &vt, .u, .n, .u, 0, 0, &a, 1, &x, 1);
    // `v` and `m` are Hermitian forms, `o` to `d` conjugate `x`: neither applies
    try e(error.InvalidUplo, call(gpa, check.io, f64, .trmv, &vt, .m, .o, .u, -1, 0, &a, 0, &x, 0));
    try e(error.InvalidTrans, call(gpa, check.io, f64, .tpsv, &vt, .l, .d, .u, -1, 0, &a, 0, &x, 0));
}
