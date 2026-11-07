//! `GEADD`: `C := alpha*op(A) + beta*C`, or the same on `C^T`. The portable
//! kernels (`*_generic`) and their tests against a plain loop.

const std = @import("std");
const blas = @import("../../__root__.zig");
const types = @import("../../types.zig");
const Vtable = @import("../../Vtable.zig");
const portable = @import("../../portable.zig");
const check = @import("../../_test/blas1/check.zig");

const BlasLong = types.BlasLong;
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;

// ---- portable kernels ------------------------------------------------------

/// `kernel/generic/geadd.c` and `zgeadd.c`: per column of `op(C)` (`rows`
/// by `cols`), `axpby` of that column of `op(A)`; for a zero `alpha`, `scal`
/// by `beta` with flag 0 instead, so a zero `beta` stores zeros. Neither reads
/// the operand its zero scalar multiplies. A transposed operand is read across
/// its rows. The vtable's `axpby` and `scal` do the work, so they follow the
/// running core, as upstream's do.
fn Geadd(comptime T: type) type {
    const F = types.Float(T);
    const cc = types.kernelCc(F);
    return struct {
        fn real(rows: BlasLong, cols: BlasLong, alpha: F, a: [*]F, lda: BlasLong, beta: F, b: [*]F, ldb: BlasLong, transa: BlasLong, transb: BlasLong) callconv(cc) i32 {
            if (rows <= 0 or cols <= 0) return 0;
            const w = Walk.of(lda, ldb, transa, transb);
            for (0..@intCast(cols)) |j| {
                const pb = b + j * w.step_b;
                if (alpha == 0) {
                    _ = portable.fetch(T, .scal_k, blas.vtable)(rows, 0, 0, beta, pb, w.inc_b, pb, 0, pb, 0);
                } else {
                    _ = portable.fetch(T, .axpby_k, blas.vtable)(rows, alpha, a + j * w.step_a, w.inc_a, beta, pb, w.inc_b);
                }
            }
            return 0;
        }

        fn complex(rows: BlasLong, cols: BlasLong, alpha_r: F, alpha_i: F, a: [*]F, lda: BlasLong, beta_r: F, beta_i: F, b: [*]F, ldb: BlasLong, transa: BlasLong, transb: BlasLong) callconv(cc) i32 {
            if (rows <= 0 or cols <= 0) return 0;
            const w = Walk.of(lda, ldb, transa, transb);
            for (0..@intCast(cols)) |j| {
                const pb = b + 2 * j * w.step_b;
                if (alpha_r == 0 and alpha_i == 0) {
                    _ = portable.fetch(T, .scal_k, blas.vtable)(rows, 0, 0, beta_r, beta_i, pb, w.inc_b, pb, 0, pb, 0);
                } else {
                    _ = portable.fetch(T, .axpby_k, blas.vtable)(rows, alpha_r, alpha_i, a + 2 * j * w.step_a, w.inc_a, beta_r, beta_i, pb, w.inc_b);
                }
            }
            return 0;
        }
    };
}

/// How the kernels walk an operand, in elements: down a column (stride 1, the
/// next column `ld` on) or, transposed, across a row (stride `ld`, the next
/// row 1 on).
const Walk = struct {
    inc_a: BlasLong,
    step_a: usize,
    inc_b: BlasLong,
    step_b: usize,

    fn of(lda: BlasLong, ldb: BlasLong, transa: BlasLong, transb: BlasLong) Walk {
        return .{
            .inc_a = if (transa != 0) lda else 1,
            .step_a = if (transa != 0) 1 else @intCast(lda),
            .inc_b = if (transb != 0) ldb else 1,
            .step_b = if (transb != 0) 1 else @intCast(ldb),
        };
    }
};

pub const sgeadd_generic = &Geadd(f32).real;
pub const dgeadd_generic = &Geadd(f64).real;
pub const cgeadd_generic = &Geadd(Complex32).complex;
pub const zgeadd_generic = &Geadd(Complex64).complex;

/// The portable kernels over `T`, any float or its `Complex`, by the slot
/// they would fill (`portable.kernel`).
pub fn Portable(comptime T: type) type {
    return struct {
        pub const geadd_k = if (types.isComplex(T)) &Geadd(T).complex else &Geadd(T).real;
    };
}

/// Points every geadd slot at the portable kernels; upstream has no `q` or
/// `x` ones.
pub fn installGeneric(vt: *Vtable) void {
    vt.sgeadd_k = sgeadd_generic;
    vt.dgeadd_k = dgeadd_generic;
    vt.cgeadd_k = cgeadd_generic;
    vt.zgeadd_k = zgeadd_generic;
}

// ---- tests -----------------------------------------------------------------

fn isZero(comptime T: type, v: T) bool {
    return if (comptime types.isComplex(T)) v.re == 0 and v.im == 0 else v == 0;
}

/// `kern` against a plain loop: every transpose pair, `alpha` and `beta`
/// zero or not, NaN where a zero scalar must not read, and padding that must
/// stay as it was.
fn checkGeadd(comptime T: type, kern: anytype) !void {
    const F = types.Float(T);
    const complex = comptime types.isComplex(T);
    const al = std.testing.allocator;
    const nan = std.math.nan(F);
    const nan_t: T = if (complex) T.init(nan, nan) else nan;
    const zero: T = if (complex) T.init(0, 0) else 0;
    const alphas = if (complex) [_]T{ T.init(1.5, -0.75), T.init(0, 0) } else [_]T{ 1.5, 0 };
    const betas = if (complex) [_]T{ T.init(-0.5, 0.25), T.init(0, 0), T.init(1, 0) } else [_]T{ -0.5, 0, 1 };
    for ([_]usize{ 1, 2, 5, 17 }) |rows| for ([_]usize{ 1, 3, 8 }) |cols| for ([_]bool{ false, true }) |ta| for ([_]bool{ false, true }) |tb| for (alphas) |alpha| for (betas) |beta| {
        // `op(A)` and `op(C)` are rows by cols; stored transposed, cols by rows
        const lda = (if (ta) cols else rows) + 2;
        const ldb = (if (tb) cols else rows) + 3;
        const a = try al.alloc(T, lda * @max(rows, cols));
        defer al.free(a);
        const b = try al.alloc(T, ldb * @max(rows, cols));
        defer al.free(b);
        check.fill(T, a, rows * 7 + cols);
        check.fill(T, b, rows + cols * 11);
        const a0 = isZero(T, alpha);
        const b0 = isZero(T, beta);
        if (a0) @memset(a, nan_t);
        if (b0) for (0..cols) |j| for (0..rows) |i| {
            b[at(tb, ldb, i, j)] = nan_t;
        };
        const want = try al.dupe(T, b);
        defer al.free(want);
        for (0..cols) |j| for (0..rows) |i| {
            const k = at(tb, ldb, i, j);
            const x = a[at(ta, lda, i, j)];
            const y = b[k];
            want[k] = if (a0 and b0) zero else if (a0) check.mulW(beta, y) else if (b0) check.mulW(alpha, x) else check.addW(check.mulW(alpha, x), check.mulW(beta, y));
        };
        const pa: [*]F = @ptrCast(a.ptr);
        const pb: [*]F = @ptrCast(b.ptr);
        const r: BlasLong = @intCast(rows);
        const c: BlasLong = @intCast(cols);
        if (complex) {
            _ = kern(r, c, alpha.re, alpha.im, pa, @intCast(lda), beta.re, beta.im, pb, @intCast(ldb), @intFromBool(ta), @intFromBool(tb));
        } else {
            _ = kern(r, c, alpha, pa, @intCast(lda), beta, pb, @intCast(ldb), @intFromBool(ta), @intFromBool(tb));
        }
        check.expectCloseSlices(T, want, b, check.tolerance(F, 4)) catch |e| {
            check.say("geadd {s} rows={d} cols={d} ta={} tb={} alpha={any} beta={any}: ", .{ @typeName(T), rows, cols, ta, tb, alpha, beta });
            return e;
        };
    };
}

/// Where `(i, j)` of `op(X)` is, `X` stored transposed where `t`.
fn at(t: bool, ld: usize, i: usize, j: usize) usize {
    return if (t) j + i * ld else i + j * ld;
}

test "portable geadd kernels match a plain loop" {
    try checkGeadd(f32, sgeadd_generic);
    try checkGeadd(f64, dgeadd_generic);
    try checkGeadd(Complex32, cgeadd_generic);
    try checkGeadd(Complex64, zgeadd_generic);
    try checkGeadd(check.Ext, Portable(check.Ext).geadd_k);
    try checkGeadd(std.math.Complex(check.Ext), Portable(std.math.Complex(check.Ext)).geadd_k);
}
