//! `xTGSEN` without condition estimates (`IJOB = 0`): reordering a
//! generalized Schur form by a selection, Reference-LAPACK 3.12.0
//! (`lapack-netlib/SRC/dtgsen.f`, `ztgsen.f`, OpenBLAS v0.3.34-317;
//! modified BSD), over `tgexc.zig`'s moves.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const ops = @import("ops.zig");
const schur = @import("schur.zig");
const tgexc_mod = @import("tgexc.zig");
const qz_2x2 = @import("qz_2x2.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const isComplex = types.isComplex;
const Float = types.Float;
const ReorderVectors = schur.ReorderVectors;
const Eigenvalues = schur.Eigenvalues;
const eigenvaluesHold = schur.eigenvaluesHold;
const Pencil = tgexc_mod.Pencil;

/// What `tgsen` returns: a slice too short or a leading dimension too small,
/// `SwapRejected` for the reference's `INFO = 1`: a swap failed the
/// stability tests, `(A, B)` part reordered, the eigenvalues written.
pub const TgsenError = tgexc_mod.TgexcError;

/// `?TGSEN` with `IJOB = 0` (`dtgsen.f`, `ztgsen.f`): reorders the `n x n`
/// generalized Schur form `(A, B)` so the eigenvalues `select` picks (`n`;
/// real: either of a pair picks both) lead it, `(A, B) := Qₘᴴ (A, B) Zₘ`,
/// `Q := Q Qₘ` and `Z := Z Zₘ` as asked; returns their count. Writes the
/// eigenvalues `alpha / beta`: real, `B`'s 1 x 1 diagonal made nonnegative
/// (its row and `Q`'s column negated) and a pair's by `xLAG2`; complex,
/// `B`'s diagonal made real and nonnegative. Complex, nothing moved or a
/// swap rejected: `alpha` and `beta` the diagonals as they are (the
/// reference leaves them unset).
pub fn tgsen(allocator: Allocator, io: Io, comptime T: type, compq: ReorderVectors, compz: ReorderVectors, select: []const bool, n: usize, a: []T, lda: usize, b: []T, ldb: usize, alpha: Eigenvalues(T), beta: []T, q: []T, ldq: usize, z: []T, ldz: usize) TgsenError!usize {
    comptime ops.requireFloat(T, "tgsen");
    try util.holds(a.len, n, n, lda);
    try util.holds(b.len, n, n, ldb);
    if (compq == .update) try util.holds(q.len, n, n, ldq);
    if (compz == .update) try util.holds(z.len, n, n, ldz);
    if (select.len < n or beta.len < n) return error.SliceTooShort;
    try eigenvaluesHold(T, alpha, n);
    const p: Pencil(T) = .{ .n = n, .a = a.ptr, .lda = lda, .b = b.ptr, .ldb = ldb, .q = if (compq == .update) q.ptr else null, .ldq = ldq, .z = if (compz == .update) z.ptr else null, .ldz = ldz };
    // the cluster's size, a real pair counted whole
    var m: usize = 0;
    var k: usize = 0;
    while (k < n) : (k += 1) {
        if (opensPair(T, p, k)) {
            if (select[k] or select[k + 1]) m += 2;
            k += 1;
        } else if (select[k]) m += 1;
    }
    var rejected = false;
    if (m != n and m != 0) rejected = try collect(allocator, io, T, p, select);
    if (comptime isComplex(T)) {
        if (m == n or m == 0 or rejected) {
            for (0..n) |i| {
                alpha[i] = a[i + i * lda];
                beta[i] = b[i + i * ldb];
            }
        } else normalizeComplex(T, p, alpha, beta);
    } else eigenvaluesReal(T, p, alpha, beta);
    if (rejected) return error.SwapRejected;
    return m;
}

/// A real `A`'s row `k` opens a 2 x 2 block.
fn opensPair(comptime T: type, p: Pencil(T), k: usize) bool {
    if (comptime isComplex(T)) return false;
    return k + 1 < p.n and p.a[k + 1 + k * p.lda] != 0;
}

/// `select`'s blocks to the top-left, in order; true when a swap failed.
fn collect(allocator: Allocator, io: Io, comptime T: type, p: Pencil(T), select: []const bool) ops.Error!bool {
    var ks: usize = 0;
    var k: usize = 0;
    while (k < p.n) : (k += 1) {
        const pair = opensPair(T, p, k);
        const swap = select[k] or (pair and select[k + 1]);
        if (swap) {
            if (k != ks) {
                if (comptime isComplex(T)) {
                    if (tgexc_mod.moveComplex(T, p, k, ks) != null) return true;
                } else {
                    // `KS` is `DTGEXC`'s `ILST`, as it leaves it
                    var f = k + 1;
                    var l = ks + 1;
                    if (!try tgexc_mod.moveReal(allocator, io, T, p, &f, &l)) return true;
                    ks = l - 1;
                }
            }
            ks += if (pair) 2 else 1;
        }
        if (pair) k += 1;
    }
    return false;
}

/// `dtgsen.f`'s eigenvalues: a 1 x 1 block with `B(k, k) < 0` negated
/// (row `k` of `A` and `B`, column `k` of `Q`), a pair's from `xLAG2`.
fn eigenvaluesReal(comptime F: type, p: Pencil(F), alpha: Eigenvalues(F), beta: []F) void {
    const n = p.n;
    var k: usize = 0;
    while (k < n) : (k += 1) {
        const a = p.a;
        const b = p.b;
        if (opensPair(F, p, k)) {
            const at = struct {
                fn f(x: [*]F, ld: usize, i: usize, j: usize) F {
                    return x[i + j * ld];
                }
            }.f;
            const e = qz_2x2.eig2x2(F, .{ at(a, p.lda, k, k), at(a, p.lda, k, k + 1), at(a, p.lda, k + 1, k), at(a, p.lda, k + 1, k + 1) }, .{ at(b, p.ldb, k, k), at(b, p.ldb, k, k + 1), at(b, p.ldb, k + 1, k), at(b, p.ldb, k + 1, k + 1) });
            beta[k] = e.scale1;
            beta[k + 1] = e.scale2;
            alpha.re[k] = e.wr1;
            alpha.re[k + 1] = e.wr2;
            alpha.im[k] = e.wi;
            alpha.im[k + 1] = -e.wi;
            k += 1;
            continue;
        }
        if (std.math.signbit(b[k + k * p.ldb])) {
            for (0..n) |i| {
                a[k + i * p.lda] = -a[k + i * p.lda];
                b[k + i * p.ldb] = -b[k + i * p.ldb];
                if (p.q) |q| q[i + k * p.ldq] = -q[i + k * p.ldq];
            }
        }
        alpha.re[k] = a[k + k * p.lda];
        alpha.im[k] = 0;
        beta[k] = b[k + k * p.ldb];
    }
}

/// `ztgsen.f`'s last loop: `B(k, k) := |B(k, k)|`, row `k` of `A` and `B`
/// (from `k`) scaled by its phase's conjugate and column `k` of `Q` by the
/// phase; 0 below the smallest normal.
fn normalizeComplex(comptime T: type, p: Pencil(T), alpha: []T, beta: []T) void {
    const F = Float(T);
    const safmin = std.math.floatMin(F);
    for (0..p.n) |k| {
        const bkk = &p.b[k + k * p.ldb];
        const dscale = ops.abs(T, bkk.*);
        if (dscale > safmin) {
            const temp2 = ops.divReal(T, bkk.*, dscale);
            const temp1 = ops.conj(T, temp2);
            bkk.* = ops.scalar(T, dscale);
            for (k + 1..p.n) |j| p.b[k + j * p.ldb] = ops.mul(T, temp1, p.b[k + j * p.ldb]);
            for (k..p.n) |j| p.a[k + j * p.lda] = ops.mul(T, temp1, p.a[k + j * p.lda]);
            if (p.q) |q| for (0..p.n) |i| {
                q[i + k * p.ldq] = ops.mul(T, temp2, q[i + k * p.ldq]);
            };
        } else bkk.* = ops.zero(T);
        alpha[k] = p.a[k + k * p.lda];
        beta[k] = bkk.*;
    }
}
