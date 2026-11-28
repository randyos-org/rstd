//! `xTGEX2`: swapping two adjacent diagonal blocks of a generalized Schur
//! form, Reference-LAPACK 3.12.0 (`lapack-netlib/SRC/dtgex2.f`, `ztgex2.f`,
//! OpenBLAS v0.3.34-317; modified BSD), the weak and strong stability tests
//! (`WANDS`) included.

const std = @import("std");
const blas = @import("../blas/__root__.zig");
const types = @import("../blas/types.zig");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");
const norm = @import("norm.zig");
const auxiliary = @import("auxiliary.zig");
const tgsy2_mod = @import("tgsy2.zig");
const qr_small = @import("qr_small.zig");
const lagv2_mod = @import("lagv2.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const Error = ops.Error;

/// `TWENTY`: the stability tests' multiple of `eps` times the norm.
const twenty = 20;

/// The Frobenius norm of the `m x m` block at `a` (`ld` apart) by `lassq`,
/// in column order as the reference's packed copy.
fn frob(comptime T: type, m: usize, a: [*]const T, ld: usize) Float(T) {
    var scale: Float(T) = 0;
    var sum: Float(T) = 1;
    for (0..m) |j| norm.lassq(T, .of(a[j * ld ..][0..m]), &scale, &sum);
    return scale * @sqrt(sum);
}

/// `ZTGEX2` (`ztgex2.f`): swaps the 1 x 1 blocks at rows `j1` and `j1 + 1`
/// (0-based) of the `n x n` triangular `(A, B)` by unitary equivalence, `Q`
/// and `Z` updated as asked; false for the reference's `INFO = 1` (the swap
/// failed the stability tests and nothing changed).
pub fn tgex2Complex(comptime T: type, wantq: bool, wantz: bool, n: usize, a: [*]T, lda: usize, b: [*]T, ldb: usize, q: [*]T, ldq: usize, z: [*]T, ldz: usize, j1: usize) bool {
    if (n <= 1) return true;
    const F = Float(T);
    const abs = ops.abs;
    const zrot = blas_calls.zrot;
    // `S` and `T`, 2 x 2 column-major
    var s = [4]T{ a[j1 + j1 * lda], a[j1 + 1 + j1 * lda], a[j1 + (j1 + 1) * lda], a[j1 + 1 + (j1 + 1) * lda] };
    var t = [4]T{ b[j1 + j1 * ldb], b[j1 + 1 + j1 * ldb], b[j1 + (j1 + 1) * ldb], b[j1 + 1 + (j1 + 1) * ldb] };
    const eps = std.math.floatEps(F);
    const smlnum = std.math.floatMin(F) / eps;
    const thresha = @max(twenty * eps * frob(T, 2, &s, 2), smlnum);
    const threshb = @max(twenty * eps * frob(T, 2, &t, 2), smlnum);
    const f = ops.sub(T, ops.mul(T, s[3], t[0]), ops.mul(T, t[3], s[0]));
    var g = ops.sub(T, ops.mul(T, s[3], t[2]), ops.mul(T, t[3], s[2]));
    const sa = abs(T, s[3]) * abs(T, t[0]);
    const sb = abs(T, s[0]) * abs(T, t[3]);
    var cz: F = undefined;
    var sz: T = undefined;
    blas.zrotg(T, &g, f, &cz, &sz);
    sz = ops.neg(T, sz);
    zrot(T, 2, &s, 1, s[2..].ptr, 1, cz, ops.conj(T, sz));
    zrot(T, 2, &t, 1, t[2..].ptr, 1, cz, ops.conj(T, sz));
    var cq: F = undefined;
    var sq: T = undefined;
    var r = if (sa >= sb) s[0] else t[0];
    blas.zrotg(T, &r, if (sa >= sb) s[1] else t[1], &cq, &sq);
    zrot(T, 2, &s, 2, s[1..].ptr, 2, cq, sq);
    zrot(T, 2, &t, 2, t[1..].ptr, 2, cq, sq);
    const weak = abs(T, s[1]) <= thresha and abs(T, t[1]) <= threshb;
    if (!weak) return false;
    // strong: undoing both rotations gives the blocks back
    var w = s ++ t;
    zrot(T, 2, &w, 1, w[2..].ptr, 1, cz, ops.neg(T, ops.conj(T, sz)));
    zrot(T, 2, w[4..].ptr, 1, w[6..].ptr, 1, cz, ops.neg(T, ops.conj(T, sz)));
    zrot(T, 2, &w, 2, w[1..].ptr, 2, cq, ops.neg(T, sq));
    zrot(T, 2, w[4..].ptr, 2, w[5..].ptr, 2, cq, ops.neg(T, sq));
    for (0..2) |i| {
        w[i] = ops.sub(T, w[i], a[j1 + i + j1 * lda]);
        w[i + 2] = ops.sub(T, w[i + 2], a[j1 + i + (j1 + 1) * lda]);
        w[i + 4] = ops.sub(T, w[i + 4], b[j1 + i + j1 * ldb]);
        w[i + 6] = ops.sub(T, w[i + 6], b[j1 + i + (j1 + 1) * ldb]);
    }
    const strong = frob(T, 2, &w, 2) <= thresha and frob(T, 2, w[4..].ptr, 2) <= threshb;
    if (!strong) return false;
    zrot(T, j1 + 2, a + j1 * lda, 1, a + (j1 + 1) * lda, 1, cz, ops.conj(T, sz));
    zrot(T, j1 + 2, b + j1 * ldb, 1, b + (j1 + 1) * ldb, 1, cz, ops.conj(T, sz));
    zrot(T, n - j1, a + j1 + j1 * lda, lda, a + j1 + 1 + j1 * lda, lda, cq, sq);
    zrot(T, n - j1, b + j1 + j1 * ldb, ldb, b + j1 + 1 + j1 * ldb, ldb, cq, sq);
    a[j1 + 1 + j1 * lda] = ops.zero(T);
    b[j1 + 1 + j1 * ldb] = ops.zero(T);
    if (wantz) zrot(T, n, z + j1 * ldz, 1, z + (j1 + 1) * ldz, 1, cz, ops.conj(T, sz));
    if (wantq) zrot(T, n, q + j1 * ldq, 1, q + (j1 + 1) * ldq, 1, cq, ops.conj(T, sq));
    return true;
}

/// `LDST`: the leading dimension of `dtgex2.f`'s 4 x 4 work blocks.
const ldst = 4;

/// `C := alpha op(A) op(B) + beta C` (`DGEMM` on the small blocks), `m x n`,
/// each sum in `k`'s order; `C` is not read when `beta` is 0.
fn gemm(comptime F: type, ta: bool, tb: bool, m: usize, n: usize, k: usize, alpha: F, a: [*]const F, lda: usize, b: [*]const F, ldb: usize, beta: F, c: [*]F, ldc: usize) void {
    for (0..n) |j| for (0..m) |i| {
        var sum: F = 0;
        for (0..k) |l| {
            const x = if (ta) a[l + i * lda] else a[i + l * lda];
            const y = if (tb) b[j + l * ldb] else b[l + j * ldb];
            sum += x * y;
        }
        c[i + j * ldc] = alpha * sum + if (beta == 0) 0 else beta * c[i + j * ldc];
    };
}

/// `m x n` from `src` (`lds` apart) to `dst` (`ldd`), `DLACPY`'s full copy.
fn lacpy(comptime F: type, m: usize, n: usize, src: [*]const F, lds: usize, dst: [*]F, ldd: usize) void {
    for (0..n) |j| @memcpy(dst[j * ldd ..][0..m], src[j * lds ..][0..m]);
}

/// `DTGEX2` (`dtgex2.f`): swaps the `n1 x n1` and `n2 x n2` diagonal blocks
/// (each 1 or 2) at row `j1` (0-based) of the `n x n` real generalized Schur
/// form `(A, B)` by orthogonal equivalence, `Q` and `Z` updated as asked;
/// false for the reference's `INFO = 1` (the swap failed the stability
/// tests and nothing changed).
pub fn tgex2Real(al: Allocator, io: Io, comptime F: type, wantq: bool, wantz: bool, n: usize, a: [*]F, lda: usize, b: [*]F, ldb: usize, q: [*]F, ldq: usize, z: [*]F, ldz: usize, j1: usize, n1: usize, n2: usize) Error!bool {
    if (n <= 1 or n1 == 0 or n2 == 0) return true;
    if (n1 > n or j1 + n1 >= n) return true;
    const m = n1 + n2;
    const work = try al.alloc(F, @max(n * m, 2 * m * m));
    defer al.free(work);
    const at = struct {
        fn p(x: [*]F, ld: usize, i: usize, j: usize) *F {
            return &x[i + j * ld];
        }
    }.p;
    var li: [ldst * ldst]F = @splat(0);
    var ir: [ldst * ldst]F = @splat(0);
    var s: [ldst * ldst]F = @splat(0);
    var t: [ldst * ldst]F = @splat(0);
    const aj = a + j1 + j1 * lda;
    const bj = b + j1 + j1 * ldb;
    lacpy(F, m, m, aj, lda, &s, ldst);
    lacpy(F, m, m, bj, ldb, &t, ldst);
    const eps = std.math.floatEps(F);
    const smlnum = std.math.floatMin(F) / eps;
    const thresha = @max(twenty * eps * frob(F, m, &s, ldst), smlnum);
    const threshb = @max(twenty * eps * frob(F, m, &t, ldst), smlnum);
    const rot = blas_calls.rot;
    if (m == 2) {
        const f = at(&s, ldst, 1, 1).* * at(&t, ldst, 0, 0).* - at(&t, ldst, 1, 1).* * at(&s, ldst, 0, 0).*;
        const g = at(&s, ldst, 1, 1).* * at(&t, ldst, 0, 1).* - at(&t, ldst, 1, 1).* * at(&s, ldst, 0, 1).*;
        const sa = @abs(at(&s, ldst, 1, 1).*) * @abs(at(&t, ldst, 0, 0).*);
        const sb = @abs(at(&s, ldst, 0, 0).*) * @abs(at(&t, ldst, 1, 1).*);
        const gr = auxiliary.lartg(F, f, g);
        at(&ir, ldst, 0, 1).* = gr.c;
        at(&ir, ldst, 0, 0).* = gr.s;
        at(&ir, ldst, 1, 0).* = -at(&ir, ldst, 0, 1).*;
        at(&ir, ldst, 1, 1).* = at(&ir, ldst, 0, 0).*;
        const irc = at(&ir, ldst, 0, 0).*;
        const irs = at(&ir, ldst, 1, 0).*;
        try rot(io, F, 2, &s, 1, s[ldst..].ptr, 1, irc, irs);
        try rot(io, F, 2, &t, 1, t[ldst..].ptr, 1, irc, irs);
        const gl = if (sa >= sb) auxiliary.lartg(F, s[0], s[1]) else auxiliary.lartg(F, t[0], t[1]);
        at(&li, ldst, 0, 0).* = gl.c;
        at(&li, ldst, 1, 0).* = gl.s;
        try rot(io, F, 2, &s, ldst, s[1..].ptr, ldst, gl.c, gl.s);
        try rot(io, F, 2, &t, ldst, t[1..].ptr, ldst, gl.c, gl.s);
        at(&li, ldst, 1, 1).* = gl.c;
        at(&li, ldst, 0, 1).* = -gl.s;
        const weak = @abs(at(&s, ldst, 1, 0).*) <= thresha and @abs(at(&t, ldst, 1, 0).*) <= threshb;
        if (!weak) return false;
        // strong: `LI S IRᵀ` gives the blocks back
        const w2 = work[m * m ..].ptr;
        lacpy(F, m, m, aj, lda, w2, m);
        gemm(F, false, false, m, m, m, 1, &li, ldst, &s, ldst, 0, work.ptr, m);
        gemm(F, false, true, m, m, m, -1, work.ptr, m, &ir, ldst, 1, w2, m);
        const da = frob(F, m, w2, m);
        lacpy(F, m, m, bj, ldb, w2, m);
        gemm(F, false, false, m, m, m, 1, &li, ldst, &t, ldst, 0, work.ptr, m);
        gemm(F, false, true, m, m, m, -1, work.ptr, m, &ir, ldst, 1, w2, m);
        const db = frob(F, m, w2, m);
        if (!(da <= thresha and db <= threshb)) return false;
        try rot(io, F, j1 + 2, a + j1 * lda, 1, a + (j1 + 1) * lda, 1, irc, irs);
        try rot(io, F, j1 + 2, b + j1 * ldb, 1, b + (j1 + 1) * ldb, 1, irc, irs);
        try rot(io, F, n - j1, aj, lda, aj + 1, lda, gl.c, gl.s);
        try rot(io, F, n - j1, bj, ldb, bj + 1, ldb, gl.c, gl.s);
        aj[1] = 0;
        bj[1] = 0;
        if (wantz) try rot(io, F, n, z + j1 * ldz, 1, z + (j1 + 1) * ldz, 1, irc, irs);
        if (wantq) try rot(io, F, n, q + j1 * ldq, 1, q + (j1 + 1) * ldq, 1, gl.c, gl.s);
        return true;
    }
    // `(A11, B11) X - Y (A22, B22) = scale (A12, B12)` for the subspaces
    lacpy(F, n1, n2, t[n1 * ldst ..].ptr, ldst, &li, ldst);
    lacpy(F, n1, n2, s[n1 * ldst ..].ptr, ldst, ir[n2 + n1 * ldst ..].ptr, ldst);
    var dif: tgsy2_mod.Dif(F) = .{};
    const sol = try tgsy2_mod.tgsy2(al, io, F, false, 0, n1, n2, &s, ldst, s[n1 + n1 * ldst ..].ptr, ldst, ir[n2 + n1 * ldst ..].ptr, ldst, &t, ldst, t[n1 + n1 * ldst ..].ptr, ldst, &li, ldst, &dif);
    if (sol.perturbed) return false;
    const scale = sol.scale;
    var taul: [ldst]F = undefined;
    var taur: [ldst]F = undefined;
    // `[-X; scale I] = Q_L R`, `[scale I, Y] = R Q_R`
    for (0..n2) |i| {
        for (0..n1) |l| at(&li, ldst, l, i).* = -at(&li, ldst, l, i).*;
        at(&li, ldst, n1 + i, i).* = scale;
    }
    try qr_small.geqr2(io, F, m, n2, &li, ldst, &taul);
    qr_small.org2r(F, m, m, n2, &li, ldst, &taul);
    for (0..n1) |i| at(&ir, ldst, n2 + i, i).* = scale;
    try qr_small.gerq2(io, F, n1, m, ir[n2..].ptr, ldst, &taur);
    qr_small.orgr2(F, m, m, n1, &ir, ldst, &taur);
    // `S := LIᵀ S IRᵀ`, `T := LIᵀ T IRᵀ`
    gemm(F, true, false, m, m, m, 1, &li, ldst, &s, ldst, 0, work.ptr, m);
    gemm(F, false, true, m, m, m, 1, work.ptr, m, &ir, ldst, 0, &s, ldst);
    gemm(F, true, false, m, m, m, 1, &li, ldst, &t, ldst, 0, work.ptr, m);
    gemm(F, false, true, m, m, m, 1, work.ptr, m, &ir, ldst, 0, &t, ldst);
    var scpy = s;
    var tcpy = t;
    const ircop = ir;
    var licop = li;
    // triangularize `T` by RQ, and the copy by QR: keep the one leaving
    // `S` closer to block upper triangular
    try qr_small.gerq2(io, F, m, m, &t, ldst, &taur);
    qr_small.ormr2(F, false, true, m, m, m, &t, ldst, &taur, &s, ldst);
    qr_small.ormr2(F, true, false, m, m, m, &t, ldst, &taur, &ir, ldst);
    var dscale: F = 0;
    var dsum: F = 1;
    for (0..n2) |i| norm.lassq(F, .of(s[n2 + i * ldst ..][0..n1]), &dscale, &dsum);
    const brqa21 = dscale * @sqrt(dsum);
    try qr_small.geqr2(io, F, m, m, &tcpy, ldst, &taul);
    qr_small.orm2r(F, true, true, m, m, m, &tcpy, ldst, &taul, &scpy, ldst);
    qr_small.orm2r(F, false, false, m, m, m, &tcpy, ldst, &taul, &licop, ldst);
    dscale = 0;
    dsum = 1;
    for (0..n2) |i| norm.lassq(F, .of(scpy[n2 + i * ldst ..][0..n1]), &dscale, &dsum);
    const bqra21 = dscale * @sqrt(dsum);
    if (bqra21 <= brqa21 and bqra21 <= thresha) {
        s = scpy;
        t = tcpy;
        ir = ircop;
        li = licop;
    } else if (brqa21 >= thresha) return false;
    for (0..m - 1) |j| for (j + 1..m) |i| {
        at(&t, ldst, i, j).* = 0;
    };
    // strong: `LI S IR` gives the blocks back
    const w2 = work[m * m ..].ptr;
    lacpy(F, m, m, aj, lda, w2, m);
    gemm(F, false, false, m, m, m, 1, &li, ldst, &s, ldst, 0, work.ptr, m);
    gemm(F, false, false, m, m, m, -1, work.ptr, m, &ir, ldst, 1, w2, m);
    const da = frob(F, m, w2, m);
    lacpy(F, m, m, bj, ldb, w2, m);
    gemm(F, false, false, m, m, m, 1, &li, ldst, &t, ldst, 0, work.ptr, m);
    gemm(F, false, false, m, m, m, -1, work.ptr, m, &ir, ldst, 1, w2, m);
    const db = frob(F, m, w2, m);
    if (!(da <= thresha and db <= threshb)) return false;
    for (0..n2) |j| for (n2..m) |i| {
        at(&s, ldst, i, j).* = 0;
    };
    lacpy(F, m, m, &s, ldst, aj, lda);
    lacpy(F, m, m, &t, ldst, bj, ldb);
    // standardize the new 2 x 2 blocks: `work` (`m x m`) the left
    // rotations, `t` the right
    t = @splat(0);
    @memset(work[0 .. m * m], 0);
    work[0] = 1;
    t[0] = 1;
    if (n2 > 1) {
        const g = lagv2_mod.lagv2(F, aj, lda, bj, ldb);
        work[0] = g.csl;
        work[1] = g.snl;
        at(&t, ldst, 0, 0).* = g.csr;
        at(&t, ldst, 1, 0).* = g.snr;
        work[m] = -work[1];
        work[m + 1] = work[0];
        at(&t, ldst, n2 - 1, n2 - 1).* = at(&t, ldst, 0, 0).*;
        at(&t, ldst, 0, 1).* = -at(&t, ldst, 1, 0).*;
    }
    work[m * m - 1] = 1;
    at(&t, ldst, m - 1, m - 1).* = 1;
    if (n1 > 1) {
        const g = lagv2_mod.lagv2(F, aj + n2 + n2 * lda, lda, bj + n2 + n2 * ldb, ldb);
        work[n2 * m + n2] = g.csl;
        work[n2 * m + n2 + 1] = g.snl;
        at(&t, ldst, n2, n2).* = g.csr;
        at(&t, ldst, m - 1, m - 2).* = g.snr;
        work[m * m - 1] = work[n2 * m + n2];
        work[m * m - 2] = -work[n2 * m + n2 + 1];
        at(&t, ldst, m - 1, m - 1).* = at(&t, ldst, n2, n2).*;
        at(&t, ldst, m - 2, m - 1).* = -at(&t, ldst, m - 1, m - 2).*;
    }
    const tmp = work[m * m ..].ptr;
    gemm(F, true, false, n2, n1, n2, 1, work.ptr, m, aj + n2 * lda, lda, 0, tmp, n2);
    lacpy(F, n2, n1, tmp, n2, aj + n2 * lda, lda);
    gemm(F, true, false, n2, n1, n2, 1, work.ptr, m, bj + n2 * ldb, ldb, 0, tmp, n2);
    lacpy(F, n2, n1, tmp, n2, bj + n2 * ldb, ldb);
    gemm(F, false, false, m, m, m, 1, &li, ldst, work.ptr, m, 0, tmp, m);
    lacpy(F, m, m, tmp, m, &li, ldst);
    gemm(F, false, false, n2, n1, n1, 1, aj + n2 * lda, lda, t[n2 + n2 * ldst ..].ptr, ldst, 0, work.ptr, n2);
    lacpy(F, n2, n1, work.ptr, n2, aj + n2 * lda, lda);
    gemm(F, false, false, n2, n1, n1, 1, bj + n2 * ldb, ldb, t[n2 + n2 * ldst ..].ptr, ldst, 0, work.ptr, n2);
    lacpy(F, n2, n1, work.ptr, n2, bj + n2 * ldb, ldb);
    gemm(F, true, false, m, m, m, 1, &ir, ldst, &t, ldst, 0, work.ptr, m);
    lacpy(F, m, m, work.ptr, m, &ir, ldst);
    // the rest of `A`, `B`, `Q`, `Z`
    if (wantq) {
        gemm(F, false, false, n, m, m, 1, q + j1 * ldq, ldq, &li, ldst, 0, work.ptr, n);
        lacpy(F, n, m, work.ptr, n, q + j1 * ldq, ldq);
    }
    if (wantz) {
        gemm(F, false, false, n, m, m, 1, z + j1 * ldz, ldz, &ir, ldst, 0, work.ptr, n);
        lacpy(F, n, m, work.ptr, n, z + j1 * ldz, ldz);
    }
    const i = j1 + m;
    if (i < n) {
        gemm(F, true, false, m, n - i, m, 1, &li, ldst, a + j1 + i * lda, lda, 0, work.ptr, m);
        lacpy(F, m, n - i, work.ptr, m, a + j1 + i * lda, lda);
        gemm(F, true, false, m, n - i, m, 1, &li, ldst, b + j1 + i * ldb, ldb, 0, work.ptr, m);
        lacpy(F, m, n - i, work.ptr, m, b + j1 + i * ldb, ldb);
    }
    if (j1 > 0) {
        gemm(F, false, false, j1, m, m, 1, a + j1 * lda, lda, &ir, ldst, 0, work.ptr, j1);
        lacpy(F, j1, m, work.ptr, j1, a + j1 * lda, lda);
        gemm(F, false, false, j1, m, m, 1, b + j1 * ldb, ldb, &ir, ldst, 0, work.ptr, j1);
        lacpy(F, j1, m, work.ptr, j1, b + j1 * ldb, ldb);
    }
    return true;
}
