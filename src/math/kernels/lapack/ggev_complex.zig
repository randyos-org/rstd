//! `?GGEV` for complex pencils: faer v0.24.4's `compute_gevd_generic` and
//! `qz_to_gevd_cplx` (`faer/src/linalg/gevd/mod.rs`, MIT): `B`'s QR, `gghrd`,
//! the QZ, then the eigenvectors from the generalized Schur form, each of
//! unit 2-norm. Its `solve_shifted_1x1` follows Reference-LAPACK's
//! `ztgevc.f` (`DMIN` perturbation) where faer departs from it (issue #766).

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const ops = @import("ops.zig");
const interface = @import("interface.zig");
const qr = @import("qr.zig");
const gghrd_reduce = @import("gghrd.zig");
const qz_complex = @import("qz_complex.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const scalar = ops.scalar;
const add = ops.add;
const sub = ops.sub;
const mul = ops.mul;
const conj = ops.conj;
const mulReal = ops.mulReal;

pub const VectorJob = interface.VectorJob;

/// What `ggevComplex` and `ggevReal` return: a slice too short or a leading dimension too
/// small; `NoConvergence` for a NaN or infinity in `A` or `B`, or when the QZ
/// does not converge; the allocator's and cancellation.
pub const GgevError = interface.ShapeError || error{NoConvergence} || Allocator.Error || Io.Cancelable;

/// `?GGEV`, complex: the generalized eigenvalues `alpha[j] / beta[j]` of
/// `(A, B)` (`n x n`; both overwritten) and, as `jobvl` and `jobvr` ask, the
/// left (`yᴴ A = λ yᴴ B`) and right (`A x = λ B x`) eigenvectors, the
/// columns of `vl` and `vr`, each of unit 2-norm. An eigenvalue with
/// `beta[j] = 0` is infinite.
pub fn ggevComplex(allocator: Allocator, io: Io, comptime T: type, jobvl: VectorJob, jobvr: VectorJob, n: usize, a: []T, lda: usize, b: []T, ldb: usize, alpha: []T, beta: []T, vl: []T, ldvl: usize, vr: []T, ldvr: usize) GgevError!void {
    comptime ops.requireComplex(T, "ggevComplex");
    try util.holds(a.len, n, n, lda);
    try util.holds(b.len, n, n, ldb);
    if (alpha.len < n or beta.len < n) return error.SliceTooShort;
    const left = jobvl == .compute;
    const right = jobvr == .compute;
    if (left) try util.holds(vl.len, n, n, ldvl);
    if (right) try util.holds(vr.len, n, n, ldvr);
    if (n == 0) return;
    try reduce(allocator, io, T, left, right, n, a, lda, b, ldb, vl, ldvl, vr, ldvr);
    const job: qz_complex.SchurJob = if (left or right) .schur_form else .eigenvalues;
    qz_complex.hgeqzComplex(T, job, if (left) .accumulate else .none, if (right) .accumulate else .none, n, 0, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr) catch |e| return ops.narrow(GgevError, e);
    if (left or right) try eigenvectors(allocator, T, n, a, lda, b, ldb, if (left) vl else null, ldvl, if (right) vr else null, ldvr);
}

/// The steps before the QZ, real or complex (shapes checked, `n > 0`):
/// `NoConvergence` for a NaN or infinity; `B := R` and `A := Qᴴ A` by `B`'s
/// QR, `vl := Q` (`left`), then `gghrd` (`vr` from the identity, `right`).
pub fn reduce(allocator: Allocator, io: Io, comptime T: type, left: bool, right: bool, n: usize, a: []T, lda: usize, b: []T, ldb: usize, vl: []T, ldvl: usize, vr: []T, ldvr: usize) GgevError!void {
    for (0..n) |j| for (0..n) |i| {
        if (!finite(T, a[i + j * lda]) or !finite(T, b[i + j * ldb])) return error.NoConvergence;
    };
    const tau = try allocator.alloc(T, n);
    defer allocator.free(tau);
    qr.geqrf(allocator, io, T, n, n, b, ldb, tau) catch |e| return ops.narrow(GgevError, e);
    qr.ormqr(allocator, io, T, .l, .c, n, n, n, b, ldb, tau, a, lda) catch |e| return ops.narrow(GgevError, e);
    if (left) {
        for (0..n) |j| for (0..n) |i| {
            vl[i + j * ldvl] = scalar(T, if (i == j) 1 else 0);
        };
        qr.ormqr(allocator, io, T, .l, .n, n, n, n, b, ldb, tau, vl, ldvl) catch |e| return ops.narrow(GgevError, e);
    }
    gghrd_reduce.gghrd(T, if (left) .accumulate else .none, if (right) .from_identity else .none, n, 0, n, a, lda, b, ldb, vl, ldvl, vr, ldvr) catch unreachable; // shapes checked by the caller
}

fn finite(comptime T: type, x: T) bool {
    if (comptime !types.isComplex(T)) return std.math.isFinite(x);
    return std.math.isFinite(x.re) and std.math.isFinite(x.im);
}

/// faer's `qz_to_gevd_cplx`: the eigenvectors of the triangular pencil `(S,
/// P)` brought back through `Q` (`u_left`) and `Z` (`u_right`).
fn eigenvectors(allocator: Allocator, comptime T: type, n: usize, a: []const T, lda: usize, b: []const T, ldb: usize, u_left: ?[]T, ldl: usize, u_right: ?[]T, ldr: usize) Allocator.Error!void {
    const F = Float(T);
    const ulp = std.math.floatEps(F);
    const safmin = std.math.floatMin(F);
    const smallnum = safmin * @as(F, @floatFromInt(n));
    const small = smallnum / ulp;
    const bignum = 1 / smallnum;
    const big = 1 / small;
    const work = try allocator.alloc(F, 2 * n);
    defer allocator.free(work);
    const acolnorm = work[0..n];
    const bcolnorm = work[n..];
    const vecs = try allocator.alloc(T, 2 * n);
    defer allocator.free(vecs);
    const rhs = vecs[0..n];
    const tmp = vecs[n..];
    const at = struct {
        fn f(m: []const T, ld: usize, i: usize, j: usize) T {
            return m[i + j * ld];
        }
    }.f;
    var anorm: F = 0;
    var bnorm: F = 0;
    for (0..n) |j| {
        var sa: F = 0;
        var sb: F = 0;
        for (0..j) |i| {
            sa += ops.abs(T, at(a, lda, i, j));
            sb += ops.abs(T, at(b, ldb, i, j));
        }
        acolnorm[j] = sa;
        bcolnorm[j] = sb;
        anorm = @max(anorm, sa + ops.abs(T, at(a, lda, j, j)));
        bnorm = @max(bnorm, sb + ops.abs(T, at(b, ldb, j, j)));
    }
    const ascale = 1 / @max(anorm, safmin);
    const bscale = 1 / @max(bnorm, safmin);
    const Coefs = struct { acoef: F, bcoef: T };
    const coefs = struct {
        /// The scaled `(beta, alpha)` of eigenvalue `je`, kept from
        /// overflow: faer's `acoef`, `bcoef`.
        fn f(ajj: T, bjj: T, as: F, bs: F, an: F, bn: F, sm: F, bg: F) Coefs {
            const temp = @max(safmin, @max(ops.abs(T, ajj) * as, ops.abs(T, bjj) * bs));
            const salfar = mulReal(T, mulReal(T, ajj, temp), as);
            const sbeta = bjj.re * temp * bs;
            var acoef = sbeta * as;
            var bcoef = mulReal(T, salfar, bs);
            var scale: F = 1;
            const lsa = @abs(sbeta) >= safmin and @abs(acoef) < sm;
            const lsb = ops.abs(T, salfar) >= safmin and ops.abs(T, bcoef) < sm;
            if (lsa) scale = (sm / @abs(sbeta)) * @min(an, bg);
            if (lsb) scale = @max(scale, (sm / ops.abs(T, salfar)) * @min(bn, bg));
            if (lsa or lsb) {
                scale = @min(scale, 1 / (safmin * @max(1, @max(@abs(acoef), ops.abs(T, bcoef)))));
                acoef = if (lsa) sbeta * scale * as else acoef * scale;
                bcoef = if (lsb) mulReal(T, mulReal(T, salfar, scale), bs) else mulReal(T, bcoef, scale);
            }
            return .{ .acoef = acoef, .bcoef = bcoef };
        }
    }.f;

    if (u_left) |u| {
        for (0..n) |je| {
            if (@max(ops.abs(T, at(a, lda, je, je)), ops.abs(T, at(b, ldb, je, je))) < safmin) {
                for (0..n) |i| u[i + je * ldl] = scalar(T, if (i == je) 1 else 0);
                continue;
            }
            const c = coefs(at(a, lda, je, je), at(b, ldb, je, je), ascale, bscale, anorm, bnorm, small, big);
            const acoefa = @abs(c.acoef);
            const bcoefa = ops.abs(T, c.bcoef);
            @memset(rhs, scalar(T, 0));
            rhs[je] = scalar(T, 1);
            var xmax: F = 1;
            const dmin = @max(safmin, @max(ulp * acoefa * anorm, ulp * bcoefa * bnorm));
            for (je + 1..n) |j| {
                const xscale = 1 / xmax;
                const temp = @max(@max(acolnorm[j], bcolnorm[j]), acoefa * acolnorm[j] + bcoefa * bcolnorm[j]);
                if (temp > bignum * xscale) {
                    for (je..j) |jr| rhs[jr] = mulReal(T, rhs[jr], xscale);
                    xmax *= xscale;
                }
                var sums = scalar(T, 0);
                var sump = scalar(T, 0);
                for (je..j) |jr| {
                    sums = add(T, sums, ops.mulConj(T, at(a, lda, jr, j), rhs[jr]));
                    sump = add(T, sump, ops.mulConj(T, at(b, ldb, jr, j), rhs[jr]));
                }
                rhs[j] = sub(T, mul(T, conj(T, c.bcoef), sump), mulReal(T, sums, c.acoef));
                solveShifted1x1(T, dmin, c.acoef, conj(T, at(a, lda, j, j)), conj(T, at(b, ldb, j, j)), &rhs[j], conj(T, c.bcoef));
            }
            // `u(:, je) = u(:, je..) rhs(je..)`, normalized
            for (0..n) |i| {
                var s = scalar(T, 0);
                for (je..n) |k| s = add(T, s, mul(T, u[i + k * ldl], rhs[k]));
                tmp[i] = s;
            }
            normalizeInto(T, u[je * ldl ..][0..n], tmp);
        }
    }
    if (u_right) |u| {
        var je = n;
        while (je > 0) {
            je -= 1;
            if (@max(ops.abs(T, at(a, lda, je, je)), ops.abs(T, at(b, ldb, je, je))) < safmin) {
                for (0..n) |i| u[i + je * ldr] = scalar(T, if (i == je) 1 else 0);
                continue;
            }
            const c = coefs(at(a, lda, je, je), at(b, ldb, je, je), ascale, bscale, anorm, bnorm, small, big);
            const acoefa = @abs(c.acoef);
            const bcoefa = ops.abs(T, c.bcoef);
            @memset(rhs, scalar(T, 0));
            rhs[je] = scalar(T, 1);
            for (0..je) |jr| rhs[jr] = sub(T, mul(T, c.bcoef, at(b, ldb, jr, je)), mulReal(T, at(a, lda, jr, je), c.acoef));
            const dmin = @max(safmin, @max(ulp * acoefa * anorm, ulp * bcoefa * bnorm));
            var j = je;
            while (j > 0) {
                j -= 1;
                solveShifted1x1(T, dmin, c.acoef, at(a, lda, j, j), at(b, ldb, j, j), &rhs[j], c.bcoef);
                const creala = mulReal(T, rhs[j], c.acoef);
                const crealb = mul(T, c.bcoef, rhs[j]);
                for (0..j) |jr| rhs[jr] = add(T, sub(T, rhs[jr], mul(T, creala, at(a, lda, jr, j))), mul(T, crealb, at(b, ldb, jr, j)));
            }
            for (0..n) |i| {
                var s = scalar(T, 0);
                for (0..je + 1) |k| s = add(T, s, mul(T, u[i + k * ldr], rhs[k]));
                tmp[i] = s;
            }
            normalizeInto(T, u[je * ldr ..][0..n], tmp);
        }
    }
}

/// `dst = src / ‖src‖₂`.
fn normalizeInto(comptime T: type, dst: []T, src: []const T) void {
    var ss: Float(T) = 0;
    for (src) |x| ss = std.math.hypot(ss, ops.abs(T, x));
    const inv = 1 / ss;
    for (dst, src) |*d, x| d.* = mulReal(T, x, inv);
}

/// `x := x / (ca a - w d0)`, the divisor raised to `smin` when smaller, as
/// `ztgevc.f`'s `DMIN` perturbation does (faer divides by both; issue #766).
fn solveShifted1x1(comptime T: type, smin_in: Float(T), ca: Float(T), a: T, d0: T, x: *T, w: T) void {
    const safmin = std.math.floatMin(Float(T));
    const smin = @max(smin_in, safmin + safmin);
    const cr = sub(T, mulReal(T, a, ca), mul(T, w, d0));
    if (ops.abs(T, cr) < smin) {
        x.* = mulReal(T, x.*, 1 / smin);
        return;
    }
    x.* = mul(T, x.*, ops.div(T, scalar(T, 1), cr));
}
