//! Reference-LAPACK's small auxiliaries the routines share: `xLASCL`,
//! `xLANGE('M')` and the drivers' scaling prologue, `xLARTG`, `xLARFG`,
//! `xLAPY2`, `xLAPY3` and `xLASET`.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");

const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const Uplo = types.Uplo;
const Mat = ops.Mat;
const Error = ops.Error;
const scalar = ops.scalar;
const zero = ops.zero;
const one = ops.one;
const re = ops.re;
const mulReal = ops.mulReal;
const div = ops.div;
const abs = ops.abs;
const requireFloat = ops.requireFloat;
const nrm2 = blas_calls.nrm2;
const scal = blas_calls.scal;
const scalComplex = blas_calls.scalComplex;

/// `xLASCL` type `G` (with `upper`, type `U`: the upper triangle): `A := A *
/// (to / from)`, in steps that neither overflow nor underflow
/// (`lapack-netlib/SRC/dlascl.f`). `from` is finite and nonzero.
pub fn lascl(comptime T: type, upper: bool, from: Float(T), to: Float(T), a: Mat(T)) void {
    const F = Float(T);
    const smlnum = std.math.floatMin(F);
    const bignum = 1 / smlnum;
    var cfromc = from;
    var ctoc = to;
    while (true) {
        const cfrom1 = cfromc * smlnum;
        var f: F = undefined;
        var done = true;
        if (cfrom1 == cfromc) {
            f = ctoc / cfromc;
        } else {
            const cto1 = ctoc / bignum;
            if (cto1 == ctoc) {
                f = ctoc;
                cfromc = 1;
            } else if (@abs(cfrom1) > @abs(ctoc) and ctoc != 0) {
                f = smlnum;
                done = false;
                cfromc = cfrom1;
            } else if (@abs(cto1) > @abs(cfromc)) {
                f = bignum;
                done = false;
                ctoc = cto1;
            } else {
                f = ctoc / cfromc;
                if (f == 1) return;
            }
        }
        for (0..a.cols) |j| for (a.col(j)[0..if (upper) @min(j + 1, a.rows) else a.rows]) |*x| {
            x.* = mulReal(T, x.*, f);
        };
        if (done) return;
    }
}

/// `xLANGE('M')`: the largest `|a(i, j)|`, NaN if any is NaN.
pub fn normMax(comptime T: type, a: Mat(T)) Float(T) {
    var x: Float(T) = 0;
    for (0..a.cols) |j| for (0..a.rows) |i| {
        const v = abs(T, a.get(i, j));
        if (std.math.isNan(v)) return v;
        x = @max(x, v);
    };
    return x;
}

/// The scale `dgesvd.f`, `dgeev.f` and `dgels.f` bring a matrix of
/// max-norm `anrm` to: `SMLNUM` below it, `BIGNUM` above, else none.
pub fn scaleTo(comptime F: type, anrm: F, smlnum: F) ?F {
    const bignum = 1 / smlnum;
    return if (anrm > 0 and anrm < smlnum) smlnum else if (anrm > bignum) bignum else null;
}

/// The drivers' prologue (`dgeev.f:341-359`, `dgees.f`, `dgesvd.f`): `A`'s
/// max-norm, `NanInA` or `InfInA` for a NaN or an infinity in it, then `A`
/// brought by `scaleTo` into `[SMLNUM, BIGNUM]`, `SMLNUM = sqrt(safmin) /
/// eps`. Undo it with `lascl(T, false, to, anrm, ...)`.
pub fn scaleIn(comptime T: type, a: Mat(T)) error{ NanInA, InfInA }!struct { anrm: Float(T), to: ?Float(T) } {
    const F = Float(T);
    const anrm = normMax(T, a);
    if (std.math.isNan(anrm)) return error.NanInA;
    if (std.math.isInf(anrm)) return error.InfInA;
    const smlnum = @sqrt(std.math.floatMin(F)) / std.math.floatEps(F);
    const to = scaleTo(F, anrm, smlnum);
    if (to) |t| lascl(T, false, anrm, t, a);
    return .{ .anrm = anrm, .to = to };
}

/// `DLARTG` (`dlartg.f90`): `[c s; -s c] [f; g] = [r; 0]`, `r` with `f`'s
/// sign.
pub fn lartg(comptime F: type, f: F, g: F) struct { c: F, s: F, r: F } {
    const safmin = std.math.floatMin(F);
    const safmax = 1 / safmin;
    const rtmin = @sqrt(safmin);
    const rtmax = @sqrt(safmax / 2);
    const f1 = @abs(f);
    const g1 = @abs(g);
    if (g == 0) return .{ .c = 1, .s = 0, .r = f };
    if (f == 0) return .{ .c = 0, .s = std.math.copysign(@as(F, 1), g), .r = g1 };
    if (f1 > rtmin and f1 < rtmax and g1 > rtmin and g1 < rtmax) {
        const d = @sqrt(f * f + g * g);
        const r = std.math.copysign(d, f);
        return .{ .c = f1 / d, .s = g / r, .r = r };
    }
    const u = @min(safmax, @max(safmin, f1, g1));
    const fs = f / u;
    const gs = g / u;
    const d = @sqrt(fs * fs + gs * gs);
    const r = std.math.copysign(d, f);
    return .{ .c = @abs(fs) / d, .s = gs / r, .r = r * u };
}

/// `?LARFG` (`dlarfg.f`, `zlarfg.f`): the reflector `H = I - tau v vᴴ`,
/// `v(0) = 1`, with `Hᴴ [alpha; x] = [beta; 0]`, `beta` real; `beta`
/// replaces `alpha`, `v`'s tail `x` (`n - 1`); returns `tau`. Complex
/// division is Smith's (`div`) where `zlarfg.f` takes `ZLADIV`.
pub fn larfg(io: Io, comptime T: type, n: usize, alpha: *T, x: [*]T, inc: usize) Error!T {
    const F = Float(T);
    const cx = comptime isComplex(T);
    if (n == 0 or (!cx and n == 1)) return zero(T);
    var xnorm = nrm2(io, T, n - 1, x, inc);
    var alphr = re(T, alpha.*);
    var alphi: F = if (cx) alpha.*.im else 0;
    if (xnorm == 0 and alphi == 0) return zero(T);
    var beta = -std.math.copysign(if (cx) lapy3(F, alphr, alphi, xnorm) else lapy2(F, alphr, xnorm), alphr);
    // `DLAMCH('S') / DLAMCH('E')`, `E` being half of `eps`
    const safmin = std.math.floatMin(F) / (std.math.floatEps(F) / 2);
    var knt: usize = 0;
    if (@abs(beta) < safmin) {
        const rsafmn = 1 / safmin;
        while (true) {
            knt += 1;
            try scal(io, T, n - 1, rsafmn, x, inc);
            beta *= rsafmn;
            alphi *= rsafmn;
            alphr *= rsafmn;
            if (!(@abs(beta) < safmin and knt < 20)) break;
        }
        xnorm = nrm2(io, T, n - 1, x, inc);
        beta = -std.math.copysign(if (cx) lapy3(F, alphr, alphi, xnorm) else lapy2(F, alphr, xnorm), alphr);
    }
    var tau: T = undefined;
    if (cx) {
        tau = .init((beta - alphr) / beta, -alphi / beta);
        const s = div(T, one(T), .init(alphr - beta, alphi));
        try scalComplex(io, T, n - 1, s, x, inc);
    } else {
        tau = (beta - alphr) / beta;
        try scal(io, T, n - 1, 1 / (alphr - beta), x, inc);
    }
    for (0..knt) |_| beta *= safmin;
    alpha.* = scalar(T, beta);
    return tau;
}

/// `DLAPY2` (`dlapy2.f`): `sqrt(x² + y²)` without overflow, NaN kept.
pub fn lapy2(comptime F: type, x: F, y: F) F {
    if (std.math.isNan(x)) return x;
    if (std.math.isNan(y)) return y;
    const w = @max(@abs(x), @abs(y));
    const z = @min(@abs(x), @abs(y));
    if (z == 0 or w > std.math.floatMax(F)) return w;
    return w * @sqrt(1 + (z / w) * (z / w));
}

/// `DLAPY3` (`dlapy3.f`): `sqrt(x² + y² + z²)` without overflow.
pub fn lapy3(comptime F: type, x: F, y: F, z: F) F {
    const xa = @abs(x);
    const ya = @abs(y);
    const za = @abs(z);
    const w = @max(xa, ya, za);
    if (w == 0 or w > std.math.floatMax(F)) return xa + ya + za;
    return w * @sqrt((xa / w) * (xa / w) + (ya / w) * (ya / w) + (za / w) * (za / w));
}

/// What `laset` returns.
pub const LasetError = error{ InvalidM, InvalidN, InvalidLda };

/// `?LASET` (`dlaset.f`, `zlaset.f`): the `m x n` `A`'s off-diagonal
/// entries set to `alpha` and its diagonal to `beta`; with `uplo`, only that
/// strict triangle and the diagonal, the rest untouched (null: all of it).
pub fn laset(comptime T: type, uplo: ?Uplo, m: isize, n: isize, alpha: T, beta: T, a: [*]T, lda: isize) LasetError!void {
    comptime requireFloat(T, "laset");
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, m)) return error.InvalidLda;
    if (m == 0 or n == 0) return;
    if (!util.fitsMatrix(m, n, lda, @sizeOf(T))) return error.InvalidLda;
    const rows: usize = @intCast(m);
    const la: usize = @intCast(lda);
    for (0..@intCast(n)) |j| {
        const from = if (uplo == .l) @min(j, rows) else 0;
        const to = if (uplo == .u) @min(j + 1, rows) else rows;
        for (from..to) |i| a[i + j * la] = if (i == j) beta else alpha;
    }
}
