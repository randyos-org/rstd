//! `xTGSY2`, `tgsyl`'s kernel (`tgsyl.zig`), with `xGETC2`, `xGESC2` and `xLATDF`, Reference-LAPACK's
//! (`lapack-netlib/SRC/dtgsy2.f`, `dgetc2.f`, `dgesc2.f`, `dlatdf.f` and their forms, OpenBLAS v0.3.34-317; modified BSD).

const std = @import("std");
const mem = @import("../../../mem.zig");
const types = @import("../blas/types.zig");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");
const norm = @import("norm.zig");
const sylvester = @import("sylvester.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const scalar = ops.scalar;
const add = ops.add;
const sub = ops.sub;
const mul = ops.mul;
const SylvesterSolution = sylvester.SylvesterSolution;
const blockOf = sylvester.blockOf;

/// The `dif` estimate's running sum of squares (`RDSUM`, `RDSCAL`, as
/// `xLASSQ` keeps them) and the count of diagonal block pairs solved
/// (`PQ`).
pub fn Dif(comptime F: type) type {
    return struct { sum: F = 1, scale: F = 0, pq: usize = 0 };
}

/// `|x|`, the modulus for complex (Fortran's `ABS`).
fn absT(comptime T: type, x: T) Float(T) {
    return if (comptime isComplex(T)) std.math.hypot(x.re, x.im) else @abs(x);
}

fn conjT(comptime T: type, x: T) T {
    return if (comptime isComplex(T)) x.conjugate() else x;
}

/// `xTGSY2` (`lapack-netlib/SRC/dtgsy2.f`, `ztgsy2.f`), ReLAPACK's base case
/// for `tgsyl`: `A R - L B = scale C`, `D R - L E = scale F` (with `trans`,
/// `Aᴴ R + Dᴴ L = scale C`, `R Bᴴ + L Eᴴ = -scale F`) for `R` in `C` and `L`
/// in `F` (`m x n`), `A`, `B` upper quasi-triangular (complex: triangular),
/// `D`, `E` upper triangular. Each pair of diagonal blocks is a system of 2,
/// 4 or 8 unknowns, `Z x = b` (`Zᴴ` with `trans`) with `Z`'s Kronecker
/// form (`dtgsyl.f`'s (2)), by `getc2` and `gesc2`; with `ijob` 1 or 2 (not
/// `trans`), `latdf` instead adds the block's part of the `dif` estimate to
/// `dif`, unscaled. The reference's four real block shapes are one loop.
pub fn tgsy2(al: Allocator, io: Io, comptime T: type, trans: bool, ijob: u2, m: usize, n: usize, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, c: [*]T, ldc: usize, d: [*]const T, ldd: usize, e: [*]const T, lde: usize, f: [*]T, ldf: usize, dif: *Dif(Float(T))) ops.Error!SylvesterSolution(Float(T)) {
    const F = Float(T);
    var out: SylvesterSolution(F) = .{ .scale = 1, .perturbed = false };
    var p: usize = 0;
    var i: usize = 0;
    while (i < m) : (i += blockOf(T, m, a, lda, i)) p += 1;
    var q: usize = 0;
    i = 0;
    while (i < n) : (i += blockOf(T, n, b, ldb, i)) q += 1;
    dif.pq += p * q;
    const zero = scalar(T, 0);
    // rows by `A`'s blocks backward (forward with `trans`), columns by
    // `B`'s forward (backward)
    var jl: usize = if (trans) n else 0;
    while (if (trans) jl > 0 else jl < n) {
        const js: usize, const nb: usize = if (!trans) .{ jl, blockOf(T, n, b, ldb, jl) } else lastBlock(T, b, ldb, jl);
        var il: usize = if (trans) 0 else m;
        while (if (trans) il < m else il > 0) {
            const is: usize, const mb: usize = if (trans) .{ il, blockOf(T, m, a, lda, il) } else lastBlock(T, a, lda, il);
            const k = mb * nb;
            const zdim = 2 * k;
            // `Z = [ kron(I, A_ii)  -kron(B_jjᵀ, I) ; kron(I, D_ii)  -kron(E_jjᵀ, I) ]`
            var z: [64]T = @splat(zero);
            for (0..nb) |jj| for (0..mb) |ii| {
                const r = ii + mb * jj;
                for (0..mb) |kk| {
                    z[r + 8 * (kk + mb * jj)] = a[is + ii + (is + kk) * lda];
                    if (kk >= ii) z[k + r + 8 * (kk + mb * jj)] = d[is + ii + (is + kk) * ldd];
                }
                for (0..nb) |kk| {
                    z[r + 8 * (k + ii + mb * kk)] = ops.neg(T, b[js + kk + (js + jj) * ldb]);
                    if (kk <= jj) z[k + r + 8 * (k + ii + mb * kk)] = ops.neg(T, e[js + kk + (js + jj) * lde]);
                }
            };
            if (trans) for (0..zdim) |col| for (col..zdim) |row| {
                const t = z[row + 8 * col];
                z[row + 8 * col] = conjT(T, z[col + 8 * row]);
                z[col + 8 * row] = conjT(T, t);
            };
            var rhs: [8]T = undefined;
            for (0..nb) |jj| for (0..mb) |ii| {
                rhs[ii + mb * jj] = c[is + ii + (js + jj) * ldc];
                rhs[k + ii + mb * jj] = f[is + ii + (js + jj) * ldf];
            };
            var ipiv: [8]usize = undefined;
            var jpiv: [8]usize = undefined;
            if (getc2(T, zdim, &z, &ipiv, &jpiv) != 0) out.perturbed = true;
            if (!trans and ijob != 0) {
                try latdf(al, io, T, ijob, zdim, &z, &rhs, dif, &ipiv, &jpiv);
            } else {
                const scaloc = gesc2(T, zdim, &z, &rhs, &ipiv, &jpiv);
                if (scaloc != 1) {
                    for (0..n) |col| for (0..m) |row| {
                        c[row + col * ldc] = ops.mulReal(T, c[row + col * ldc], scaloc);
                        f[row + col * ldf] = ops.mulReal(T, f[row + col * ldf], scaloc);
                    };
                    out.scale *= scaloc;
                }
            }
            for (0..nb) |jj| for (0..mb) |ii| {
                c[is + ii + (js + jj) * ldc] = rhs[ii + mb * jj];
                f[is + ii + (js + jj) * ldf] = rhs[k + ii + mb * jj];
            };
            // R(i, j) and L(i, j) into the rest
            const ie = is + mb;
            const je = js + nb;
            if (!trans) {
                for (0..nb) |jj| for (0..is) |row| {
                    var sc = zero;
                    var sf = zero;
                    for (0..mb) |kk| {
                        sc = add(T, sc, mul(T, a[row + (is + kk) * lda], rhs[kk + mb * jj]));
                        sf = add(T, sf, mul(T, d[row + (is + kk) * ldd], rhs[kk + mb * jj]));
                    }
                    c[row + (js + jj) * ldc] = sub(T, c[row + (js + jj) * ldc], sc);
                    f[row + (js + jj) * ldf] = sub(T, f[row + (js + jj) * ldf], sf);
                };
                for (je..n) |col| for (0..mb) |ii| {
                    var sc = zero;
                    var sf = zero;
                    for (0..nb) |kk| {
                        const l = rhs[k + ii + mb * kk];
                        sc = add(T, sc, mul(T, l, b[js + kk + col * ldb]));
                        sf = add(T, sf, mul(T, l, e[js + kk + col * lde]));
                    }
                    c[is + ii + col * ldc] = add(T, c[is + ii + col * ldc], sc);
                    f[is + ii + col * ldf] = add(T, f[is + ii + col * ldf], sf);
                };
            } else {
                for (0..js) |col| for (0..mb) |ii| {
                    var s = zero;
                    for (0..nb) |kk| {
                        s = add(T, s, mul(T, rhs[ii + mb * kk], conjT(T, b[col + (js + kk) * ldb])));
                        s = add(T, s, mul(T, rhs[k + ii + mb * kk], conjT(T, e[col + (js + kk) * lde])));
                    }
                    f[is + ii + col * ldf] = add(T, f[is + ii + col * ldf], s);
                };
                for (0..nb) |jj| for (ie..m) |row| {
                    var s = zero;
                    for (0..mb) |kk| {
                        s = add(T, s, mul(T, conjT(T, a[is + kk + row * lda]), rhs[kk + mb * jj]));
                        s = add(T, s, mul(T, conjT(T, d[is + kk + row * ldd]), rhs[k + kk + mb * jj]));
                    }
                    c[row + (js + jj) * ldc] = sub(T, c[row + (js + jj) * ldc], s);
                };
            }
            il = if (trans) ie else is;
        }
        jl = if (trans) js else js + nb;
    }
    return out;
}

/// The diagonal block of the quasi-triangular `a` that ends before `end`:
/// its first index and size.
fn lastBlock(comptime T: type, a: [*]const T, lda: usize, end: usize) [2]usize {
    const s: usize = if (comptime !isComplex(T)) (if (end >= 2 and a[end - 1 + (end - 2) * lda] != 0) 2 else 1) else 1;
    return .{ end - s, s };
}

/// `xGETC2` (`dgetc2.f`, `zgetc2.f`): `z` (`n x n`, leading dimension 8) to
/// `P L U Q` by complete pivoting, pivots below `SMIN` raised to it; the
/// last such pivot's index (from 1), or 0.
fn getc2(comptime T: type, n: usize, z: *[64]T, ipiv: *[8]usize, jpiv: *[8]usize) usize {
    const F = Float(T);
    const eps = std.math.floatEps(F);
    const smlnum = std.math.floatMin(F) / eps;
    var info: usize = 0;
    if (n == 1) {
        ipiv[0] = 0;
        jpiv[0] = 0;
        if (absT(T, z[0]) < smlnum) {
            info = 1;
            z[0] = scalar(T, smlnum);
        }
        return info;
    }
    var smin: F = 0;
    for (0..n - 1) |i| {
        var xmax: F = 0;
        var ipv = i;
        var jpv = i;
        for (i..n) |jp| for (i..n) |ip| if (absT(T, z[ip + 8 * jp]) >= xmax) {
            xmax = absT(T, z[ip + 8 * jp]);
            ipv = ip;
            jpv = jp;
        };
        if (i == 0) smin = @max(eps * xmax, smlnum);
        if (ipv != i) for (0..n) |col| mem.swap(T, &z[ipv + 8 * col], &z[i + 8 * col]);
        ipiv[i] = ipv;
        if (jpv != i) for (0..n) |row| mem.swap(T, &z[row + 8 * jpv], &z[row + 8 * i]);
        jpiv[i] = jpv;
        if (absT(T, z[i + 8 * i]) < smin) {
            info = i + 1;
            z[i + 8 * i] = scalar(T, smin);
        }
        for (i + 1..n) |j| z[j + 8 * i] = ops.div(T, z[j + 8 * i], z[i + 8 * i]);
        for (i + 1..n) |col| for (i + 1..n) |row| {
            z[row + 8 * col] = sub(T, z[row + 8 * col], mul(T, z[row + 8 * i], z[i + 8 * col]));
        };
    }
    if (absT(T, z[n - 1 + 8 * (n - 1)]) < smin) {
        info = n;
        z[n - 1 + 8 * (n - 1)] = scalar(T, smin);
    }
    ipiv[n - 1] = n - 1;
    jpiv[n - 1] = n - 1;
    return info;
}

/// `xGESC2` (`dgesc2.f`, `zgesc2.f`): `rhs` to `getc2`'s `z`'s solution
/// times the returned scale, `<= 1`, chosen against overflow.
fn gesc2(comptime T: type, n: usize, z: *const [64]T, rhs: *[8]T, ipiv: *const [8]usize, jpiv: *const [8]usize) Float(T) {
    const F = Float(T);
    const smlnum = std.math.floatMin(F) / std.math.floatEps(F);
    for (0..n - 1) |i| mem.swap(T, &rhs[i], &rhs[ipiv[i]]);
    for (0..n - 1) |i| for (i + 1..n) |j| {
        rhs[j] = sub(T, rhs[j], mul(T, z[j + 8 * i], rhs[i]));
    };
    var scale: F = 1;
    // `I?AMAX`: `|re| + |im|` for complex
    var big: usize = 0;
    for (0..n) |i| if (ops.abs1(T, rhs[i]) > ops.abs1(T, rhs[big])) {
        big = i;
    };
    if (2 * smlnum * absT(T, rhs[big]) > absT(T, z[n - 1 + 8 * (n - 1)])) {
        const temp = 0.5 / absT(T, rhs[big]);
        for (rhs[0..n]) |*x| x.* = ops.mulReal(T, x.*, temp);
        scale *= temp;
    }
    backSolve(T, n, z, rhs);
    unswap(T, n, rhs, jpiv);
    return scale;
}

/// `xLASWP` with `INCX = -1`: the swaps undone, last first.
fn unswap(comptime T: type, n: usize, x: *[8]T, piv: *const [8]usize) void {
    var i = n - 1;
    while (i > 0) {
        i -= 1;
        mem.swap(T, &x[i], &x[piv[i]]);
    }
}

/// `U x = rhs` in place, as `xGESC2` divides: `rhs(i) *= 1 / U(i, i)`,
/// then `rhs(i) -= rhs(j) (U(i, j) / U(i, i))`.
fn backSolve(comptime T: type, n: usize, z: *const [64]T, rhs: *[8]T) void {
    var i = n;
    while (i > 0) {
        i -= 1;
        const temp = ops.div(T, scalar(T, 1), z[i + 8 * i]);
        rhs[i] = mul(T, rhs[i], temp);
        for (i + 1..n) |j| rhs[i] = sub(T, rhs[i], mul(T, rhs[j], mul(T, z[i + 8 * j], temp)));
    }
}

/// The real part of `conj(x)ᵀ y` over `n` entries, `?DOT`'s or `ZDOTC`'s.
fn redot(comptime T: type, n: usize, x: [*]const T, y: [*]const T) Float(T) {
    var s: Float(T) = 0;
    for (0..n) |i| s += if (comptime isComplex(T)) x[i].re * y[i].re + x[i].im * y[i].im else x[i] * y[i];
    return s;
}

/// `xLATDF` (`dlatdf.f`, `zlatdf.f`): `Z x = b` for `getc2`'s `z` with `b`
/// chosen, from `rhs`, to make `x` large, then `x`'s sum of squares into
/// `dif`. `ijob` 1: `b`'s entries `±1` by look-ahead; 2: `rhs` plus or
/// minus an approximate null vector of `Z` (`xGECON`'s last `lacn2` vector).
fn latdf(al: Allocator, io: Io, comptime T: type, ijob: u2, n: usize, z: *const [64]T, rhs: *[8]T, dif: *Dif(Float(T)), ipiv: *const [8]usize, jpiv: *const [8]usize) ops.Error!void {
    const F = Float(T);
    const one = scalar(T, 1);
    if (ijob == 2) {
        var xm: [8]T = @splat(scalar(T, 0));
        try nullVector(al, io, T, n, z, &xm);
        unswap(T, n, &xm, ipiv);
        const temp = 1 / @sqrt(redot(T, n, &xm, &xm));
        var xp: [8]T = undefined;
        for (0..n) |k| {
            xm[k] = ops.mulReal(T, xm[k], temp);
            xp[k] = add(T, xm[k], rhs[k]);
            rhs[k] = sub(T, rhs[k], xm[k]);
        }
        _ = gesc2(T, n, z, rhs, ipiv, jpiv);
        _ = gesc2(T, n, z, &xp, ipiv, jpiv);
        var sp: F = 0;
        var sr: F = 0;
        for (0..n) |k| {
            sp += ops.abs1(T, xp[k]);
            sr += ops.abs1(T, rhs[k]);
        }
        if (sp > sr) rhs.* = xp;
        norm.lassq(T, .of(rhs[0..n]), &dif.scale, &dif.sum);
        return;
    }
    for (0..n - 1) |i| mem.swap(T, &rhs[i], &rhs[ipiv[i]]);
    // the L part, each `b(j)` +1 or -1 by look-ahead
    var pmone: F = -1;
    for (0..n - 1) |j| {
        const col = z[8 * j ..].ptr;
        const bp = add(T, rhs[j], one);
        const bm = sub(T, rhs[j], one);
        var splus: F = 1 + redot(T, n - j - 1, col + j + 1, col + j + 1);
        const sminu = redot(T, n - j - 1, col + j + 1, rhs[j + 1 ..].ptr);
        splus *= if (comptime isComplex(T)) rhs[j].re else rhs[j];
        if (splus > sminu) {
            rhs[j] = bp;
        } else if (sminu > splus) {
            rhs[j] = bm;
        } else {
            // equal sums: -1 the first time, +1 after
            rhs[j] = add(T, rhs[j], scalar(T, pmone));
            pmone = 1;
        }
        const temp = ops.neg(T, rhs[j]);
        for (j + 1..n) |k| rhs[k] = add(T, rhs[k], mul(T, temp, col[k]));
    }
    // the U part, with look-ahead for b(n)
    var work = rhs.*;
    work[n - 1] = add(T, rhs[n - 1], one);
    rhs[n - 1] = sub(T, rhs[n - 1], one);
    var splus: F = 0;
    var sminu: F = 0;
    var i = n;
    while (i > 0) {
        i -= 1;
        const temp = ops.div(T, one, z[i + 8 * i]);
        work[i] = mul(T, work[i], temp);
        rhs[i] = mul(T, rhs[i], temp);
        for (i + 1..n) |k| {
            work[i] = sub(T, work[i], mul(T, work[k], mul(T, z[i + 8 * k], temp)));
            rhs[i] = sub(T, rhs[i], mul(T, rhs[k], mul(T, z[i + 8 * k], temp)));
        }
        splus += absT(T, work[i]);
        sminu += absT(T, rhs[i]);
    }
    if (splus > sminu) rhs.* = work;
    unswap(T, n, rhs, jpiv);
    norm.lassq(T, .of(rhs[0..n]), &dif.scale, &dif.sum);
}

/// `xGECON('I', N, Z, LDZ, ONE, ...)` as `xLATDF` uses it: the vector
/// `lacn2` leaves in `xGECON`'s `WORK(N+1:2N)`, `(L U)⁻¹ w` for the `w` it
/// found, from `getc2`'s factors (its swaps left out, as there). An
/// overflowing solve ends it early, as `xGECON`'s `GO TO 20`.
fn nullVector(al: Allocator, io: Io, comptime T: type, n: usize, z: *const [64]T, v: *[8]T) ops.Error!void {
    const F = Float(T);
    const Inv = struct {
        al: Allocator,
        io: Io,
        n: usize,
        z: [*]const T,
        cnorm: [16]F = undefined,
        normin: bool = false,

        const Self = @This();
        const Stop = ops.Error || error{Overflow};

        // the infinity norm's `KASE1 = 2`: `KASE = 2` is `U⁻¹ L⁻¹`
        fn apply(s: *Self, adjoint: bool, x: [*]T) Stop!void {
            const op: types.Trans = if (comptime isComplex(T)) .c else .t;
            const lower: [*]F = &s.cnorm;
            const upper: [*]F = s.cnorm[8..].ptr;
            var sl: F = undefined;
            var su: F = undefined;
            if (adjoint) {
                sl = try norm.scaledSolve(s.al, s.io, T, false, .n, true, s.normin, s.n, s.z, 8, x, lower);
                su = try norm.scaledSolve(s.al, s.io, T, true, .n, false, s.normin, s.n, s.z, 8, x, upper);
            } else {
                su = try norm.scaledSolve(s.al, s.io, T, true, op, false, s.normin, s.n, s.z, 8, x, upper);
                sl = try norm.scaledSolve(s.al, s.io, T, false, op, true, s.normin, s.n, s.z, 8, x, lower);
            }
            const scale = sl * su;
            s.normin = true;
            if (scale != 1) {
                const ix = blas_calls.iamax(T, s.n, x, 1);
                if (scale < ops.abs1(T, x[ix]) * std.math.floatMin(F) or scale == 0) return error.Overflow;
                try norm.rscl(s.io, T, scale, .of(x[0..s.n]));
            }
        }
    };
    var ctx: Inv = .{ .al = al, .io = io, .n = n, .z = z };
    var x: [8]T = undefined;
    _ = norm.lacn2(al, T, Inv.Stop, n, v, &x, &ctx, Inv.apply) catch |e| switch (e) {
        error.Overflow => {},
        error.OutOfMemory, error.Canceled => |y| return y,
        // `n <= 8`
        error.SliceTooShort => unreachable,
    };
}
