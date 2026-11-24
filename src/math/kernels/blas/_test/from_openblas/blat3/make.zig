//! `SMAKE`, `LSERES` and `SMMCH` of `test/[sdcz]blat3.f`: operands, the untouched-part check and the reference product.

const std = @import("std");
const ob = @import("../__root__.zig");

const Float = ob.Float;
const isComplex = ob.isComplex;

// ---- SMAKE, LSERES, SMMCH ------------------------------------------------------

pub const Kind = enum { ge, sy, he, tr };

pub fn rogue(comptime T: type) T {
    return ob.val(T, -1e10, 1e10);
}

pub fn re(v: anytype) Float(@TypeOf(v)) {
    return if (comptime isComplex(@TypeOf(v))) v.re else v;
}

/// `SMAKE`/`CMAKE` (level 3): an `m` by `n` matrix into `a` (leading
/// dimension `ld`) and stored as `kind` in `aa`, unused elements rogue; off
/// the diagonal, column `n/2` is zero when `n > 3`.
pub fn make(comptime T: type, kind: Kind, uplo: u8, diag: u8, m: usize, n: usize, a: []T, ld: usize, aa: []T, lda: usize, gen: *ob.Gen) void {
    const gen_ = kind == .ge;
    const her = kind == .he;
    const sym = kind == .sy;
    const tri = kind == .tr;
    const up = (her or sym or tri) and uplo == 'U';
    const lo = (her or sym or tri) and uplo == 'L';
    const unit = tri and diag == 'U';
    const z = ob.zero(T);
    for (1..n + 1) |j| {
        for (1..m + 1) |i| {
            if (gen_ or (up and i <= j) or (lo and i >= j)) {
                a[(i - 1) + (j - 1) * ld] = gen.next(T);
                if (i != j) {
                    if (n > 3 and j == n / 2) a[(i - 1) + (j - 1) * ld] = z;
                    const v = a[(i - 1) + (j - 1) * ld];
                    if (her) {
                        a[(j - 1) + (i - 1) * ld] = ob.conj(v);
                    } else if (sym) {
                        a[(j - 1) + (i - 1) * ld] = v;
                    } else if (tri) a[(j - 1) + (i - 1) * ld] = z;
                }
            }
        }
        const d = &a[(j - 1) + (j - 1) * ld];
        if (her) d.* = ob.val(T, re(d.*), 0);
        if (tri) d.* = ob.add(d.*, ob.val(T, 1, 0));
        if (unit) d.* = ob.val(T, 1, 0);
    }
    for (1..n + 1) |j| {
        const ibeg: usize = if (gen_ or up) 1 else if (unit) j + 1 else j;
        const iend: usize = if (gen_) m else if (up) (if (unit) j - 1 else j) else n;
        for (1..lda + 1) |i| {
            aa[(i - 1) + (j - 1) * lda] = if (i >= ibeg and i <= iend) a[(i - 1) + (j - 1) * ld] else rogue(T);
        }
        if (her) {
            const p = &aa[(j - 1) + (j - 1) * lda];
            if (comptime isComplex(T)) p.im = -1e10;
        }
    }
}

/// `LSERES`: the storage outside the result is unchanged: rows past `m` for
/// `ge`, the other triangle for `sy` and `he`.
pub fn lseres(comptime T: type, kind: Kind, uplo: u8, m: usize, n: usize, aa: []const T, as: []const T, lda: usize) bool {
    for (0..n) |j| for (0..lda) |i| {
        const inside = switch (kind) {
            .ge => i < m,
            else => if (uplo == 'U') i <= j else i >= j and i < n,
        };
        if (!inside and !ob.eql(aa[i + j * lda], as[i + j * lda])) return false;
    };
    return true;
}

pub const Result = ob.Result;

/// `SMMCH`/`CMMCH`: `alpha*op(A)*op(B) + beta*C` for `m` by `n` and `kk`
/// from the dense `a`, `b`, `c`, against `cc`; gauges are sums of `ABS1`.
pub fn mmch(comptime T: type, transa: u8, transb: u8, m: usize, n: usize, kk: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []const T, ldc: usize, ct: []ob.Ref(T), g: []Float(ob.Ref(T)), cc: []const T, ldcc: usize) Result {
    return mmchLoose(T, .none, transa, transb, m, n, kk, alpha, a, lda, b, ldb, beta, c, ldc, ct, g, cc, ldcc);
}

/// Which part of a result is held to three times the bound.
pub const Loose = enum { none, re, im };

/// `mmch` for `gemm3m`: the part it forms from three products, a difference
/// of them, gauged at three times the bound (Higham, SIAM J. Matrix Anal.
/// Appl. 13(3), 1992): the imaginary one, or the real one where exactly one
/// operand is conjugated.
pub fn mmch3m(comptime T: type, transa: u8, transb: u8, m: usize, n: usize, kk: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []const T, ldc: usize, ct: []ob.Ref(T), g: []Float(ob.Ref(T)), cc: []const T, ldcc: usize) Result {
    const loose: Loose = if ((transa == 'C') == (transb == 'C')) .im else .re;
    return mmchLoose(T, loose, transa, transb, m, n, kk, alpha, a, lda, b, ldb, beta, c, ldc, ct, g, cc, ldcc);
}

/// An error `d` with its `loose` part a third.
pub fn loosen(d: anytype, loose: Loose) @TypeOf(d) {
    if (comptime !isComplex(@TypeOf(d))) return d;
    var v = d;
    switch (loose) {
        .none => {},
        .re => v.re /= 3,
        .im => v.im /= 3,
    }
    return v;
}

pub fn mmchLoose(comptime T: type, loose: Loose, transa: u8, transb: u8, m: usize, n: usize, kk: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []const T, ldc: usize, ct: []ob.Ref(T), g: []Float(ob.Ref(T)), cc: []const T, ldcc: usize) Result {
    const R = ob.Ref(T);
    const RF = Float(R);
    const tra = transa == 'T' or transa == 'C';
    const trb = transb == 'T' or transb == 'C';
    const ra = ob.cast(R, alpha);
    const rb = ob.cast(R, beta);
    const eps: RF = std.math.floatEps(Float(T));
    var err: RF = 0;
    for (0..n) |j| {
        for (0..m) |i| {
            var acc = ob.zero(R);
            var gauge: RF = 0;
            for (0..kk) |k| {
                var x = ob.cast(R, if (tra) a[k + i * lda] else a[i + k * lda]);
                var y = ob.cast(R, if (trb) b[j + k * ldb] else b[k + j * ldb]);
                if (transa == 'C') x = ob.conj(x);
                if (transb == 'C') y = ob.conj(y);
                acc = ob.add(acc, ob.mul(x, y));
                gauge += ob.abs1(x) * ob.abs1(y);
            }
            const cv = ob.cast(R, c[i + j * ldc]);
            ct[i] = ob.add(ob.mul(ra, acc), ob.mul(rb, cv));
            g[i] = ob.abs1(ra) * gauge + ob.abs1(rb) * ob.abs1(cv);
        }
        for (0..m) |i| {
            var erri = ob.abs1(loosen(ob.sub(ct[i], ob.cast(R, cc[i + j * ldcc])), loose)) / eps;
            if (g[i] != 0) erri /= g[i];
            err = @max(err, erri);
            if (err * @sqrt(eps) >= 1) {
                ob.say("({d},{d}): expected {e}, computed {e}\n", .{ i + 1, j + 1, re(ct[i]), re(cc[i + j * ldcc]) });
                return .{ .err = @floatCast(err), .fatal = true };
            }
        }
    }
    return .{ .err = @floatCast(err), .fatal = false };
}
