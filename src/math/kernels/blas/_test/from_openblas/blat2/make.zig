//! `SMAKE`, `LSERES` and `SMVCH` of `test/[sdcz]blat2.f`: operands, the untouched-part check and the reference product.

const std = @import("std");
const ob = @import("../__root__.zig");

const Float = ob.Float;
const isComplex = ob.isComplex;

pub const Kind = enum { ge, gb, sy, sb, sp, tr, tb, tp };

pub fn rogue(comptime T: type) T {
    return ob.val(T, -1e10, 1e10);
}

/// `SMAKE`/`CMAKE`: an `m` by `n` matrix within the band `kl`, `ku`, into
/// `a` (leading dimension `ld`) and stored as `kind` in `aa`, unused
/// elements rogue. `sy`, `sb`, `sp` are Hermitian for complex `T`.
pub fn make(comptime T: type, kind: Kind, uplo: u8, diag: u8, m: usize, n: usize, a: []T, ld: usize, aa: []T, lda: usize, kl: isize, ku: isize, gen: *ob.Gen, transl: f64) void {
    const gen_ = kind == .ge or kind == .gb;
    const sym = kind == .sy or kind == .sb or kind == .sp;
    const tri_ = kind == .tr or kind == .tb or kind == .tp;
    const up = (sym or tri_) and uplo == 'U';
    const lo = (sym or tri_) and uplo == 'L';
    const unit = tri_ and diag == 'U';
    const herm = sym and comptime isComplex(T);
    const z = ob.zero(T);
    const im: isize = @intCast(m);
    const in: isize = @intCast(n);
    const A = struct {
        fn at(buf: []T, l: usize, i: isize, j: isize) *T {
            return &buf[@as(usize, @intCast(i - 1)) + @as(usize, @intCast(j - 1)) * l];
        }
    }.at;
    var j: isize = 1;
    while (j <= in) : (j += 1) {
        var i: isize = 1;
        while (i <= im) : (i += 1) {
            if (gen_ or (up and i <= j) or (lo and i >= j)) {
                A(a, ld, i, j).* = if ((i <= j and j - i <= ku) or (i >= j and i - j <= kl))
                    ob.add(gen.next(T), ob.val(T, transl, 0))
                else
                    z;
                if (i != j) {
                    if (sym) {
                        A(a, ld, j, i).* = ob.conj(A(a, ld, i, j).*);
                    } else if (tri_) A(a, ld, j, i).* = z;
                }
            }
        }
        if (herm) A(a, ld, j, j).* = ob.val(T, re(A(a, ld, j, j).*), 0);
        if (tri_) A(a, ld, j, j).* = ob.add(A(a, ld, j, j).*, ob.val(T, 1, 0));
        if (unit) A(a, ld, j, j).* = ob.val(T, 1, 0);
    }
    const ilda: isize = @intCast(lda);
    const AA = struct {
        fn at(buf: []T, k: isize) *T {
            return &buf[@intCast(k - 1)];
        }
    }.at;
    switch (kind) {
        .ge => {
            j = 1;
            while (j <= in) : (j += 1) {
                var i: isize = 1;
                while (i <= ilda) : (i += 1) AA(aa, i + (j - 1) * ilda).* = if (i <= im) A(a, ld, i, j).* else rogue(T);
            }
        },
        .gb => {
            j = 1;
            while (j <= in) : (j += 1) {
                const first = @max(1, ku + 2 - j);
                const last = @min(kl + ku + 1, ku + 1 + im - j);
                var i: isize = 1;
                while (i <= ilda) : (i += 1) {
                    AA(aa, i + (j - 1) * ilda).* = if (i >= first and i <= last) A(a, ld, i + j - ku - 1, j).* else rogue(T);
                }
            }
        },
        .sy, .tr => {
            j = 1;
            while (j <= in) : (j += 1) {
                const ibeg: isize = if (up) 1 else if (unit) j + 1 else j;
                const iend: isize = if (up) (if (unit) j - 1 else j) else in;
                var i: isize = 1;
                while (i <= ilda) : (i += 1) {
                    AA(aa, i + (j - 1) * ilda).* = if (i >= ibeg and i <= iend) A(a, ld, i, j).* else rogue(T);
                }
                if (herm) setIm(T, AA(aa, j + (j - 1) * ilda), -1e10);
            }
        },
        .sb, .tb => {
            j = 1;
            while (j <= in) : (j += 1) {
                const kk: isize = if (up) kl + 1 else 1;
                const ibeg: isize = if (up) @max(1, kl + 2 - j) else if (unit) 2 else 1;
                const iend: isize = if (up) (if (unit) kl else kl + 1) else @min(kl + 1, 1 + im - j);
                var i: isize = 1;
                while (i <= ilda) : (i += 1) {
                    AA(aa, i + (j - 1) * ilda).* = if (i >= ibeg and i <= iend) A(a, ld, i + j - kk, j).* else rogue(T);
                }
                if (herm) setIm(T, AA(aa, kk + (j - 1) * ilda), -1e10);
            }
        },
        .sp, .tp => {
            var ioff: isize = 0;
            j = 1;
            while (j <= in) : (j += 1) {
                const ibeg: isize = if (up) 1 else j;
                const iend: isize = if (up) j else in;
                var i = ibeg;
                while (i <= iend) : (i += 1) {
                    ioff += 1;
                    AA(aa, ioff).* = A(a, ld, i, j).*;
                    if (i == j) {
                        if (unit) AA(aa, ioff).* = rogue(T);
                        if (herm) setIm(T, AA(aa, ioff), -1e10);
                    }
                }
            }
        },
    }
}

pub fn re(v: anytype) Float(@TypeOf(v)) {
    return if (comptime isComplex(@TypeOf(v))) v.re else v;
}

pub fn setIm(comptime T: type, p: *T, v: Float(T)) void {
    if (comptime isComplex(T)) p.im = v;
}

/// `SMAKE('GE', ' ', ' ', 1, n, x, 1, xx, inc, ...)`: a vector, and its
/// storage at stride `inc` with rogue gaps.
pub fn makeVec(comptime T: type, n: usize, x: []T, xx: []T, inc: usize, gen: *ob.Gen, transl: f64) void {
    make(T, .ge, ' ', ' ', 1, n, x, 1, xx, inc, 0, @as(isize, @intCast(n)) - 1, gen, transl);
}

/// `makeVec`, and the middle element zero in both, as the Fortran sets it
/// for `n > 1`.
pub fn makeVecZeroMid(comptime T: type, n: usize, x: []T, xx: []T, inc: usize, gen: *ob.Gen, transl: f64) void {
    makeVec(T, n, x, xx, inc, gen, transl);
    if (n > 1) {
        x[n / 2 - 1] = ob.zero(T);
        xx[inc * (n / 2 - 1)] = ob.zero(T);
    }
}

/// `LSERES`: the elements of the `m` by `n` storage outside the result agree,
/// rows past `m` for `ge`, the other triangle for `sy`.
pub fn lseres(comptime T: type, kind: Kind, uplo: u8, m: usize, n: usize, aa: []const T, as: []const T, lda: usize) bool {
    for (0..n) |j| for (0..lda) |i| {
        const inside = switch (kind) {
            .ge => i < m,
            .sy => if (uplo == 'U') i <= j else i >= j and i < n,
            else => true,
        };
        if (!inside and !ob.eql(aa[i + j * lda], as[i + j * lda])) return false;
    };
    return true;
}

// ---- SMVCH ----------------------------------------------------------------

pub const Result = ob.Result;

/// `SMVCH`/`CMVCH`: `alpha*op(A)*x + beta*y` from the dense `a` (leading
/// dimension `ld`) and compacted `x`, `y` into `yt`, gauges into `g`, and the
/// test ratio of `yy` (at stride `|incy|`) against it.
pub fn mvch(comptime T: type, trans: u8, m: usize, n: usize, alpha: T, a: []const T, ld: usize, x: []const T, incx: isize, beta: T, y: []const T, incy: isize, yt: []ob.Ref(T), g: []Float(ob.Ref(T)), yy: []const T) Result {
    const R = ob.Ref(T);
    const RF = Float(R);
    const tran = trans == 'T' or trans == 'C';
    const ctran = trans == 'C';
    const ml = if (tran) n else m;
    const nl = if (tran) m else n;
    const ra = ob.cast(R, alpha);
    const rb = ob.cast(R, beta);
    for (0..ml) |i| {
        const iy = if (incy < 0) ml - 1 - i else i;
        var acc = ob.zero(R);
        var gauge: RF = 0;
        for (0..nl) |j| {
            const jx = if (incx < 0) nl - 1 - j else j;
            var e = ob.cast(R, if (tran) a[j + i * ld] else a[i + j * ld]);
            if (ctran) e = ob.conj(e);
            const xv = ob.cast(R, x[jx]);
            acc = ob.add(acc, ob.mul(e, xv));
            gauge += ob.abs1(e) * ob.abs1(xv);
        }
        const yv = ob.cast(R, y[iy]);
        yt[iy] = ob.add(ob.mul(ra, acc), ob.mul(rb, yv));
        g[iy] = ob.abs1(ra) * gauge + ob.abs1(rb) * ob.abs1(yv);
    }
    const eps: RF = std.math.floatEps(Float(T));
    var err: RF = 0;
    const ay: usize = @abs(incy);
    for (0..ml) |i| {
        var erri = ob.absv(ob.sub(yt[i], ob.cast(R, yy[i * ay]))) / eps;
        if (g[i] != 0) erri /= g[i];
        err = @max(err, erri);
        if (err * @sqrt(eps) >= 1) {
            ob.say("row {d}: expected {e}, computed {e}\n", .{ i + 1, re(yt[i]), re(yy[i * ay]) });
            return .{ .err = @floatCast(err), .fatal = true };
        }
    }
    return .{ .err = @floatCast(err), .fatal = false };
}
