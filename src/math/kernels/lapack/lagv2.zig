//! `DLAGV2`: the standardized generalized Schur form of a 2 x 2 pencil,
//! Reference-LAPACK 3.12.0 (`lapack-netlib/SRC/dlagv2.f`, OpenBLAS
//! v0.3.34-317; modified BSD), over `qz_2x2.zig`'s `xLAG2` and `xLASV2`.

const std = @import("std");
const auxiliary = @import("auxiliary.zig");
const qz_2x2 = @import("qz_2x2.zig");

const lartg = auxiliary.lartg;

pub fn Lagv2(comptime F: type) type {
    return struct {
        alphar: [2]F,
        alphai: [2]F,
        beta: [2]F,
        /// The left rotation `[csl snl; -snl csl]` and the right one.
        csl: F,
        snl: F,
        csr: F,
        snr: F,
    };
}

/// `DLAGV2` (`dlagv2.f`): rotations from the left and right that make the
/// 2 x 2 `(A, B)` (`B` upper triangular), in place, upper triangular for
/// real eigenvalues, else `B` diagonal with positive entries; the
/// eigenvalues `(alphar + i alphai) / beta`. `xLAG2` is faer's, with its
/// `safmin` (the reference passes `DLAMCH('S')`).
pub fn lagv2(comptime F: type, a: [*]F, lda: usize, b: [*]F, ldb: usize) Lagv2(F) {
    const safmin = std.math.floatMin(F);
    const ulp = std.math.floatEps(F);
    const A = struct {
        fn p(m: [*]F, ld: usize, i: usize, j: usize) *F {
            return &m[i + j * ld];
        }
    }.p;
    // `[x; y] := [c s; -s c] [x; y]` for rows (`ld` apart) or columns
    const rot2 = struct {
        fn f(x: [*]F, y: [*]F, inc: usize, c: F, s: F) void {
            for (0..2) |k| {
                const u = x[k * inc];
                const v = y[k * inc];
                x[k * inc] = c * u + s * v;
                y[k * inc] = c * v - s * u;
            }
        }
    }.f;
    const anorm = @max(@abs(A(a, lda, 0, 0).*) + @abs(A(a, lda, 1, 0).*), @abs(A(a, lda, 0, 1).*) + @abs(A(a, lda, 1, 1).*), safmin);
    const ascale = 1 / anorm;
    for ([_][2]usize{ .{ 0, 0 }, .{ 0, 1 }, .{ 1, 0 }, .{ 1, 1 } }) |ij| A(a, lda, ij[0], ij[1]).* *= ascale;
    const bnorm = @max(@abs(A(b, ldb, 0, 0).*), @abs(A(b, ldb, 0, 1).*) + @abs(A(b, ldb, 1, 1).*), safmin);
    const bscale = 1 / bnorm;
    for ([_][2]usize{ .{ 0, 0 }, .{ 0, 1 }, .{ 1, 1 } }) |ij| A(b, ldb, ij[0], ij[1]).* *= bscale;
    var r: Lagv2(F) = undefined;
    var wi: F = 0;
    var e: qz_2x2.Eig2(F) = undefined;
    if (@abs(A(a, lda, 1, 0).*) <= ulp) {
        r.csl = 1;
        r.snl = 0;
        r.csr = 1;
        r.snr = 0;
        A(a, lda, 1, 0).* = 0;
        A(b, ldb, 1, 0).* = 0;
    } else if (@abs(A(b, ldb, 0, 0).*) <= ulp) {
        const g = lartg(F, A(a, lda, 0, 0).*, A(a, lda, 1, 0).*);
        r.csl = g.c;
        r.snl = g.s;
        r.csr = 1;
        r.snr = 0;
        rot2(a, a + 1, lda, r.csl, r.snl);
        rot2(b, b + 1, ldb, r.csl, r.snl);
        A(a, lda, 1, 0).* = 0;
        A(b, ldb, 0, 0).* = 0;
        A(b, ldb, 1, 0).* = 0;
    } else if (@abs(A(b, ldb, 1, 1).*) <= ulp) {
        const g = lartg(F, A(a, lda, 1, 1).*, A(a, lda, 1, 0).*);
        r.csr = g.c;
        r.snr = -g.s;
        rot2(a, a + lda, 1, r.csr, r.snr);
        rot2(b, b + ldb, 1, r.csr, r.snr);
        r.csl = 1;
        r.snl = 0;
        A(a, lda, 1, 0).* = 0;
        A(b, ldb, 1, 0).* = 0;
        A(b, ldb, 1, 1).* = 0;
    } else {
        e = qz_2x2.eig2x2(F, .{ A(a, lda, 0, 0).*, A(a, lda, 0, 1).*, A(a, lda, 1, 0).*, A(a, lda, 1, 1).* }, .{ A(b, ldb, 0, 0).*, A(b, ldb, 0, 1).*, 0, A(b, ldb, 1, 1).* });
        wi = e.wi;
        if (wi == 0) {
            const s1 = e.scale1;
            var h1 = s1 * A(a, lda, 0, 0).* - e.wr1 * A(b, ldb, 0, 0).*;
            var h2 = s1 * A(a, lda, 0, 1).* - e.wr1 * A(b, ldb, 0, 1).*;
            const h3 = s1 * A(a, lda, 1, 1).* - e.wr1 * A(b, ldb, 1, 1).*;
            const rr = std.math.hypot(h1, h2);
            const qq = std.math.hypot(s1 * A(a, lda, 1, 0).*, h3);
            const g = if (rr > qq) lartg(F, h2, h1) else lartg(F, h3, s1 * A(a, lda, 1, 0).*);
            r.csr = g.c;
            r.snr = -g.s;
            rot2(a, a + lda, 1, r.csr, r.snr);
            rot2(b, b + ldb, 1, r.csr, r.snr);
            h1 = @max(@abs(A(a, lda, 0, 0).*) + @abs(A(a, lda, 0, 1).*), @abs(A(a, lda, 1, 0).*) + @abs(A(a, lda, 1, 1).*));
            h2 = @max(@abs(A(b, ldb, 0, 0).*) + @abs(A(b, ldb, 0, 1).*), @abs(A(b, ldb, 1, 0).*) + @abs(A(b, ldb, 1, 1).*));
            const gl = if (s1 * h1 >= @abs(e.wr1) * h2) lartg(F, A(b, ldb, 0, 0).*, A(b, ldb, 1, 0).*) else lartg(F, A(a, lda, 0, 0).*, A(a, lda, 1, 0).*);
            r.csl = gl.c;
            r.snl = gl.s;
            rot2(a, a + 1, lda, r.csl, r.snl);
            rot2(b, b + 1, ldb, r.csl, r.snl);
            A(a, lda, 1, 0).* = 0;
            A(b, ldb, 1, 0).* = 0;
        } else {
            const sv = qz_2x2.svdTriu2x2(F, A(b, ldb, 0, 0).*, A(b, ldb, 0, 1).*, A(b, ldb, 1, 1).*);
            r.snr = sv.snr;
            r.csr = sv.csr;
            r.snl = sv.snl;
            r.csl = sv.csl;
            rot2(a, a + 1, lda, r.csl, r.snl);
            rot2(b, b + 1, ldb, r.csl, r.snl);
            rot2(a, a + lda, 1, r.csr, r.snr);
            rot2(b, b + ldb, 1, r.csr, r.snr);
            A(b, ldb, 1, 0).* = 0;
            A(b, ldb, 0, 1).* = 0;
        }
    }
    for ([_][2]usize{ .{ 0, 0 }, .{ 1, 0 }, .{ 0, 1 }, .{ 1, 1 } }) |ij| {
        A(a, lda, ij[0], ij[1]).* *= anorm;
        A(b, ldb, ij[0], ij[1]).* *= bnorm;
    }
    if (wi == 0) {
        r.alphar = .{ A(a, lda, 0, 0).*, A(a, lda, 1, 1).* };
        r.alphai = .{ 0, 0 };
        r.beta = .{ A(b, ldb, 0, 0).*, A(b, ldb, 1, 1).* };
    } else {
        const ar = anorm * e.wr1 / e.scale1 / bnorm;
        const ai = anorm * wi / e.scale1 / bnorm;
        r.alphar = .{ ar, ar };
        r.alphai = .{ ai, -ai };
        r.beta = .{ 1, 1 };
    }
    return r;
}
