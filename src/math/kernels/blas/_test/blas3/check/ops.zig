//! The reference: each routine's operators as dense matrices in the wider type.

const std = @import("std");
const types = @import("../../../types.zig");
const reference = @import("../../reference.zig");
const c3 = @import("../check.zig");
const c3one = @import("one.zig");

const isComplex = types.isComplex;
const testing = std.testing;
const Case = c3.Case;
const Ref = c3.Ref;
const RefFloat = c3.RefFloat;
const Trans = c3.Trans;
const conjIf = c3.conjIf;
const stored = c3.stored;
const toRef = c3.toRef;
const zeroW = c3.zeroW;
const compare = c3one.compare;
const fail = c3one.fail;

/// The reference: the routine's operators as dense matrices, by element.
pub fn Ops(comptime T: type) type {
    const W = Ref(T);
    return struct {
        const Self = @This();
        case: Case,
        a: []const T,
        lda: usize,
        b: []const T,
        ldb: usize,
        c0: []const T,
        ldc: usize,

        pub fn at(x: []const T, ld: usize, i: usize, j: usize) W {
            return toRef(T, x[i + j * ld]);
        }

        /// `op(X)(i, l)` of a general `X`.
        pub fn general(x: []const T, ld: usize, t: Trans, i: usize, l: usize) W {
            const v = if (t.transposed()) at(x, ld, l, i) else at(x, ld, i, l);
            return conjIf(v, t.conjugated());
        }

        /// `A`'s full symmetric or Hermitian matrix at `(i, j)`.
        pub fn sym(s: Self, i: usize, j: usize) W {
            const c = s.case;
            const v = if (stored(c, i, j)) at(s.a, s.lda, i, j) else conjIf(at(s.a, s.lda, j, i), c.herm());
            if (comptime isComplex(T)) if (c.herm() and i == j) return .init(v.re, 0);
            return v;
        }

        /// `op(A)(i, l)` of the triangular `A`, null where it is zero.
        pub fn tri(s: Self, i: usize, l: usize) ?W {
            const c = s.case;
            const r = if (c.ta.transposed()) l else i;
            const q = if (c.ta.transposed()) i else l;
            if (!stored(c, r, q)) return null;
            if (r == q and c.unit) return if (comptime isComplex(T)) .init(1, 0) else 1;
            return conjIf(at(s.a, s.lda, r, q), c.ta.conjugated());
        }

        /// One factor pair's `(i, l)` and `(l, j)`, null for a structural zero.
        pub fn factors(s: Self, second: bool, i: usize, j: usize, l: usize) [2]?W {
            const c = s.case;
            const y = if (second) s.a else s.b;
            const x = if (second) s.b else s.a;
            const ldy = if (second) s.lda else s.ldb;
            const ldx = if (second) s.ldb else s.lda;
            return switch (c.r.base()) {
                .gemm, .gemmt => .{ general(s.a, s.lda, c.ta, i, l), general(s.b, s.ldb, c.tb, l, j) },
                .symm, .hemm => if (c.left) .{ s.sym(i, l), at(s.b, s.ldb, l, j) } else .{ at(s.b, s.ldb, i, l), s.sym(l, j) },
                // `op(X)` is `n x k`; the second factor is `op(Y)'`
                .syrk, .herk, .syr2k, .her2k => blk: {
                    const two = c.r == .syr2k or c.r == .her2k;
                    const yy = if (two) y else x;
                    const ldyy = if (two) ldy else ldx;
                    // `trans` `c` conjugates `op(X)` too
                    break :blk .{ general(x, ldx, c.ta, i, l), conjIf(general(yy, ldyy, c.ta, j, l), c.herm()) };
                },
                .trmm, .trsm => if (c.left) .{ s.tri(i, l), at(s.c0, s.ldc, l, j) } else .{ at(s.c0, s.ldc, i, l), s.tri(l, j) },
                .gemm3m, .symm3m, .hemm3m => unreachable,
            };
        }

        const Elem = struct { v: W, bound: f128, terms: usize };

        /// `C(i, j)` as the routine defines it, the sum of the magnitudes in it,
        /// and how many products it sums.
        pub fn element(s: Self, i: usize, j: usize) Elem {
            const c = s.case;
            const depth = switch (c.r.base()) {
                .gemm, .gemmt, .syrk, .herk, .syr2k, .her2k => c.k,
                else => c.ka(),
            };
            const alpha = scalar(c.alpha);
            const alpha0 = c.alpha[0] == 0 and c.alpha[1] == 0;
            var sum = zeroW(W);
            var bound: f128 = 0;
            if (!alpha0) {
                const passes: usize = if (c.r == .syr2k or c.r == .her2k) 2 else 1;
                for (0..passes) |pass| {
                    var acc = zeroW(W);
                    var mags: f128 = 0;
                    for (0..depth) |l| {
                        const f = s.factors(pass == 1, i, j, l);
                        const x = f[0] orelse continue;
                        const y = f[1] orelse continue;
                        acc = reference.addW(acc, reference.mulW(x, y));
                        mags += reference.mag(x) * reference.mag(y);
                    }
                    // the second pass of `her2k` takes `conj(alpha)`
                    const al = if (pass == 1 and c.r == .her2k) reference.conjW(alpha) else alpha;
                    sum = reference.addW(sum, mulScalar(al, acc));
                    bound += magScalar(al) * mags;
                }
            }
            if (!c.triangular()) {
                const beta = scalar(c.beta);
                if (!(c.beta[0] == 0 and c.beta[1] == 0)) {
                    const cw = at(s.c0, s.ldc, i, j);
                    // `herk`'s real `beta` scales each part
                    sum = reference.addW(sum, if (c.herm() and c.rank()) partwise(c.beta[0], cw) else mulScalar(beta, cw));
                    bound += magScalar(beta) * reference.mag(cw);
                }
            }
            return .{ .v = sum, .bound = bound, .terms = depth * if (c.r == .syr2k or c.r == .her2k) @as(usize, 2) else 1 };
        }

        pub fn partwise(f: f64, x: W) W {
            const g: RefFloat(T) = @floatCast(f);
            return if (comptime isComplex(T)) .init(g * x.re, g * x.im) else g * x;
        }

        pub fn scalar(v: [2]f64) W {
            return if (comptime isComplex(T)) .init(v[0], v[1]) else v[0];
        }
        pub fn mulScalar(s: W, x: W) W {
            // `alpha` and `beta` multiply as the type's scalars
            return reference.mulW(s, x);
        }
        pub fn magScalar(s: W) f128 {
            return reference.mag(s);
        }

        /// `trsm`: the solution by substitution in the wide type, element by
        /// element against `got`, within a tolerance from the magnitudes of
        /// each step.
        pub fn checkSolve(s: Self, got: []const T, rows: usize) !void {
            const c = s.case;
            const al = testing.allocator;
            const n = c.n;
            const ka = c.ka();
            const x = try al.alloc(W, rows * n);
            defer al.free(x);
            const scale = try al.alloc(f128, rows * n);
            defer al.free(scale);
            const alpha = scalar(c.alpha);
            // `op(A)` lower: forward from the first row (left) or back from
            // the last column (right)
            const lower = c.upper == c.ta.transposed();
            for (0..n) |jj| for (0..rows) |ii| {
                // the element solved at step (ii, jj)
                const i = if (c.left and !lower) rows - 1 - ii else ii;
                const j = if (!c.left and lower) n - 1 - jj else jj;
                var r = mulScalar(alpha, at(s.c0, s.ldc, i, j));
                var mags = magScalar(alpha) * reference.mag(at(s.c0, s.ldc, i, j));
                for (0..ka) |l| {
                    // the solved neighbors: earlier in the order of substitution
                    const before = if (c.left) (if (lower) l < i else l > i) else (if (lower) l > j else l < j);
                    if (!before) continue;
                    const t = (if (c.left) s.tri(i, l) else s.tri(l, j)) orelse continue;
                    const xv = if (c.left) x[l + j * rows] else x[i + l * rows];
                    r = reference.subW(r, reference.mulW(t, xv));
                    mags += reference.mag(t) * reference.mag(xv);
                }
                const dg = (if (c.left) s.tri(i, i) else s.tri(j, j)).?;
                x[i + j * rows] = reference.divW(r, dg);
                scale[i + j * rows] = mags / @max(reference.mag(dg), reference.tiny_mag);
            };
            for (0..n) |j| for (0..rows) |i| {
                compare(T, x[i + j * rows], got[i + j * s.ldc], 4 * ka, scale[i + j * rows], true) catch return fail(c, i, j, "wrong");
            };
        }
    };
}
