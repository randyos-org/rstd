//! The level-3 test harness: each routine's dense operators in a wider type
//! as the reference, a forward-error tolerance from epsilon and the depth,
//! and the cases every set of kernels must pass: every side, triangle,
//! transpose and diagonal, shapes around every block and panel edge (on a
//! vtable with its blocking shrunk, so the edges come at small sizes, and
//! once at its own), padded leading dimensions, the `alpha` and `beta`
//! special values (`beta` 0 meets a NaN `C`, `alpha` 0 NaN operands), and the
//! special values of `check.Pattern` in each operand.
//!
//! What a routine must not read holds NaN: the padding rows, the triangle a
//! symmetric or triangular operand does not store, a unit diagonal, the
//! imaginary parts of a Hermitian diagonal. What it must not write holds a
//! sentinel: rows of `C` past `m`, what follows its last column, and the
//! triangle a rank-k update leaves.

const std = @import("std");
const types = @import("../../types.zig");
const Vtable = @import("../../Vtable.zig");
const blas = @import("../../__root__.zig");
const driver = @import("../../blas3/driver.zig");
const portable = @import("../../portable.zig");
const check1 = @import("../blas1/check.zig");
const guard = @import("../blas2/guard.zig");

const XDouble = types.XDouble;
const Float = types.Float;
const isComplex = types.isComplex;
const testing = std.testing;

/// The `Io` the tests pass (`check1.io`).
pub const io = check1.io;

/// `check1.f80_unslotted`.
pub const f80_unslotted = check1.f80_unslotted;

/// The reference's float: twice the precision.
fn WideFloat(comptime T: type) type {
    return if (@bitSizeOf(Float(T)) <= 32) f64 else f128;
}

/// The reference's type for `T`.
pub fn Wide(comptime T: type) type {
    return if (isComplex(T)) std.math.Complex(WideFloat(T)) else WideFloat(T);
}

pub const Routine = enum { gemm, symm, hemm, syrk, herk, syr2k, her2k, trmm, trsm, gemmt };

/// An operand's `op`: as it is, transposed, conjugated (`r`) or both (`c`).
pub const Trans = enum {
    n,
    t,
    r,
    c,

    pub fn transposed(t: Trans) bool {
        return t == .t or t == .c;
    }
    pub fn conjugated(t: Trans) bool {
        return t == .r or t == .c;
    }
};

pub const Where = enum { a, b, c };

pub const Case = struct {
    r: Routine = .gemm,
    /// `op(A)`, `op(B)` of `gemm` and `gemmt`; `ta` is the rank-k updates' `trans` and
    /// `trmm`'s and `trsm`'s `transa`.
    ta: Trans = .n,
    tb: Trans = .n,
    left: bool = true,
    upper: bool = true,
    unit: bool = false,
    m: usize,
    n: usize,
    k: usize = 0,
    alpha: [2]f64 = .{ 1, 0 },
    beta: [2]f64 = .{ 0, 0 },
    /// Special values in operand `where` (`check.shape`, by column-major
    /// index over its rows and columns); `c` is `trmm`'s and `trsm`'s `B`.
    p: ?check1.Pattern = null,
    where: Where = .a,

    fn rank(case: Case) bool {
        return switch (case.r) {
            .syrk, .herk, .syr2k, .her2k => true,
            else => false,
        };
    }
    /// Whether the case computes one triangle of a square `C`: the rank-k
    /// updates and `gemmt`.
    fn triangle(case: Case) bool {
        return case.rank() or case.r == .gemmt;
    }
    fn herm(case: Case) bool {
        return case.r == .hemm or case.r == .herk or case.r == .her2k;
    }
    fn triangular(case: Case) bool {
        return case.r == .trmm or case.r == .trsm;
    }
    /// The order of the symmetric or triangular `A`.
    fn ka(case: Case) usize {
        return if (case.left) case.m else case.n;
    }
};

/// Rows and columns of each operand as stored.
const Dims = struct { ra: usize, ca: usize, rb: usize, cb: usize };

fn dims(case: Case) Dims {
    const t = case.ta.transposed();
    return switch (case.r) {
        .gemm, .gemmt => .{
            .ra = if (t) case.k else case.m,
            .ca = if (t) case.m else case.k,
            .rb = if (case.tb.transposed()) case.n else case.k,
            .cb = if (case.tb.transposed()) case.k else case.n,
        },
        .symm, .hemm => .{ .ra = case.ka(), .ca = case.ka(), .rb = case.m, .cb = case.n },
        .syrk, .herk => .{ .ra = if (t) case.k else case.n, .ca = if (t) case.n else case.k, .rb = 0, .cb = 0 },
        .syr2k, .her2k => .{ .ra = if (t) case.k else case.n, .ca = if (t) case.n else case.k, .rb = if (t) case.k else case.n, .cb = if (t) case.n else case.k },
        .trmm, .trsm => .{ .ra = case.ka(), .ca = case.ka(), .rb = 0, .cb = 0 },
    };
}

/// Rows of padding below every matrix.
const pad = 3;
/// What padding holds, so a stray write shows.
const sentinel = 1234.5;

fn val(comptime T: type, re: f64, im: f64) T {
    return if (comptime isComplex(T)) .init(@floatCast(re), @floatCast(im)) else @floatCast(re);
}

fn nanOf(comptime T: type) T {
    const nan = std.math.nan(Float(T));
    return if (comptime isComplex(T)) .init(nan, nan) else nan;
}

fn wide(comptime T: type, v: T) Wide(T) {
    return if (comptime isComplex(T)) .init(v.re, v.im) else v;
}

fn conjIf(x: anytype, c: bool) @TypeOf(x) {
    return if (c) check1.conjW(x) else x;
}

fn zeroW(comptime W: type) W {
    return if (comptime isComplex(W)) .init(0, 0) else 0;
}

/// Whether `(i, j)` of `C` is one the case computes.
fn kept(case: Case, i: usize, j: usize) bool {
    if (!case.triangle()) return true;
    return if (case.upper) i <= j else i >= j;
}

/// Whether element `(i, j)` of the symmetric or triangular `A` is stored.
fn stored(case: Case, i: usize, j: usize) bool {
    return if (case.upper) i <= j else i >= j;
}

/// Runs `case` through `runner.run` and compares with the reference.
pub fn one(comptime T: type, case: Case, runner: anytype) !void {
    const al = testing.allocator;
    const d = dims(case);
    const lda = d.ra + pad;
    const ldb = d.rb + pad;
    const rc = if (case.triangle()) case.n else case.m;
    const ldc = rc + pad;
    const a = try al.alloc(T, lda * @max(d.ca, 1));
    defer al.free(a);
    const b = try al.alloc(T, ldb * @max(d.cb, 1));
    defer al.free(b);
    const c = try al.alloc(T, ldc * case.n + check1.tail);
    defer al.free(c);
    const c0 = try al.alloc(T, c.len);
    defer al.free(c0);
    check1.fill(T, a, case.m * 31 + case.k + @backingInt(case.r));
    check1.fill(T, b, case.n * 17 + case.k);
    check1.fill(T, c, case.m + case.n);
    const nan = nanOf(T);
    if (case.r == .trsm) wellConditioned(T, case, a, lda);
    if (case.p) |p| switch (case.where) {
        .a => shape(T, a, d.ra, d.ca, lda, p),
        .b => shape(T, b, d.rb, d.cb, ldb, p),
        .c => shape(T, c, rc, case.n, ldc, p),
    };
    // padding is never read, nor what a structured `A` does not store
    padNan(T, a, d.ra, lda);
    padNan(T, b, d.rb, ldb);
    if (case.r == .symm or case.r == .hemm or case.triangular()) {
        for (0..d.ca) |j| for (0..d.ra) |i| {
            if (!stored(case, i, j) or (i == j and case.triangular() and case.unit)) a[i + j * lda] = nan;
        };
        if (comptime isComplex(T)) {
            if (case.r == .hemm) for (0..d.ra) |i| {
                a[i + i * lda].im = nan.im;
            };
        }
    }
    const alpha0 = case.alpha[0] == 0 and case.alpha[1] == 0;
    const beta0 = case.beta[0] == 0 and case.beta[1] == 0;
    if (alpha0) {
        @memset(a, nan);
        @memset(b, nan);
    }
    @memset(c[ldc * case.n ..], val(T, sentinel, sentinel));
    for (0..case.n) |j| for (0..ldc) |i| {
        if (i >= rc or !kept(case, i, j)) {
            c[i + j * ldc] = val(T, sentinel, sentinel);
        } else if ((beta0 and !case.triangular()) or (alpha0 and case.triangular())) {
            c[i + j * ldc] = nan;
        } else if (comptime isComplex(T)) {
            // a Hermitian `C`'s diagonal is real: its imaginary part unread
            if ((case.r == .herk or case.r == .her2k) and i == j) c[i + j * ldc].im = nan.im;
        }
    };
    @memcpy(c0, c);

    const alpha = val(T, case.alpha[0], case.alpha[1]);
    const beta = val(T, case.beta[0], case.beta[1]);
    try runner.run(T, case, a.ptr, lda, b.ptr, ldb, alpha, beta, c.ptr, ldc);

    for (c[ldc * case.n ..]) |v| if (!check1.same(T, v, val(T, sentinel, sentinel), true)) return fail(case, ldc, case.n, "written past the end");
    const ops = Ops(T){ .case = case, .a = a, .lda = lda, .b = b, .ldb = ldb, .c0 = c0, .ldc = ldc };
    if (case.r == .trsm and !alpha0) return ops.checkSolve(c, rc);
    for (0..case.n) |j| for (0..ldc) |i| {
        const got = c[i + j * ldc];
        if (i >= rc or !kept(case, i, j)) {
            if (!check1.same(T, got, c0[i + j * ldc], true)) return fail(case, i, j, "written outside");
            continue;
        }
        const ref = ops.element(i, j);
        var want = ref.v;
        if (comptime isComplex(T)) {
            if ((case.r == .herk or case.r == .her2k) and i == j) {
                // made real unless nothing happens (`beta` 1 and no product)
                const nothing = case.beta[0] == 1 and (alpha0 or case.k == 0);
                want.im = if (nothing) c0[i + j * ldc].im else 0;
            }
        }
        // a rank-k update of special values meets them in sums of its own order
        const lax = isComplex(T) or (case.triangle() and case.p != null);
        compare(T, want, got, ref.terms, ref.bound, lax) catch return fail(case, i, j, "wrong");
    };
}

/// `got` against `want`; with `lax`, a part that should not be finite may
/// come out NaN or infinite (complex products, and solves, meet infinities
/// in an order of their own).
fn compare(comptime T: type, want: Wide(T), got: T, terms: usize, bound: f128, lax: bool) !void {
    const ulps: f128 = @floatFromInt(2 * (terms + 2));
    if (comptime isComplex(T)) {
        try part(Float(T), want.re, got.re, ulps, bound, lax);
        try part(Float(T), want.im, got.im, ulps, bound, lax);
    } else try part(T, want, got, ulps, bound, lax);
}

fn part(comptime F: type, want: anytype, got: F, ulps: f128, bound: f128, lax: bool) !void {
    const w: f128 = want;
    // past `F`'s range is infinite in `F`, of a sign the order of the sum picks
    if (lax and (!std.math.isFinite(w) or @abs(w) > std.math.floatMax(F))) {
        if (std.math.isFinite(got)) return error.TestUnexpectedResult;
        return;
    }
    try check1.component(F, w, got, ulps, bound);
}

/// The reference: the routine's operators as dense matrices, by element.
fn Ops(comptime T: type) type {
    const W = Wide(T);
    return struct {
        const Self = @This();
        case: Case,
        a: []const T,
        lda: usize,
        b: []const T,
        ldb: usize,
        c0: []const T,
        ldc: usize,

        fn at(x: []const T, ld: usize, i: usize, j: usize) W {
            return wide(T, x[i + j * ld]);
        }

        /// `op(X)(i, l)` of a general `X`.
        fn general(x: []const T, ld: usize, t: Trans, i: usize, l: usize) W {
            const v = if (t.transposed()) at(x, ld, l, i) else at(x, ld, i, l);
            return conjIf(v, t.conjugated());
        }

        /// `A`'s full symmetric or Hermitian matrix at `(i, j)`.
        fn sym(s: Self, i: usize, j: usize) W {
            const c = s.case;
            const v = if (stored(c, i, j)) at(s.a, s.lda, i, j) else conjIf(at(s.a, s.lda, j, i), c.herm());
            if (comptime isComplex(T)) if (c.herm() and i == j) return .init(v.re, 0);
            return v;
        }

        /// `op(A)(i, l)` of the triangular `A`, null where it is zero.
        fn tri(s: Self, i: usize, l: usize) ?W {
            const c = s.case;
            const r = if (c.ta.transposed()) l else i;
            const q = if (c.ta.transposed()) i else l;
            if (!stored(c, r, q)) return null;
            if (r == q and c.unit) return if (comptime isComplex(T)) .init(1, 0) else 1;
            return conjIf(at(s.a, s.lda, r, q), c.ta.conjugated());
        }

        /// One factor pair's `(i, l)` and `(l, j)`, null for a structural zero.
        fn factors(s: Self, second: bool, i: usize, j: usize, l: usize) [2]?W {
            const c = s.case;
            const y = if (second) s.a else s.b;
            const x = if (second) s.b else s.a;
            const ldy = if (second) s.lda else s.ldb;
            const ldx = if (second) s.ldb else s.lda;
            return switch (c.r) {
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
            };
        }

        const Ref = struct { v: W, bound: f128, terms: usize };

        /// `C(i, j)` as the routine defines it, the sum of the magnitudes in it,
        /// and how many products it sums.
        fn element(s: Self, i: usize, j: usize) Ref {
            const c = s.case;
            const depth = switch (c.r) {
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
                        acc = check1.addW(acc, check1.mulW(x, y));
                        mags += check1.mag(x) * check1.mag(y);
                    }
                    // the second pass of `her2k` takes `conj(alpha)`
                    const al = if (pass == 1 and c.r == .her2k) check1.conjW(alpha) else alpha;
                    sum = check1.addW(sum, mulScalar(al, acc));
                    bound += magScalar(al) * mags;
                }
            }
            if (!c.triangular()) {
                const beta = scalar(c.beta);
                if (!(c.beta[0] == 0 and c.beta[1] == 0)) {
                    const cw = at(s.c0, s.ldc, i, j);
                    // `herk`'s real `beta` scales each part
                    sum = check1.addW(sum, if (c.herm() and c.rank()) partwise(c.beta[0], cw) else mulScalar(beta, cw));
                    bound += magScalar(beta) * check1.mag(cw);
                }
            }
            return .{ .v = sum, .bound = bound, .terms = depth * if (c.r == .syr2k or c.r == .her2k) @as(usize, 2) else 1 };
        }

        fn partwise(f: f64, x: W) W {
            const g: WideFloat(T) = @floatCast(f);
            return if (comptime isComplex(T)) .init(g * x.re, g * x.im) else g * x;
        }

        fn scalar(v: [2]f64) W {
            return if (comptime isComplex(T)) .init(v[0], v[1]) else v[0];
        }
        fn mulScalar(s: W, x: W) W {
            // `alpha` and `beta` multiply as the type's scalars
            return check1.mulW(s, x);
        }
        fn magScalar(s: W) f128 {
            return check1.mag(s);
        }

        /// `trsm`: the solution by substitution in the wide type, element by
        /// element against `got`, within a tolerance from the magnitudes of
        /// each step.
        fn checkSolve(s: Self, got: []const T, rows: usize) !void {
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
                var mags = magScalar(alpha) * check1.mag(at(s.c0, s.ldc, i, j));
                for (0..ka) |l| {
                    // the solved neighbors: earlier in the order of substitution
                    const before = if (c.left) (if (lower) l < i else l > i) else (if (lower) l > j else l < j);
                    if (!before) continue;
                    const t = (if (c.left) s.tri(i, l) else s.tri(l, j)) orelse continue;
                    const xv = if (c.left) x[l + j * rows] else x[i + l * rows];
                    r = check1.subW(r, check1.mulW(t, xv));
                    mags += check1.mag(t) * check1.mag(xv);
                }
                const dg = (if (c.left) s.tri(i, i) else s.tri(j, j)).?;
                x[i + j * rows] = check1.divW(r, dg);
                scale[i + j * rows] = mags / @max(check1.mag(dg), 1e-300);
            };
            for (0..n) |j| for (0..rows) |i| {
                compare(T, x[i + j * rows], got[i + j * s.ldc], 4 * ka, scale[i + j * rows], true) catch return fail(c, i, j, "wrong");
            };
        }
    };
}

/// `trsm`'s `A`: off the diagonal a quarter of the order's share, so the
/// solve is well conditioned whatever the triangle; the diagonal at least 1.
fn wellConditioned(comptime T: type, case: Case, a: []T, lda: usize) void {
    const ka = case.ka();
    const f: Float(T) = @floatCast(1.0 / @as(f64, @floatFromInt(4 * ka)));
    for (0..ka) |j| for (0..ka) |i| {
        const p = &a[i + j * lda];
        if (i != j) {
            p.* = if (comptime isComplex(T)) .init(p.re * f, p.im * f) else p.* * f;
        } else if (comptime isComplex(T)) {
            p.* = .init(p.re / 4 + std.math.copysign(@as(Float(T), 1), p.re), p.im / 4);
        } else p.* = p.* / 4 + std.math.copysign(@as(T, 1), p.*);
    };
}

/// Operand `p` shapes, the `rows` by `cols` of `x` at `ld`.
fn shape(comptime T: type, x: []T, rows: usize, cols: usize, ld: usize, p: check1.Pattern) void {
    for (0..cols) |j| for (0..rows) |i| {
        x[i + j * ld] = check1.shape(T, x[i + j * ld], i + j * rows, p);
    };
}

/// NaN in the rows past `rows` of each column.
fn padNan(comptime T: type, x: []T, rows: usize, ld: usize) void {
    var j: usize = 0;
    while (j * ld < x.len) : (j += 1) @memset(x[j * ld + rows .. @min(x.len, (j + 1) * ld)], nanOf(T));
}

/// Whether `fail` describes the case (`accelerateCases` stops after a few).
var say_fails = true;

fn fail(case: Case, i: usize, j: usize, what: []const u8) anyerror {
    if (!say_fails) return error.TestUnexpectedResult;
    if (case.p) |p| {
        p.describe();
        check1.say("in {s}, ", .{@tagName(case.where)});
    }
    check1.say("{s} at ({d},{d}): {s} ta={s} tb={s} left={} upper={} unit={} m={d} n={d} k={d} alpha={d},{d} beta={d},{d}\n", .{ what, i, j, @tagName(case.r), @tagName(case.ta), @tagName(case.tb), case.left, case.upper, case.unit, case.m, case.n, case.k, case.alpha[0], case.alpha[1], case.beta[0], case.beta[1] });
    return error.TestUnexpectedResult;
}

// ---- runners --------------------------------------------------------------------

/// Calls the driver of `r` on a fixed vtable through a poisoned allocator, so
/// a packing buffer read before it is written shows as NaN. The routine is
/// fixed at compile time, so a test builds only its own driver.
pub fn Driver(comptime r: Routine) type {
    return struct {
        vt: *const Vtable,

        pub fn run(self: @This(), comptime T: type, case: Case, a: [*]T, lda: usize, b: [*]T, ldb: usize, alpha: T, beta: T, c: [*]T, ldc: usize) !void {
            var poison: guard.Poison = .{ .child = testing.allocator };
            const al = poison.allocator();
            const vt = self.vt;
            const herm = comptime isComplex(T);
            switch (r) {
                .gemm => try driver.gemm(al, io, T, vt, case.ta.transposed(), case.tb.transposed(), .of(case.ta.conjugated(), case.tb.conjugated()), case.m, case.n, case.k, alpha, a, lda, b, ldb, beta, c, ldc),
                .gemmt => try driver.gemmt(al, io, T, vt, case.upper, case.ta.transposed(), case.tb.transposed(), .of(case.ta.conjugated(), case.tb.conjugated()), case.n, case.k, alpha, a, lda, b, ldb, beta, c, ldc),
                .symm => try driver.symm(al, io, T, vt, false, case.left, case.upper, case.m, case.n, alpha, a, lda, b, ldb, beta, c, ldc),
                .hemm => if (herm) try driver.symm(al, io, T, vt, true, case.left, case.upper, case.m, case.n, alpha, a, lda, b, ldb, beta, c, ldc) else unreachable,
                .syrk => try driver.syrk(al, io, T, vt, false, case.upper, case.ta.transposed(), case.n, case.k, alpha, a, lda, null, 0, beta, c, ldc),
                .syr2k => try driver.syrk(al, io, T, vt, false, case.upper, case.ta.transposed(), case.n, case.k, alpha, a, lda, b, ldb, beta, c, ldc),
                .herk => if (herm) try driver.syrk(al, io, T, vt, true, case.upper, case.ta.transposed(), case.n, case.k, alpha, a, lda, null, 0, beta, c, ldc) else unreachable,
                .her2k => if (herm) try driver.syrk(al, io, T, vt, true, case.upper, case.ta.transposed(), case.n, case.k, alpha, a, lda, b, ldb, beta, c, ldc) else unreachable,
                .trmm, .trsm => if (@hasDecl(driver, "trmm")) {
                    const args = .{ al, io, T, vt, case.left, case.upper, case.ta.transposed(), case.ta.conjugated(), case.unit, case.m, case.n, alpha, a, lda, c, ldc };
                    if (r == .trmm) try @call(.auto, driver.trmm, args) else try @call(.auto, driver.trsm, args);
                } else return error.SkipZigTest,
            }
        }
    };
}

/// Calls the public routines (`blas.gemm` and the rest) through a poisoned
/// allocator; real `r` and `c` stand in for `n` and `t` half the time. As
/// `Driver`, for one routine `r`.
pub fn Iface(comptime r: Routine) type {
    return struct {
        pub fn run(_: @This(), comptime T: type, case: Case, a: [*]T, lda: usize, b: [*]T, ldb: usize, alpha: T, beta: T, c: [*]T, ldc: usize) !void {
            var poison: guard.Poison = .{ .child = testing.allocator };
            const al = poison.allocator();
            const odd = (case.m + case.n + case.k) % 2 == 1;
            const ta = trans(T, case.ta, odd, case.rank());
            const tb = trans(T, case.tb, !odd, false);
            const m: isize = @intCast(case.m);
            const n: isize = @intCast(case.n);
            const k: isize = @intCast(case.k);
            const ia: isize = @intCast(lda);
            const ib: isize = @intCast(ldb);
            const ic: isize = @intCast(ldc);
            const side: blas.Side = if (case.left) .l else .r;
            const uplo: blas.Uplo = if (case.upper) .u else .l;
            const diag: blas.Diag = if (case.unit) .u else .n;
            const herm = comptime isComplex(T);
            switch (r) {
                .gemm => try blas.gemm(al, io, T, ta, tb, m, n, k, alpha, a, ia, b, ib, beta, c, ic),
                .gemmt => try blas.gemmt(al, io, T, uplo, ta, tb, n, k, alpha, a, ia, b, ib, beta, c, ic),
                .symm => try blas.symm(al, io, T, side, uplo, m, n, alpha, a, ia, b, ib, beta, c, ic),
                .hemm => if (herm) try blas.hemm(al, io, T, side, uplo, m, n, alpha, a, ia, b, ib, beta, c, ic) else unreachable,
                .syrk => try blas.syrk(al, io, T, uplo, ta, n, k, alpha, a, ia, beta, c, ic),
                .syr2k => try blas.syr2k(al, io, T, uplo, ta, n, k, alpha, a, ia, b, ib, beta, c, ic),
                .herk => if (herm) try blas.herk(al, io, T, uplo, ta, n, k, alpha.re, a, ia, beta.re, c, ic) else unreachable,
                .her2k => if (herm) try blas.her2k(al, io, T, uplo, ta, n, k, alpha, a, ia, b, ib, beta.re, c, ic) else unreachable,
                .trmm => try blas.trmm(al, io, T, side, uplo, ta, diag, m, n, alpha, a, ia, c, ic),
                .trsm => try blas.trsm(al, io, T, side, uplo, ta, diag, m, n, alpha, a, ia, c, ic),
            }
        }

        /// `t` as `blas.Trans`; real `r` and `c` for `n` and `t` where `alt`.
        fn trans(comptime T: type, t: Trans, alt: bool, rank: bool) blas.Trans {
            if (comptime isComplex(T)) return switch (t) {
                inline else => |v| @field(blas.Trans, @tagName(v)),
            };
            // the rank-k updates take no `r`
            return if (t.transposed()) (if (alt) .c else .t) else (if (alt and !rank) .r else .n);
        }
    };
}

// ---- cases -----------------------------------------------------------------------

/// `vt` with its `T` blocking shrunk so every block edge comes at small
/// sizes; static, as a vtable is large for a rig's stack.
var shrunk: Vtable = undefined;

/// The options of `r` over `T`, each a case to fill in: every side,
/// triangle, transpose and diagonal it takes.
fn variants(comptime T: type, r: Routine, buf: []Case) []Case {
    const complex = comptime isComplex(T);
    const all_t: []const Trans = if (complex) &.{ .n, .t, .r, .c } else &.{ .n, .t };
    var len: usize = 0;
    switch (r) {
        // the pairs of the shape sweep; the rest come with the scalars
        .gemm => for (all_t, 0..) |ta, i| {
            buf[len] = .{ .r = r, .ta = ta, .tb = all_t[(i + 1) % all_t.len], .m = 0, .n = 0 };
            buf[len + 1] = .{ .r = r, .ta = ta, .tb = ta, .m = 0, .n = 0 };
            len += 2;
        },
        .gemmt => for ([_]bool{ true, false }) |upper| for (all_t, 0..) |ta, i| {
            buf[len] = .{ .r = r, .ta = ta, .tb = all_t[(i + 1) % all_t.len], .upper = upper, .m = 0, .n = 0 };
            buf[len + 1] = .{ .r = r, .ta = ta, .tb = ta, .upper = upper, .m = 0, .n = 0 };
            len += 2;
        },
        .symm, .hemm => for ([_]bool{ true, false }) |left| for ([_]bool{ true, false }) |upper| {
            buf[len] = .{ .r = r, .left = left, .upper = upper, .m = 0, .n = 0 };
            len += 1;
        },
        .syrk, .herk, .syr2k, .her2k => for ([_]bool{ true, false }) |upper| for ([_]bool{ false, true }) |t| {
            const tt: Trans = if (!t) .n else if (r == .herk or r == .her2k) .c else .t;
            buf[len] = .{ .r = r, .upper = upper, .ta = tt, .m = 0, .n = 0 };
            len += 1;
        },
        .trmm, .trsm => for ([_]bool{ true, false }) |left| for ([_]bool{ true, false }) |upper| for (all_t) |t| for ([_]bool{ false, true }) |unit| {
            buf[len] = .{ .r = r, .left = left, .upper = upper, .ta = t, .unit = unit, .m = 0, .n = 0 };
            len += 1;
        },
    }
    return buf[0..len];
}

fn sized(v: Case, m: usize, n: usize, k: usize, alpha: [2]f64, beta: [2]f64) Case {
    var c = v;
    c.m = if (v.triangle()) n else m;
    c.n = n;
    c.k = k;
    c.alpha = alpha;
    c.beta = if (v.r == .herk or v.r == .her2k) .{ beta[0], 0 } else beta;
    if (v.r == .herk) c.alpha = .{ alpha[0], 0 };
    return c;
}

/// The scalars: `alpha` and `beta` special and plain, complex ones too.
fn scalarsOf(comptime T: type) []const [2][2]f64 {
    const real = [_][2][2]f64{
        .{ .{ -0.75, 0 }, .{ 0.5, 0 } }, .{ .{ 0, 0 }, .{ 0.5, 0 } },   .{ .{ 1, 0 }, .{ 0, 0 } }, .{ .{ 0, 0 }, .{ 0, 0 } },
        .{ .{ 1, 0 }, .{ 1, 0 } },       .{ .{ -0.75, 0 }, .{ 0, 0 } }, .{ .{ 0, 0 }, .{ 1, 0 } },
    };
    const complex = real ++ [_][2][2]f64{ .{ .{ 0.7, -0.9 }, .{ 1.3, -1.1 } }, .{ .{ 0, 1 }, .{ 0, -1 } } };
    return if (comptime isComplex(T)) &complex else &real;
}

/// Every case of `r` over `T` on `vt`'s kernels.
pub fn routineCases(comptime T: type, vt: *const Vtable, comptime r: Routine) !void {
    const pre = comptime portable.letter(T).? ++ "gemm_";
    const own = driver.Blocking.of(T, vt);
    const um = own.um;
    const un = own.un;

    shrunk = vt.*;
    @field(shrunk, pre ++ "p") = @intCast(2 * um);
    @field(shrunk, pre ++ "q") = 5;
    @field(shrunk, pre ++ "r") = @intCast(6 * un + 1);
    const p = 2 * um;
    const q = 5;
    const rr = 6 * un + 1;
    const small: Driver(r) = .{ .vt = &shrunk };
    const full: Driver(r) = .{ .vt = vt };
    const rig = check1.rig();
    var vbuf: [32]Case = undefined;
    const vs = variants(T, r, &vbuf);

    // every block and panel edge, on the shrunk blocking
    const ms = [_]usize{ 1, @max(1, um - 1), um + 1, p + 1, 2 * p + um + 1 };
    const ns = [_]usize{ 1, un + 1, 3 * un + 1, rr + 1, 2 * rr + 3 };
    const ks = [_]usize{ 1, q + 1, 2 * q + 3 };
    for (vs, 0..) |v, vi| for (ms, 0..) |m, mi| for (ns, 0..) |n, ni| for (ks, 0..) |k, ki| {
        // on the rigs, and past the first variants, a diagonal of the sizes
        const thin = rig or vi >= 4;
        if (thin and (mi + ni + ki + vi) % 4 != 0) continue;
        // the depth of the triangular and symmetric routines is `m` or `n`
        if (ki > 0 and !(r == .gemm or v.triangle())) continue;
        // a triangle of `C` has one order, from either list
        if (v.triangle() and mi != 0 and ni != 0) continue;
        try one(T, sized(v, m, if (v.triangle() and ni == 0) m else n, k, .{ -0.75, 0.25 }, .{ 0.5, -0.5 }), small);
    };

    // the special scalars, on a shape with every remainder; every transpose
    // pair of `gemm`
    const sm = 2 * um + 3;
    const sn = 2 * un + 1;
    for (scalarsOf(T)) |s| {
        if (r == .gemm) {
            const all_t: []const Trans = if (comptime isComplex(T)) &.{ .n, .t, .r, .c } else &.{ .n, .t };
            for (all_t) |ta| for (all_t) |tb| {
                try one(T, sized(.{ .r = r, .ta = ta, .tb = tb, .m = 0, .n = 0 }, sm, sn, 7, s[0], s[1]), full);
            };
        } else for (vs) |v| try one(T, sized(v, sm, sn, 7, s[0], s[1]), full);
    }

    // special values in each operand, and NaN scalars
    for (vs, 0..) |v, vi| for ([_]Where{ .a, .b, .c }) |where| {
        if (rig and vi % 4 != 0) continue;
        const has_b = r == .gemm or r == .gemmt or r == .symm or r == .hemm or r == .syr2k or r == .her2k;
        if (where == .b and !has_b) continue;
        // a solve's special values go in its right side
        if (r == .trsm and where == .a) continue;
        const c = sized(v, sm, sn, 9, .{ -0.75, 0.25 }, .{ 0.5, 0 });
        const d = dims(c);
        const len = switch (where) {
            .a => d.ra * d.ca,
            .b => d.rb * d.cb,
            .c => c.m * c.n,
        };
        var buf: [13]check1.Pattern = undefined;
        for (check1.patternsFor(len, &buf)) |pat| {
            var cp = c;
            cp.p = pat;
            cp.where = where;
            try one(T, cp, full);
        }
    };
    const nan = std.math.nan(f64);
    for ([_][2][2]f64{ .{ .{ nan, 0 }, .{ 0.5, 0 } }, .{ .{ 1.5, 0 }, .{ nan, 0 } } }) |s| {
        if (vs[0].triangular() and std.math.isNan(s[1][0])) continue;
        try one(T, sized(vs[vs.len - 1], sm, sn, 7, s[0], s[1]), full);
    }

    // past each of the vtable's own blocks once (row blocks split in two)
    try one(T, sized(vs[0], own.p + um + 1, 3 * un + 2, own.q + 1, .{ 1.5, 0 }, .{ -1, 0 }), full);
    // past `R` (a left side, so `A` is small; a rank-k update's `C` would be
    // `R` square, and the shrunk blocking crosses it; `trmm` and `trsm` pack
    // `Q x R`, the reference's whole buffer, too much beside a rig's tests);
    // natively, whole blocks of `P` and `Q` and then halves
    if (!vs[0].triangle() and !(rig and vs[0].triangular())) try one(T, sized(if (r == .gemm) vs[vs.len - 1] else vs[0], um + 1, own.r + 2 * un + 1, 3, .{ 1.25, 0 }, .{ 0.5, 0 }), full);
    if (!rig) try one(T, sized(vs[vs.len / 2], 2 * own.p + um + 3, un + 1, 2 * own.q + 3, .{ -0.5, 0 }, .{ 1, 0 }), full);
}

/// The routines over `T`.
fn routinesOf(comptime T: type) []const Routine {
    return if (comptime isComplex(T)) std.enums.values(Routine) else &.{ .gemm, .symm, .syrk, .syr2k, .trmm, .trsm, .gemmt };
}

/// Whether the drivers have `r` yet.
fn ported(r: Routine) bool {
    return switch (r) {
        .trmm, .trsm => @hasDecl(driver, "trmm"),
        else => true,
    };
}

/// `routineCases` for every routine of every precision; `XDouble` only where
/// it is a type of its own.
pub fn coreCases(vt: *const Vtable) !void {
    inline for (.{ f32, f64, types.Complex32, types.Complex64 } ++ check1.x_real ++ check1.x_complex) |T| {
        inline for (comptime routinesOf(T)) |r| if (comptime ported(r)) try routineCases(T, vt, r);
    }
}

/// `routineCases` of one routine for every precision it has.
pub fn routineAll(vt: *const Vtable, comptime r: Routine) !void {
    inline for (.{ f32, f64, types.Complex32, types.Complex64 } ++ check1.x_real ++ check1.x_complex) |T| {
        if (comptime std.mem.findScalar(Routine, routinesOf(T), r) != null) try routineCases(T, vt, r);
    }
}

/// Cases of `r` through an interface (`Iface`), at sizes small and past the
/// generic blocking.
pub fn interfaceCases(comptime T: type, comptime r: Routine) !void {
    var vbuf: [32]Case = undefined;
    for (variants(T, r, &vbuf)) |v| for (scalarsOf(T)) |s| for ([_][3]usize{ .{ 1, 1, 1 }, .{ 9, 7, 5 }, .{ 37, 29, 70 } }) |d| {
        try one(T, sized(v, d[0], d[1], d[2], s[0], s[1]), Iface(r){});
    };
}

/// Through the public routine (`Iface`), as Accelerate runs it: the special
/// scalars, the special values of `check.Pattern` in each operand and NaN
/// scalars, at small sizes and past its small-problem paths. Carries on past
/// a failure, saying each, up to `max_said` per routine.
pub fn accelerateCases(comptime T: type, comptime r: Routine) !void {
    var vbuf: [32]Case = undefined;
    const vs = variants(T, r, &vbuf);
    const all_t: []const Trans = if (comptime isComplex(T)) &.{ .n, .t, .r, .c } else &.{ .n, .t };
    var bad: Bad = .{};
    defer say_fails = true;
    for ([_][3]usize{ .{ 19, 9, 9 }, .{ 70, 45, 40 } }) |d| {
        for (scalarsOf(T)) |s| {
            if (r == .gemm) {
                for (all_t) |ta| for (all_t) |tb| bad.note(one(T, sized(.{ .r = r, .ta = ta, .tb = tb, .m = 0, .n = 0 }, d[0], d[1], d[2], s[0], s[1]), Iface(r){}));
            } else for (vs) |v| bad.note(one(T, sized(v, d[0], d[1], d[2], s[0], s[1]), Iface(r){}));
        }
        for (vs) |v| for ([_]Where{ .a, .b, .c }) |where| {
            const has_b = r == .gemm or r == .gemmt or r == .symm or r == .hemm or r == .syr2k or r == .her2k;
            if (where == .b and !has_b) continue;
            if (r == .trsm and where == .a) continue;
            const c = sized(v, d[0], d[1], d[2], .{ -0.75, 0.25 }, .{ 0.5, 0 });
            const dd = dims(c);
            const len = switch (where) {
                .a => dd.ra * dd.ca,
                .b => dd.rb * dd.cb,
                .c => c.m * c.n,
            };
            var buf: [13]check1.Pattern = undefined;
            for (check1.patterns(len, &buf)) |pat| {
                var cp = c;
                cp.p = pat;
                cp.where = where;
                bad.note(one(T, cp, Iface(r){}));
            }
        };
        const nan = std.math.nan(f64);
        for (vs) |v| for ([_][2][2]f64{ .{ .{ nan, 0 }, .{ 0.5, 0 } }, .{ .{ 1.5, 0 }, .{ nan, 0 } } }) |s| {
            if (v.triangular() and std.math.isNan(s[1][0])) continue;
            bad.note(one(T, sized(v, d[0], d[1], d[2], s[0], s[1]), Iface(r){}));
        };
    }
    if (bad.err) |e| return e;
}

/// How many failing cases `accelerateCases` describes per routine.
const max_said = 8;

const Bad = struct {
    err: ?anyerror = null,
    count: usize = 0,

    fn note(b: *Bad, r: anyerror!void) void {
        r catch |e| {
            b.err = e;
            b.count += 1;
            if (b.count == max_said) {
                check1.say("(more not shown)\n", .{});
                say_fails = false;
            }
        };
    }
};
