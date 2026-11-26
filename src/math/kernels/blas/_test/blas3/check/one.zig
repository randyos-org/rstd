//! One case against the reference: the run, the comparison and its operands' shapes.

const std = @import("std");
const types = @import("../../../types.zig");
const check = @import("../../check.zig");
const bitwise = @import("../../bitwise.zig");
const pattern = @import("../../pattern.zig");
const reference = @import("../../reference.zig");
const c3 = @import("../check.zig");
const c3ops = @import("ops.zig");

const Float = types.Float;
const isComplex = types.isComplex;
const testing = std.testing;
const Case = c3.Case;
const Ref = c3.Ref;
const dims = c3.dims;
const kept = c3.kept;
const nanOf = c3.nanOf;
const pad = c3.pad;
const sentinel = c3.sentinel;
const stored = c3.stored;
const val = c3.val;
const Ops = c3ops.Ops;

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
    const c = try al.alloc(T, ldc * case.n + check.tail);
    defer al.free(c);
    const c0 = try al.alloc(T, c.len);
    defer al.free(c0);
    check.fill(T, a, case.m * 31 + case.k + @backingInt(case.r));
    check.fill(T, b, case.n * 17 + case.k);
    check.fill(T, c, case.m + case.n);
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
    if (case.r.base() == .symm or case.r.base() == .hemm or case.triangular()) {
        for (0..d.ca) |j| for (0..d.ra) |i| {
            if (!stored(case, i, j) or (i == j and case.triangular() and case.unit)) a[i + j * lda] = nan;
        };
        if (comptime isComplex(T)) {
            if (case.r.base() == .hemm) for (0..d.ra) |i| {
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

    for (c[ldc * case.n ..]) |v| if (!bitwise.same(T, v, val(T, sentinel, sentinel), true)) return fail(case, ldc, case.n, "written past the end");
    const ops = Ops(T){ .case = case, .a = a, .lda = lda, .b = b, .ldb = ldb, .c0 = c0, .ldc = ldc };
    if (case.r == .trsm and !alpha0) return ops.checkSolve(c, rc);
    for (0..case.n) |j| for (0..ldc) |i| {
        const got = c[i + j * ldc];
        if (i >= rc or !kept(case, i, j)) {
            if (!bitwise.same(T, got, c0[i + j * ldc], true)) return fail(case, i, j, "written outside");
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
        if (case.three()) {
            compare3m(T, case, want, got, ref.terms, ref.bound) catch return fail(case, i, j, "wrong");
        } else compare(T, want, got, ref.terms, ref.bound, lax) catch return fail(case, i, j, "wrong");
    };
}

/// `got` against `want`; with `lax`, a part that should not be finite may
/// come out NaN or infinite (complex products, and solves, meet infinities
/// in an order of their own).
pub fn compare(comptime T: type, want: Ref(T), got: T, terms: usize, bound: f128, lax: bool) !void {
    const ulps = reference.sumUlps(terms);
    if (comptime isComplex(T)) {
        try part(Float(T), want.re, got.re, ulps, bound, lax);
        try part(Float(T), want.im, got.im, ulps, bound, lax);
    } else try part(T, want, got, ulps, bound, lax);
}

/// `compare` for a 3M routine: the part it forms from three products, a
/// difference of them (the imaginary one, or the real one where exactly one
/// operand is conjugated), takes the error of three sums of `terms`, and of
/// the copies' `alpha * op(B)` and part sums (Higham 1992).
pub fn compare3m(comptime T: type, case: Case, want: Ref(T), got: T, terms: usize, bound: f128) !void {
    if (comptime !isComplex(T)) unreachable;
    const ulps = reference.sumUlps(terms);
    const ulps3 = reference.sumUlps(3 * terms + 2);
    const real3 = case.r == .gemm3m and case.ta.conjugated() != case.tb.conjugated();
    try part(Float(T), want.re, got.re, if (real3) ulps3 else ulps, bound, true);
    try part(Float(T), want.im, got.im, if (real3) ulps else ulps3, bound, true);
}

pub fn part(comptime F: type, want: anytype, got: F, ulps: f128, bound: f128, lax: bool) !void {
    const w: f128 = want;
    // past `F`'s range is infinite in `F`, of a sign the order of the sum picks
    if (lax and (!std.math.isFinite(w) or @abs(w) > std.math.floatMax(F))) {
        if (std.math.isFinite(got)) return error.TestUnexpectedResult;
        return;
    }
    try pattern.component(F, w, got, ulps, bound);
}

/// `trsm`'s `A`: off the diagonal a quarter of the order's share, so the
/// solve is well conditioned whatever the triangle; the diagonal at least 1.
pub fn wellConditioned(comptime T: type, case: Case, a: []T, lda: usize) void {
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
pub fn shape(comptime T: type, x: []T, rows: usize, cols: usize, ld: usize, p: pattern.Pattern) void {
    for (0..cols) |j| for (0..rows) |i| {
        x[i + j * ld] = pattern.shape(T, x[i + j * ld], i + j * rows, p);
    };
}

/// NaN in the rows past `rows` of each column.
pub fn padNan(comptime T: type, x: []T, rows: usize, ld: usize) void {
    var j: usize = 0;
    while (j * ld < x.len) : (j += 1) @memset(x[j * ld + rows .. @min(x.len, (j + 1) * ld)], nanOf(T));
}

pub fn fail(case: Case, i: usize, j: usize, what: []const u8) anyerror {
    if (case.p) |p| {
        p.describe();
        check.say("in {s}, ", .{@tagName(case.where)});
    }
    check.say("{s} at ({d},{d}): {s} ta={s} tb={s} left={} upper={} unit={} m={d} n={d} k={d} alpha={d},{d} beta={d},{d}: ", .{ what, i, j, @tagName(case.r), @tagName(case.ta), @tagName(case.tb), case.left, case.upper, case.unit, case.m, case.n, case.k, case.alpha[0], case.alpha[1], case.beta[0], case.beta[1] });
    return error.TestUnexpectedResult;
}
