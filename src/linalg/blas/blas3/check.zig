//! The level-3 test harness: a plain gemm in a wider type as the reference,
//! a forward-error tolerance from epsilon and `k`, and the cases every set of
//! kernels must pass: all transposes, shapes around every block and panel
//! edge (on a vtable with its blocking shrunk, so the edges come at small
//! sizes, and once at its own), padded leading dimensions, and the `alpha`
//! and `beta` special values (`beta` 0 meets a NaN `C`, `alpha` 0 NaN `A`
//! and `B`), and the special values of `check.Pattern` in each operand. The
//! padding rows of `A` and `B` hold NaN; rows of `C` past `m`, and what
//! follows its last column, must stay untouched.

const std = @import("std");
const types = @import("../types.zig");
const Vtable = @import("../Vtable.zig");
const level3 = @import("level3.zig");
const check1 = @import("../blas1/check.zig");
const guard = @import("../blas2/guard.zig");

const XDouble = types.XDouble;
const testing = std.testing;

/// The reference's type: twice the precision.
pub fn Wide(comptime T: type) type {
    return if (@bitSizeOf(T) <= 32) f64 else f128;
}

pub const Case = struct {
    ta: bool,
    tb: bool,
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    beta: f64,
    /// Special values in operand `where` (`check.shape`, by column-major
    /// index over its rows and columns).
    p: ?check1.Pattern = null,
    where: enum { a, b, c } = .a,
};

/// Rows of padding below every matrix.
const pad = 3;
/// What padding holds, so a stray write shows.
const sentinel = 1234.5;

/// Runs `case` through `runner.run` and compares with the reference.
pub fn one(comptime T: type, case: Case, runner: anytype) !void {
    const W = Wide(T);
    const al = testing.allocator;
    const ra = if (case.ta) case.k else case.m;
    const ca = if (case.ta) case.m else case.k;
    const rb = if (case.tb) case.n else case.k;
    const cb = if (case.tb) case.k else case.n;
    const lda = ra + pad;
    const ldb = rb + pad;
    const ldc = case.m + pad;
    const a = try al.alloc(T, lda * @max(ca, 1));
    defer al.free(a);
    const b = try al.alloc(T, ldb * @max(cb, 1));
    defer al.free(b);
    const c = try al.alloc(T, ldc * case.n + check1.tail);
    defer al.free(c);
    const c0 = try al.alloc(T, c.len);
    defer al.free(c0);
    check1.fill(T, a, case.m * 31 + case.k);
    check1.fill(T, b, case.n * 17 + case.k);
    check1.fill(T, c, case.m + case.n);
    const nan = std.math.nan(T);
    if (case.p) |p| switch (case.where) {
        .a => shape(T, a, ra, ca, lda, p),
        .b => shape(T, b, rb, cb, ldb, p),
        .c => shape(T, c, case.m, case.n, ldc, p),
    };
    // padding is never read
    padNan(T, a, ra, lda);
    padNan(T, b, rb, ldb);
    if (case.alpha == 0) {
        @memset(a, nan);
        @memset(b, nan);
    }
    @memset(c[ldc * case.n ..], sentinel);
    for (0..case.n) |j| {
        for (0..ldc) |i| {
            if (i >= case.m) c[i + j * ldc] = sentinel else if (case.beta == 0) c[i + j * ldc] = nan;
        }
    }
    @memcpy(c0, c);

    const alpha: T = @floatCast(case.alpha);
    const beta: T = @floatCast(case.beta);
    try runner.run(T, case.ta, case.tb, case.m, case.n, case.k, alpha, a.ptr, lda, b.ptr, ldb, beta, c.ptr, ldc);

    for (c[ldc * case.n ..]) |v| if (v != sentinel) return fail(case, ldc, case.n, "written past the end");
    for (0..case.n) |j| for (0..ldc) |i| {
        const got = c[i + j * ldc];
        if (i >= case.m) {
            if (got != sentinel) return fail(case, i, j, "padding written");
            continue;
        }
        var s: W = 0;
        var bound: f128 = 0;
        if (case.alpha != 0) for (0..case.k) |l| {
            const x: W = a[if (case.ta) l + i * lda else i + l * lda];
            const y: W = b[if (case.tb) j + l * ldb else l + j * ldb];
            s += x * y;
            bound += check1.mag(@as(f128, x)) * check1.mag(@as(f128, y));
        };
        const cw: W = c0[i + j * ldc];
        const aw: W = alpha;
        const bw: W = beta;
        var want = aw * s;
        bound = check1.mag(@as(f128, aw)) * bound;
        if (case.beta != 0) {
            want += bw * cw;
            bound += check1.mag(@as(f128, bw)) * check1.mag(@as(f128, cw));
        }
        check1.component(T, want, got, @floatFromInt(2 * (case.k + 2)), bound) catch return fail(case, i, j, "wrong");
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
    while (j * ld < x.len) : (j += 1) @memset(x[j * ld + rows .. @min(x.len, (j + 1) * ld)], std.math.nan(T));
}

fn fail(case: Case, i: usize, j: usize, what: []const u8) anyerror {
    if (case.p) |p| {
        p.describe();
        check1.say("in {s}, ", .{@tagName(case.where)});
    }
    check1.say("{s} at ({d},{d}): ta={} tb={} m={d} n={d} k={d} alpha={d} beta={d}\n", .{ what, i, j, case.ta, case.tb, case.m, case.n, case.k, case.alpha, case.beta });
    return error.TestUnexpectedResult;
}

const transposes = [_][2]bool{ .{ false, false }, .{ false, true }, .{ true, false }, .{ true, true } };
const scalars = [_][2]f64{ .{ -0.75, 0.5 }, .{ 0, 0.5 }, .{ 1, 0 }, .{ 0, 0 }, .{ 1, 1 }, .{ -0.75, 0 }, .{ 0, 1 } };

/// Calls the driver on a fixed vtable through a poisoned allocator, so a
/// packing buffer read before it is written shows as NaN.
const Driver = struct {
    vt: *const Vtable,

    pub fn run(self: Driver, comptime T: type, ta: bool, tb: bool, m: usize, n: usize, k: usize, alpha: T, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize) !void {
        var poison: guard.Poison = .{ .child = testing.allocator };
        try level3.gemm(poison.allocator(), T, self.vt, ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
};

/// `vt` with its `T` blocking shrunk so every block edge comes at small
/// sizes; static, as a vtable is large for a rig's stack.
var shrunk: Vtable = undefined;

/// Every case for `T` on `vt`'s kernels.
pub fn gemmCases(comptime T: type, vt: *const Vtable) !void {
    const pre = comptime @import("../portable.zig").letter(T).? ++ "gemm_";
    const own = level3.Blocking.of(T, vt);
    const um = own.um;
    const un = own.un;

    shrunk = vt.*;
    @field(shrunk, pre ++ "p") = @intCast(2 * um);
    @field(shrunk, pre ++ "q") = 5;
    @field(shrunk, pre ++ "r") = @intCast(6 * un + 1);
    const p = 2 * um;
    const q = 5;
    const r = 6 * un + 1;
    const small: Driver = .{ .vt = &shrunk };
    const ms = [_]usize{ 1, @max(1, um - 1), um + 1, p + 1, 2 * p + um + 1 };
    const ns = [_]usize{ 1, un + 1, 3 * un + 1, r + 1, 2 * r + 3 };
    const ks = [_]usize{ 1, q + 1, 2 * q + 3 };
    for (transposes) |t| for (ms) |m| for (ns) |n| for (ks) |k| {
        try one(T, .{ .ta = t[0], .tb = t[1], .m = m, .n = n, .k = k, .alpha = -0.75, .beta = 0.5 }, small);
    };

    // the special scalars, on a shape with every remainder
    const full: Driver = .{ .vt = vt };
    for (transposes) |t| for (scalars) |s| {
        try one(T, .{ .ta = t[0], .tb = t[1], .m = 2 * um + 3, .n = 2 * un + 1, .k = 7, .alpha = s[0], .beta = s[1] }, full);
    };

    // special values in each operand, and NaN scalars
    for (transposes) |t| for ([_]@FieldType(Case, "where"){ .a, .b, .c }) |where| {
        const m = 2 * um + 3;
        const n = 2 * un + 1;
        const k = 9;
        const len = switch (where) {
            .a => m * k,
            .b => k * n,
            .c => m * n,
        };
        var buf: [13]check1.Pattern = undefined;
        for (check1.patternsFor(len, &buf)) |pat| {
            try one(T, .{ .ta = t[0], .tb = t[1], .m = m, .n = n, .k = k, .alpha = -0.75, .beta = 0.5, .p = pat, .where = where }, full);
        }
    };
    for ([_][2]f64{ .{ std.math.nan(f64), 0.5 }, .{ 1.5, std.math.nan(f64) } }) |s| {
        try one(T, .{ .ta = false, .tb = true, .m = 2 * um + 3, .n = 2 * un + 1, .k = 7, .alpha = s[0], .beta = s[1] }, full);
    }

    // past each of the vtable's own blocks once (row blocks split in two)
    for ([_][2]bool{ .{ false, false }, .{ true, true } }) |t| {
        try one(T, .{ .ta = t[0], .tb = t[1], .m = own.p + um + 1, .n = 3 * un + 2, .k = own.q + 1, .alpha = 1.5, .beta = -1 }, full);
    }
    // past `R`; natively, whole blocks of `P` and `Q` and then halves
    try one(T, .{ .ta = false, .tb = true, .m = um + 1, .n = own.r + 2 * un + 1, .k = 3, .alpha = 1.25, .beta = 0.5 }, full);
    if (!check1.rig()) try one(T, .{ .ta = true, .tb = false, .m = 2 * own.p + um + 3, .n = un + 1, .k = 2 * own.q + 3, .alpha = -0.5, .beta = 1 }, full);
}

/// `gemmCases` for every real precision; `XDouble` only where it is a type
/// of its own.
pub fn coreCases(vt: *const Vtable) !void {
    try gemmCases(f32, vt);
    try gemmCases(f64, vt);
    if (XDouble != f64) try gemmCases(XDouble, vt);
}

/// Cases through an interface `gemm` (`interface.zig`'s signature), with
/// real `r` and `c` standing in for `n` and `t` half the time.
pub fn interfaceCases(comptime T: type, comptime Trans: type, gemm: anytype) !void {
    const Iface = struct {
        pub fn run(_: @This(), comptime U: type, ta: bool, tb: bool, m: usize, n: usize, k: usize, alpha: U, a: [*]const U, lda: usize, b: [*]const U, ldb: usize, beta: U, c: [*]U, ldc: usize) !void {
            const odd = (m + n + k) % 2 == 1;
            const tra: Trans = if (ta) (if (odd) .c else .t) else (if (odd) .r else .n);
            const trb: Trans = if (tb) (if (odd) .t else .c) else (if (odd) .n else .r);
            var poison: guard.Poison = .{ .child = testing.allocator };
            try gemm(poison.allocator(), U, tra, trb, @intCast(m), @intCast(n), @intCast(k), alpha, a, @intCast(lda), b, @intCast(ldb), beta, c, @intCast(ldc));
        }
    };
    for (transposes) |t| for (scalars) |s| for ([_][3]usize{ .{ 1, 1, 1 }, .{ 9, 7, 5 }, .{ 37, 29, 70 } }) |d| {
        try one(T, .{ .ta = t[0], .tb = t[1], .m = d[0], .n = d[1], .k = d[2], .alpha = s[0], .beta = s[1] }, Iface{});
    };
}
