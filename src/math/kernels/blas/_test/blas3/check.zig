//! The level-3 test harness: each routine's dense operators in a wider type
//! as the reference, and the cases every set of kernels must pass, around every
//! block edge, with NaN where a routine must not read and sentinels where it must not write.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../../types.zig");
const Vtable = @import("../../Vtable.zig");
const driver = @import("../../blas3/driver.zig");
const driver3m = @import("../../blas3/driver3m.zig");
const portable = @import("../../portable.zig");
const check = @import("../check.zig");
const pattern = @import("../pattern.zig");
const reference = @import("../reference.zig");
const c3one = @import("check/one.zig");
const c3run = @import("check/runners.zig");

const Float = types.Float;
const isComplex = types.isComplex;
const one = c3one.one;
const Driver = c3run.Driver;
const Iface = c3run.Iface;
pub const callAny = c3run.callAny;

/// The `Io` the tests pass (`check.io`).
pub const io = check.io;

/// `check.f80_unslotted`.
pub const f80_unslotted = check.f80_unslotted;

/// The reference's float: twice the precision.
pub fn RefFloat(comptime T: type) type {
    return if (@bitSizeOf(Float(T)) <= 32) f64 else f128;
}

/// The reference's type for `T`.
pub fn Ref(comptime T: type) type {
    return if (isComplex(T)) std.math.Complex(RefFloat(T)) else RefFloat(T);
}

pub const Routine = enum {
    gemm,
    symm,
    hemm,
    syrk,
    herk,
    syr2k,
    her2k,
    trmm,
    trsm,
    gemmt,
    gemm3m,
    symm3m,
    hemm3m,

    /// The routine a 3M one computes, by three real products; others, itself.
    pub fn base(r: Routine) Routine {
        return switch (r) {
            .gemm3m => .gemm,
            .symm3m => .symm,
            .hemm3m => .hemm,
            else => r,
        };
    }
};

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
    pub fn of(trans: bool, conj: bool) Trans {
        return if (trans) (if (conj) .c else .t) else (if (conj) .r else .n);
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
    /// Special values in operand `where` (`pattern.shape`, by column-major
    /// index over its rows and columns); `c` is `trmm`'s and `trsm`'s `B`.
    p: ?pattern.Pattern = null,
    where: Where = .a,

    pub fn rank(case: Case) bool {
        return switch (case.r) {
            .syrk, .herk, .syr2k, .her2k => true,
            else => false,
        };
    }
    /// Whether the case computes one triangle of a square `C`: the rank-k
    /// updates and `gemmt`.
    pub fn triangle(case: Case) bool {
        return case.rank() or case.r == .gemmt;
    }
    pub fn herm(case: Case) bool {
        return case.r.base() == .hemm or case.r == .herk or case.r == .her2k;
    }
    pub fn three(case: Case) bool {
        return case.r != case.r.base();
    }
    pub fn triangular(case: Case) bool {
        return case.r == .trmm or case.r == .trsm;
    }
    /// The order of the symmetric or triangular `A`.
    pub fn ka(case: Case) usize {
        return if (case.left) case.m else case.n;
    }
};

/// Rows and columns of each operand as stored.
pub const Dims = struct { ra: usize, ca: usize, rb: usize, cb: usize };

pub fn dims(case: Case) Dims {
    const t = case.ta.transposed();
    return switch (case.r.base()) {
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
        .gemm3m, .symm3m, .hemm3m => unreachable,
    };
}

/// Rows of padding below every matrix.
pub const pad = 3;
/// What padding holds, so a stray write shows.
pub const sentinel = 1234.5;

pub const val = check.val;

pub fn nanOf(comptime T: type) T {
    const nan = std.math.nan(Float(T));
    return if (comptime isComplex(T)) .init(nan, nan) else nan;
}

pub fn toRef(comptime T: type, v: T) Ref(T) {
    return if (comptime isComplex(T)) .init(v.re, v.im) else v;
}

pub fn conjIf(x: anytype, c: bool) @TypeOf(x) {
    return if (c) reference.conjW(x) else x;
}

pub fn zeroW(comptime W: type) W {
    return if (comptime isComplex(W)) .init(0, 0) else 0;
}

/// Whether `(i, j)` of `C` is one the case computes.
pub fn kept(case: Case, i: usize, j: usize) bool {
    if (!case.triangle()) return true;
    return if (case.upper) i <= j else i >= j;
}

/// Whether element `(i, j)` of the symmetric or triangular `A` is stored.
pub fn stored(case: Case, i: usize, j: usize) bool {
    return if (case.upper) i <= j else i >= j;
}

// ---- cases -----------------------------------------------------------------------

/// `vt` with its `T` blocking shrunk so every block edge comes at small
/// sizes; static, as a vtable is large for a rig's stack.
pub var shrunk: Vtable = undefined;

/// The options of `r` over `T`, each a case to fill in: every side,
/// triangle, transpose and diagonal it takes.
pub fn variants(comptime T: type, r: Routine, buf: []Case) []Case {
    const complex = comptime isComplex(T);
    const all_t: []const Trans = if (complex) &.{ .n, .t, .r, .c } else &.{ .n, .t };
    var len: usize = 0;
    switch (r) {
        // the pairs of the shape sweep; the rest come with the scalars
        .gemm, .gemm3m => for (all_t, 0..) |ta, i| {
            buf[len] = .{ .r = r, .ta = ta, .tb = all_t[(i + 1) % all_t.len], .m = 0, .n = 0 };
            buf[len + 1] = .{ .r = r, .ta = ta, .tb = ta, .m = 0, .n = 0 };
            len += 2;
        },
        .gemmt => for ([_]bool{ true, false }) |upper| for (all_t, 0..) |ta, i| {
            buf[len] = .{ .r = r, .ta = ta, .tb = all_t[(i + 1) % all_t.len], .upper = upper, .m = 0, .n = 0 };
            buf[len + 1] = .{ .r = r, .ta = ta, .tb = ta, .upper = upper, .m = 0, .n = 0 };
            len += 2;
        },
        .symm, .hemm, .symm3m, .hemm3m => for ([_]bool{ true, false }) |left| for ([_]bool{ true, false }) |upper| {
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

pub fn sized(v: Case, m: usize, n: usize, k: usize, alpha: [2]f64, beta: [2]f64) Case {
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
pub fn scalarsOf(comptime T: type) []const [2][2]f64 {
    const real = [_][2][2]f64{
        .{ .{ -0.75, 0 }, .{ 0.5, 0 } }, .{ .{ 0, 0 }, .{ 0.5, 0 } },   .{ .{ 1, 0 }, .{ 0, 0 } }, .{ .{ 0, 0 }, .{ 0, 0 } },
        .{ .{ 1, 0 }, .{ 1, 0 } },       .{ .{ -0.75, 0 }, .{ 0, 0 } }, .{ .{ 0, 0 }, .{ 1, 0 } },
    };
    const complex = real ++ [_][2][2]f64{ .{ .{ 0.7, -0.9 }, .{ 1.3, -1.1 } }, .{ .{ 0, 1 }, .{ 0, -1 } } };
    return if (comptime isComplex(T)) &complex else &real;
}

/// Shrinks `vt`'s `T` blocking (the 3M one for `three`) around the unrolls
/// `um` and `un`, so every block and panel edge comes at small sizes.
pub fn shrinkBlocking(comptime T: type, vt: *Vtable, comptime three: bool, um: usize, un: usize) void {
    const g = comptime portable.Gemm.of(T, three);
    portable.paramPtr(vt, g, .p).* = @intCast(2 * um);
    portable.paramPtr(vt, g, .q).* = 5;
    portable.paramPtr(vt, g, .r).* = @intCast(6 * un + 1);
}

/// Every case of `r` over `T` on `vt`'s kernels.
pub fn routineCases(comptime T: type, vt: *const Vtable, comptime r: Routine) !void {
    // a 3M routine's blocking is its own
    const three = comptime r != r.base();
    const own = if (three) driver3m.blocking(T, vt) else driver.Blocking.of(T, vt);
    const um = own.um;
    const un = own.un;

    shrunk = vt.*;
    shrinkBlocking(T, &shrunk, three, um, un);
    const p = 2 * um;
    const q = 5;
    const rr = 6 * un + 1;
    const small: Driver(r) = .{ .vt = &shrunk };
    const full: Driver(r) = .{ .vt = vt };
    const rig = check.rig();
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
        if (ki > 0 and !(r.base() == .gemm or v.triangle())) continue;
        // a triangle of `C` has one order, from either list
        if (v.triangle() and mi != 0 and ni != 0) continue;
        try one(T, sized(v, m, if (v.triangle() and ni == 0) m else n, k, .{ -0.75, 0.25 }, .{ 0.5, -0.5 }), small);
    };

    // the special scalars, special values and NaN scalars, on a shape with
    // every remainder
    try edgeCases(T, r, vs, .{ .m = 2 * um + 3, .n = 2 * un + 1, .k = 7, .k_special = 9, .thin = rig }, full);

    // past each of the vtable's own blocks once (row blocks split in two)
    try one(T, sized(vs[0], own.p + um + 1, 3 * un + 2, own.q + 1, .{ 1.5, 0 }, .{ -1, 0 }), full);
    // past `R` (a left side, so `A` is small; a rank-k update's `C` would be
    // `R` square, and the shrunk blocking crosses it; `trmm` and `trsm` pack
    // `Q x R`, the reference's whole buffer, too much beside a rig's tests);
    // one row on the freestanding and UEFI rigs, where `B`, `C` and its copy
    // take about a MiB a row each, too much for their heap at `unroll_m + 1`
    // rows (the shrunk sweep crosses `R` with whole panels); natively, whole
    // blocks of `P` and `Q` and then halves
    const past_r_rows = switch (sysinfo.target.os.tag) {
        .freestanding, .uefi => 1,
        else => um + 1,
    };
    if (!vs[0].triangle() and !(rig and vs[0].triangular())) try one(T, sized(if (r.base() == .gemm) vs[vs.len - 1] else vs[0], past_r_rows, own.r + 2 * un + 1, 3, .{ 1.25, 0 }, .{ 0.5, 0 }), full);
    if (!rig) try one(T, sized(vs[vs.len / 2], 2 * own.p + um + 3, un + 1, 2 * own.q + 3, .{ -0.5, 0 }, .{ 1, 0 }), full);
}

/// The routines over `T`.
pub fn routinesOf(comptime T: type) []const Routine {
    return if (comptime isComplex(T)) std.enums.values(Routine) else &.{ .gemm, .symm, .syrk, .syr2k, .trmm, .trsm, .gemmt };
}

/// Whether the drivers have `r` yet.
pub fn ported(r: Routine) bool {
    return switch (r) {
        .trmm, .trsm => @hasDecl(driver, "trmm"),
        else => true,
    };
}

/// `routineCases` for every routine of every precision; `XDouble` only where
/// it is a type of its own.
pub fn coreCases(vt: *const Vtable) !void {
    inline for (check.all_types) |T| {
        inline for (comptime routinesOf(T)) |r| if (comptime ported(r)) try routineCases(T, vt, r);
    }
}

/// `routineCases` of one routine for every precision it has.
pub fn routineAll(vt: *const Vtable, comptime r: Routine) !void {
    inline for (check.all_types) |T| {
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
/// scalars, the special values of `pattern.Pattern` in each operand and NaN
/// scalars, at small sizes and past its small-problem paths. Carries on past
/// a failure, saying each, up to `max_said` per routine.
pub fn accelerateCases(comptime T: type, comptime r: Routine) !void {
    var vbuf: [32]Case = undefined;
    const vs = variants(T, r, &vbuf);
    var failed: ?anyerror = null;
    check.capKeep(max_said);
    defer check.uncap();
    for ([_][3]usize{ .{ 19, 9, 9 }, .{ 70, 45, 40 } }) |d| {
        try edgeCases(T, r, vs, .{ .m = d[0], .n = d[1], .k = d[2], .k_special = d[2], .every_pattern = true, .nan_every_variant = true, .keep = &failed }, Iface(r){});
    }
    if (failed) |e| return e;
}

/// How many failing cases `accelerateCases` describes per routine.
pub const max_said = 8;

/// The shapes and reach of `edgeCases`.
pub const Edges = struct {
    m: usize,
    n: usize,
    /// The depth of the scalar cases, and of the special values'.
    k: usize,
    k_special: usize,
    /// Every fourth variant's special values only.
    thin: bool = false,
    /// `pattern.patterns`, not the rigs' few (`patternsFor`).
    every_pattern: bool = false,
    /// The NaN scalars on every variant, not the last alone.
    nan_every_variant: bool = false,
    /// Carries on past a failure, kept here (`check.keep`), not returned.
    keep: ?*?anyerror = null,
};

/// The special scalars (every transpose pair for `gemm`), the special
/// values of `pattern.Pattern` in each operand, and NaN scalars, of `r`'s
/// variants `vs` through `runner`.
pub fn edgeCases(comptime T: type, comptime r: Routine, vs: []const Case, e: Edges, runner: anytype) !void {
    for (scalarsOf(T)) |s| {
        if (r.base() == .gemm) {
            const all_t: []const Trans = if (comptime isComplex(T)) &.{ .n, .t, .r, .c } else &.{ .n, .t };
            for (all_t) |ta| for (all_t) |tb| {
                try edge(e, one(T, sized(.{ .r = r, .ta = ta, .tb = tb, .m = 0, .n = 0 }, e.m, e.n, e.k, s[0], s[1]), runner));
            };
        } else for (vs) |v| try edge(e, one(T, sized(v, e.m, e.n, e.k, s[0], s[1]), runner));
    }
    const has_b = r.base() == .gemm or r == .gemmt or r.base() == .symm or r.base() == .hemm or r == .syr2k or r == .her2k;
    for (vs, 0..) |v, vi| for ([_]Where{ .a, .b, .c }) |where| {
        if (e.thin and vi % 4 != 0) continue;
        if (where == .b and !has_b) continue;
        // a solve's special values go in its right side
        if (r == .trsm and where == .a) continue;
        const c = sized(v, e.m, e.n, e.k_special, .{ -0.75, 0.25 }, .{ 0.5, 0 });
        const d = dims(c);
        const len = switch (where) {
            .a => d.ra * d.ca,
            .b => d.rb * d.cb,
            .c => c.m * c.n,
        };
        var buf: [13]pattern.Pattern = undefined;
        for (if (e.every_pattern) pattern.patterns(len, &buf) else pattern.patternsFor(len, &buf)) |pat| {
            var cp = c;
            cp.p = pat;
            cp.where = where;
            try edge(e, one(T, cp, runner));
        }
    };
    const nan = std.math.nan(f64);
    for (if (e.nan_every_variant) vs else vs[vs.len - 1 ..]) |v| for ([_][2][2]f64{ .{ .{ nan, 0 }, .{ 0.5, 0 } }, .{ .{ 1.5, 0 }, .{ nan, 0 } } }) |s| {
        if (v.triangular() and std.math.isNan(s[1][0])) continue;
        try edge(e, one(T, sized(v, e.m, e.n, e.k, s[0], s[1]), runner));
    };
}

/// One of `edgeCases`'s results: kept or returned, as `e` says.
pub fn edge(e: Edges, r: anyerror!void) anyerror!void {
    if (e.keep) |f| return check.keep(f, r);
    return r;
}
