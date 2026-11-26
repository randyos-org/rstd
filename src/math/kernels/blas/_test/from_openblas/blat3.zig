//! The level-3 programs `test/[sdcz]blat3.f`, `[cz]blat3_3m.f` and their CBLAS
//! `ctest/` versions with their `.dat` tables, OpenBLAS v0.3.34-317; `SMMCH`
//! keeps the largest ratio over all columns and gauges 3M at three times the bound.

const std = @import("std");
const ob = @import("__root__.zig");
const blas = @import("../../__root__.zig");
const l3 = @import("../../blas3/core.zig");
const check = @import("../check.zig");
const mk = @import("blat3/make.zig");
const cl = @import("blat3/calls.zig");
const sc = @import("blat3/schk.zig");

const Float = ob.Float;
const isComplex = ob.isComplex;
const testing = std.testing;
const C32 = std.math.Complex(f32);
const C64 = std.math.Complex(f64);
const mmch = mk.mmch;
const asSide = cl.asSide;
const asTrans = cl.asTrans;
const asUplo = cl.asUplo;
const chk1 = sc.chk1;
const chk2 = sc.chk2;
const chk3 = sc.chk3;
const chk4 = sc.chk4;
const chk5 = sc.chk5;

pub const nmax = 65;

// ---- parameters ----------------------------------------------------------------

/// A `.dat` file's sweep: values of `N`, `ALPHA`, `BETA` and the threshold.
pub const Params = struct {
    thresh: f64 = 16,
    n: []const usize,
    alpha: []const [2]f64,
    beta: []const [2]f64,
};

/// `test/[sd]blat3.dat`: the first 6 of its `N` (`NIDIM` is 6).
pub const real_dat: Params = .{
    .n = &.{ 0, 1, 2, 3, 7, 31 },
    .alpha = &.{ .{ 0, 0 }, .{ 1, 0 }, .{ 0.7, 0 } },
    .beta = &.{ .{ 0, 0 }, .{ 1, 0 }, .{ 1.3, 0 } },
};

/// `test/[cz]blat3.dat`, `[cz]blat3_3m.dat`.
pub const complex_dat: Params = .{
    .n = &.{ 0, 1, 2, 3, 7, 31 },
    .alpha = &.{ .{ 0, 0 }, .{ 1, 0 }, .{ 0.7, -0.9 } },
    .beta = &.{ .{ 0, 0 }, .{ 1, 0 }, .{ 1.3, -1.1 } },
};

/// `ctest/[sdcz]in3` (and `[cz]in3_3m`, the same as `[cz]in3`): `N` by type.
pub fn ctestN(comptime T: type) []const usize {
    return switch (T) {
        f32 => &.{ 0, 1, 2, 3, 5, 9, 35 },
        f64 => &.{ 1, 2, 3, 5, 7, 9, 35 },
        C32 => &.{ 0, 1, 2, 3, 5, 9 },
        else => &.{ 0, 1, 2, 3, 5, 9, 35 },
    };
}

/// The rigs' `N`: the example data file in each program's header.
pub const rig_n = [_]usize{ 0, 1, 2, 3, 5, 9 };

pub fn dat(comptime T: type) Params {
    return if (comptime isComplex(T)) complex_dat else real_dat;
}

/// The `.dat` table natively, ctest's after it; the header's on the rigs.
pub fn tables(comptime T: type) [2]?Params {
    var p = dat(T);
    if (check.rig()) {
        p.n = &rig_n;
        return .{ p, null };
    }
    var c = p;
    c.n = ctestN(T);
    return .{ p, c };
}

pub const Routine = enum { gemm, gemm3m, symm, hemm, trmm, trsm, syrk, herk, syr2k, her2k };

pub fn letter(comptime T: type) u8 {
    return switch (Float(T)) {
        f32 => if (isComplex(T)) 'C' else 'S',
        else => if (isComplex(T)) 'Z' else 'D',
    };
}

pub fn chk(comptime T: type, comptime r: Routine, p: Params) !void {
    const f = switch (r) {
        .gemm, .gemm3m => chk1,
        .symm, .hemm => chk2,
        .trmm, .trsm => chk3,
        .syrk, .herk => chk4,
        .syr2k, .her2k => chk5,
    };
    try f(T, r, p);
}

/// `r` over `types`, each on its tables, after its error exits.
pub fn sweep(comptime r: Routine, comptime types: anytype) !void {
    inline for (types) |T| {
        try errorExits(T, r);
        for (tables(T)) |t| if (t) |p| try chk(T, r, p);
    }
}

/// Real `gemm` on the `N` of `n`, after its error exits, for `check.eachCore`.
pub fn all(n: []const usize) !void {
    inline for (.{ f32, f64 }) |T| {
        try errorExits(T, .gemm);
        var p = dat(T);
        p.n = n;
        try chk(T, .gemm, p);
    }
}

// ---- SCHKE ----------------------------------------------------------------------

pub const E = ob.expectError;

/// `SCHKE`/`CCHKE` for `r`: one bad argument per call, the first a bad letter
/// (not translatable) or one `r` does not take.
pub fn errorExits(comptime T: type, comptime r: Routine) !void {
    const al = testing.failing_allocator;
    var a: [4]T = @splat(ob.zero(T));
    var b: [4]T = @splat(ob.zero(T));
    var c: [4]T = @splat(ob.zero(T));
    const z = ob.zero(T);
    const rz: Float(T) = 0;
    const nm = @tagName(r);
    switch (r) {
        .gemm, .gemm3m => {
            const G = struct { ta: blas.Trans, tb: blas.Trans, m: isize, n: isize, k: isize, lda: isize, ldb: isize, ldc: isize, e: anyerror };
            const cases = [_]G{
                .{ .ta = .n, .tb = .n, .m = -1, .n = 0, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidM },
                .{ .ta = .n, .tb = .t, .m = -1, .n = 0, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidM },
                .{ .ta = .t, .tb = .n, .m = -1, .n = 0, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidM },
                .{ .ta = .t, .tb = .t, .m = -1, .n = 0, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidM },
                .{ .ta = .n, .tb = .n, .m = 0, .n = -1, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidN },
                .{ .ta = .n, .tb = .t, .m = 0, .n = -1, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidN },
                .{ .ta = .t, .tb = .n, .m = 0, .n = -1, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidN },
                .{ .ta = .t, .tb = .t, .m = 0, .n = -1, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidN },
                .{ .ta = .n, .tb = .n, .m = 0, .n = 0, .k = -1, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidK },
                .{ .ta = .n, .tb = .t, .m = 0, .n = 0, .k = -1, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidK },
                .{ .ta = .t, .tb = .n, .m = 0, .n = 0, .k = -1, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidK },
                .{ .ta = .t, .tb = .t, .m = 0, .n = 0, .k = -1, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidK },
                .{ .ta = .n, .tb = .n, .m = 2, .n = 0, .k = 0, .lda = 1, .ldb = 1, .ldc = 2, .e = error.InvalidLda },
                .{ .ta = .n, .tb = .t, .m = 2, .n = 0, .k = 0, .lda = 1, .ldb = 1, .ldc = 2, .e = error.InvalidLda },
                .{ .ta = .t, .tb = .n, .m = 0, .n = 0, .k = 2, .lda = 1, .ldb = 2, .ldc = 1, .e = error.InvalidLda },
                .{ .ta = .t, .tb = .t, .m = 0, .n = 0, .k = 2, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidLda },
                .{ .ta = .n, .tb = .n, .m = 0, .n = 0, .k = 2, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidLdb },
                .{ .ta = .t, .tb = .n, .m = 0, .n = 0, .k = 2, .lda = 2, .ldb = 1, .ldc = 1, .e = error.InvalidLdb },
                .{ .ta = .n, .tb = .t, .m = 0, .n = 2, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidLdb },
                .{ .ta = .t, .tb = .t, .m = 0, .n = 2, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidLdb },
                .{ .ta = .n, .tb = .n, .m = 2, .n = 0, .k = 0, .lda = 2, .ldb = 1, .ldc = 1, .e = error.InvalidLdc },
                .{ .ta = .n, .tb = .t, .m = 2, .n = 0, .k = 0, .lda = 2, .ldb = 1, .ldc = 1, .e = error.InvalidLdc },
                .{ .ta = .t, .tb = .n, .m = 2, .n = 0, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidLdc },
                .{ .ta = .t, .tb = .t, .m = 2, .n = 0, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidLdc },
            };
            for (cases) |g| try E(g.e, @field(l3, nm)(al, ob.io, T, g.ta, g.tb, g.m, g.n, g.k, z, &a, g.lda, &b, g.ldb, z, &c, g.ldc), nm);
        },
        .symm, .hemm => {
            const S = struct { s: u8, u: u8, m: isize, n: isize, lda: isize, ldb: isize, ldc: isize, e: anyerror };
            const cases = [_]S{
                .{ .s = 'L', .u = 'U', .m = -1, .n = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidM },
                .{ .s = 'R', .u = 'U', .m = -1, .n = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidM },
                .{ .s = 'L', .u = 'L', .m = -1, .n = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidM },
                .{ .s = 'R', .u = 'L', .m = -1, .n = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidM },
                .{ .s = 'L', .u = 'U', .m = 0, .n = -1, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidN },
                .{ .s = 'R', .u = 'U', .m = 0, .n = -1, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidN },
                .{ .s = 'L', .u = 'L', .m = 0, .n = -1, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidN },
                .{ .s = 'R', .u = 'L', .m = 0, .n = -1, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidN },
                .{ .s = 'L', .u = 'U', .m = 2, .n = 0, .lda = 1, .ldb = 2, .ldc = 2, .e = error.InvalidLda },
                .{ .s = 'R', .u = 'U', .m = 0, .n = 2, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidLda },
                .{ .s = 'L', .u = 'L', .m = 2, .n = 0, .lda = 1, .ldb = 2, .ldc = 2, .e = error.InvalidLda },
                .{ .s = 'R', .u = 'L', .m = 0, .n = 2, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidLda },
                .{ .s = 'L', .u = 'U', .m = 2, .n = 0, .lda = 2, .ldb = 1, .ldc = 2, .e = error.InvalidLdb },
                .{ .s = 'R', .u = 'U', .m = 2, .n = 0, .lda = 1, .ldb = 1, .ldc = 2, .e = error.InvalidLdb },
                .{ .s = 'L', .u = 'L', .m = 2, .n = 0, .lda = 2, .ldb = 1, .ldc = 2, .e = error.InvalidLdb },
                .{ .s = 'R', .u = 'L', .m = 2, .n = 0, .lda = 1, .ldb = 1, .ldc = 2, .e = error.InvalidLdb },
                .{ .s = 'L', .u = 'U', .m = 2, .n = 0, .lda = 2, .ldb = 2, .ldc = 1, .e = error.InvalidLdc },
                .{ .s = 'R', .u = 'U', .m = 2, .n = 0, .lda = 1, .ldb = 2, .ldc = 1, .e = error.InvalidLdc },
                .{ .s = 'L', .u = 'L', .m = 2, .n = 0, .lda = 2, .ldb = 2, .ldc = 1, .e = error.InvalidLdc },
                .{ .s = 'R', .u = 'L', .m = 2, .n = 0, .lda = 1, .ldb = 2, .ldc = 1, .e = error.InvalidLdc },
            };
            for (cases) |g| try E(g.e, @field(l3, nm)(al, ob.io, T, asSide(g.s), asUplo(g.u), g.m, g.n, z, &a, g.lda, &b, g.ldb, z, &c, g.ldc), nm);
        },
        .trmm, .trsm => {
            const S = struct { s: u8, u: u8, t: u8, m: isize, n: isize, lda: isize, ldb: isize, e: anyerror };
            var cases: [32]S = undefined;
            var i: usize = 0;
            // INFOT 5, 6, 9, 11 over every side, triangle and trans N, T
            for ("UL") |u| for ("LR") |s| for ("NT") |t| {
                const left = s == 'L';
                cases[i] = .{ .s = s, .u = u, .t = t, .m = -1, .n = 0, .lda = 1, .ldb = 1, .e = error.InvalidM };
                cases[i + 1] = .{ .s = s, .u = u, .t = t, .m = 0, .n = -1, .lda = 1, .ldb = 1, .e = error.InvalidN };
                cases[i + 2] = if (left)
                    .{ .s = s, .u = u, .t = t, .m = 2, .n = 0, .lda = 1, .ldb = 2, .e = error.InvalidLda }
                else
                    .{ .s = s, .u = u, .t = t, .m = 0, .n = 2, .lda = 1, .ldb = 1, .e = error.InvalidLda };
                cases[i + 3] = .{ .s = s, .u = u, .t = t, .m = 2, .n = 0, .lda = if (left) 2 else 1, .ldb = 1, .e = error.InvalidLdb };
                i += 4;
            };
            for (cases) |g| try E(g.e, @field(l3, nm)(al, ob.io, T, asSide(g.s), asUplo(g.u), asTrans(g.t), .n, g.m, g.n, z, &a, g.lda, &b, g.ldb), nm);
        },
        .syrk, .herk, .syr2k, .her2k => {
            const two = r == .syr2k or r == .her2k;
            // the second letter of each pair is the transpose this routine takes
            const tt: u8 = if (r == .herk or r == .her2k) 'C' else 'T';
            // a transpose it does not take: `herk`'s T, complex `syrk`'s C
            if (comptime isComplex(T)) try E(error.InvalidTrans, call2(T, r, al, 'U', if (tt == 'C') 'T' else 'C', 0, 0, 1, 1, 1, &a, &b, &c, z, rz), nm);
            const S = struct { u: u8, t: u8, n: isize, k: isize, lda: isize, ldb: isize, ldc: isize, e: anyerror };
            var cases: [20]S = undefined;
            var i: usize = 0;
            for ("UL") |u| for ([_]u8{ 'N', tt }) |t| {
                const nt = t == 'N';
                cases[i] = .{ .u = u, .t = t, .n = -1, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidN };
                cases[i + 1] = .{ .u = u, .t = t, .n = 0, .k = -1, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidK };
                cases[i + 2] = if (nt)
                    .{ .u = u, .t = t, .n = 2, .k = 0, .lda = 1, .ldb = 1, .ldc = 2, .e = error.InvalidLda }
                else
                    .{ .u = u, .t = t, .n = 0, .k = 2, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidLda };
                cases[i + 3] = if (nt)
                    .{ .u = u, .t = t, .n = 2, .k = 0, .lda = 2, .ldb = 1, .ldc = 2, .e = error.InvalidLdb }
                else
                    .{ .u = u, .t = t, .n = 0, .k = 2, .lda = 2, .ldb = 1, .ldc = 1, .e = error.InvalidLdb };
                cases[i + 4] = if (nt)
                    .{ .u = u, .t = t, .n = 2, .k = 0, .lda = 2, .ldb = 2, .ldc = 1, .e = error.InvalidLdc }
                else
                    .{ .u = u, .t = t, .n = 2, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidLdc };
                i += 5;
            };
            for (cases) |g| {
                if (!two and g.e == error.InvalidLdb) continue;
                try E(g.e, call2(T, r, al, g.u, g.t, g.n, g.k, g.lda, g.ldb, g.ldc, &a, &b, &c, z, rz), nm);
            }
        },
    }
}

/// The rank-k updates with `alpha` and `beta` zero.
pub fn call2(comptime T: type, comptime r: Routine, al: std.mem.Allocator, u: u8, t: u8, n: isize, k: isize, lda: isize, ldb: isize, ldc: isize, a: [*]T, b: [*]T, c: [*]T, z: T, rz: Float(T)) l3.Error!void {
    return switch (r) {
        .syrk => l3.syrk(al, ob.io, T, asUplo(u), asTrans(t), n, k, z, a, lda, z, c, ldc),
        .herk => l3.herk(al, ob.io, T, asUplo(u), asTrans(t), n, k, rz, a, lda, rz, c, ldc),
        .syr2k => l3.syr2k(al, ob.io, T, asUplo(u), asTrans(t), n, k, z, a, lda, b, ldb, z, c, ldc),
        .her2k => l3.her2k(al, ob.io, T, asUplo(u), asTrans(t), n, k, z, a, lda, b, ldb, rz, c, ldc),
        else => unreachable,
    };
}

// ---- tests ----------------------------------------------------------------------

test "blat3: SMMCH self-check on exact data" {
    inline for (.{ f32, f64, C32, C64 }) |T| {
        const R = ob.Ref(T);
        const n = 32;
        var ab: [n * (n + 1)]T = undefined;
        var c: [n]T = @splat(ob.zero(T));
        var cc: [n]T = undefined;
        var ct: [n]R = undefined;
        var g: [n]Float(R) = undefined;
        for (1..n + 1) |j| for (1..n + 1) |i| {
            ab[(i - 1) + (j - 1) * n] = ob.val(T, @floatFromInt(@as(i64, @intCast(i)) - @as(i64, @intCast(j)) + 1), 0);
            if (i < j) ab[(i - 1) + (j - 1) * n] = ob.zero(T);
        };
        const bcol = ab[n * n ..];
        for ([_][2]u8{ "NN".*, "NT".*, "TN".*, "TT".* }) |t| {
            const rev = t[0] == 'T';
            for (1..n + 1) |jj| {
                const j: i64 = @intCast(jj);
                bcol[jj - 1] = ob.val(T, @floatFromInt(if (rev) n - j + 1 else j), 0);
                const v = @divTrunc(j * ((j + 1) * j), 2) - @divTrunc((j + 1) * j * (j - 1), 3);
                cc[if (rev) n - jj else jj - 1] = ob.val(T, @floatFromInt(v), 0);
            }
            // B is `bcol` as a column (N) or a row with leading dimension 1 (T)
            const res = mmch(T, t[0], t[1], n, 1, n, ob.val(T, 1, 0), &ab, n, bcol, if (t[1] == 'T') 1 else n, ob.zero(T), &c, n, &ct, &g, &cc, n);
            for (ct, cc) |x, y| try ob.expect(ob.eql(x, ob.cast(R, y)), "SMMCH: wrong");
            try ob.expect(res.err == 0 and !res.fatal, "SMMCH: ratio not 0");
        }
    }
}

test "blat3: SGEMM, DGEMM" {
    try sweep(.gemm, .{ f32, f64 });
}

test "blat3: CGEMM, ZGEMM" {
    try sweep(.gemm, .{ C32, C64 });
}

test "blat3_3m: CGEMM3M, ZGEMM3M" {
    try sweep(.gemm3m, .{ C32, C64 });
}

test "blat3: xSYMM" {
    try sweep(.symm, .{ f32, f64, C32, C64 });
}

test "blat3: CHEMM, ZHEMM" {
    try sweep(.hemm, .{ C32, C64 });
}

test "blat3: xTRMM" {
    try sweep(.trmm, .{ f32, f64, C32, C64 });
}

test "blat3: xTRSM" {
    try sweep(.trsm, .{ f32, f64, C32, C64 });
}

test "blat3: xSYRK" {
    try sweep(.syrk, .{ f32, f64, C32, C64 });
}

test "blat3: CHERK, ZHERK" {
    try sweep(.herk, .{ C32, C64 });
}

test "blat3: xSYR2K" {
    try sweep(.syr2k, .{ f32, f64, C32, C64 });
}

test "blat3: CHER2K, ZHER2K" {
    try sweep(.her2k, .{ C32, C64 });
}
