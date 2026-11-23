//! The level-2 programs `test/[sdcz]blat2.f` with their `.dat` tables and
//! their CBLAS `ctest/` versions, OpenBLAS v0.3.34-317; `gbmv`'s negative `kl`
//! or `ku` is `InvalidK`, and an error on a valid call fails at once.

const std = @import("std");
const ob = @import("__root__.zig");
const blas = @import("../../__root__.zig");
const l2 = @import("../../blas2/core.zig");
const check = @import("../check.zig");
const guard = @import("../guard.zig");

const Float = ob.Float;
const isComplex = ob.isComplex;
const testing = std.testing;

/// `NMAX` and `INCMAX`: the largest order and stride the programs take.
pub const nmax = 65;
pub const incmax = 2;

/// A `.dat` file's sweep: values of `N`, `K`, `INCX` and `INCY`, `ALPHA`
/// and `BETA` (complex as re, im), and the threshold of the test ratio.
pub const Params = struct {
    thresh: f64 = 16,
    n: []const usize,
    k: []const usize = &.{ 0, 1, 2, 4 },
    inc: []const isize = &.{ 1, 2, -1, -2 },
    alpha: []const [2]f64,
    beta: []const [2]f64,
};

/// `test/[sd]blat2.dat`.
pub const real_dat: Params = .{
    .n = &.{ 0, 1, 2, 3, 7, 31, 63 },
    .alpha = &.{ .{ 0, 0 }, .{ 1, 0 }, .{ 0.7, 0 } },
    .beta = &.{ .{ 0, 0 }, .{ 1, 0 }, .{ 0.9, 0 } },
};

/// `test/[cz]blat2.dat`.
pub const complex_dat: Params = .{
    .n = &.{ 0, 1, 2, 3, 7, 31, 63 },
    .alpha = &.{ .{ 0, 0 }, .{ 1, 0 }, .{ 0.7, -0.9 } },
    .beta = &.{ .{ 0, 0 }, .{ 1, 0 }, .{ 1.3, -1.1 } },
};

/// `ctest/[sdcz]in2`: the same but for `N`.
pub const ctest_n = [_]usize{ 0, 1, 2, 3, 5, 9, 63 };

/// The rigs' `N`: the example data file in each program's header.
pub const rig_n = [_]usize{ 0, 1, 2, 3, 5, 9 };

pub fn dat(comptime T: type) Params {
    return if (comptime isComplex(T)) complex_dat else real_dat;
}

pub fn withN(p: Params, n: []const usize) Params {
    var q = p;
    q.n = n;
    return q;
}

/// The routines, by their real names; complex runs `HEMV`, `HBMV`, `HPMV`,
/// `GERU`, `HER`, `HPR`, `HER2`, `HPR2` in their place.
pub const Routine = enum { gemv, gbmv, symv, sbmv, spmv, trmv, tbmv, tpmv, trsv, tbsv, tpsv, ger, gerc, syr, spr, syr2, spr2 };

pub fn name(comptime T: type, r: Routine) []const u8 {
    const c = comptime isComplex(T);
    const s = switch (r) {
        .symv => if (c) "HEMV" else "SYMV",
        .sbmv => if (c) "HBMV" else "SBMV",
        .spmv => if (c) "HPMV" else "SPMV",
        .ger => if (c) "GERU" else "GER",
        .syr => if (c) "HER" else "SYR",
        .spr => if (c) "HPR" else "SPR",
        .syr2 => if (c) "HER2" else "SYR2",
        .spr2 => if (c) "HPR2" else "SPR2",
        inline else => |t| comptime upper(@tagName(t)),
    };
    return s;
}

pub fn upper(comptime s: []const u8) []const u8 {
    var out: [s.len]u8 = undefined;
    for (s, 0..) |ch, i| out[i] = std.ascii.toUpper(ch);
    const final = out;
    return &final;
}

pub fn letter(comptime T: type) u8 {
    return switch (Float(T)) {
        f32 => if (isComplex(T)) 'C' else 'S',
        else => if (isComplex(T)) 'Z' else 'D',
    };
}

const mk = @import("blat2/make.zig");
const schk = @import("blat2/schk.zig");

const mvch = mk.mvch;
const chk1 = schk.chk1;
const chk2 = schk.chk2;
const chk3 = schk.chk3;
const chk4 = schk.chk4;
const chk5 = schk.chk5;
const chk6 = schk.chk6;

fn chk(comptime T: type, comptime r: Routine, p: Params) !void {
    const f = switch (r) {
        .gemv, .gbmv => chk1,
        .symv, .sbmv, .spmv => chk2,
        .trmv, .tbmv, .tpmv, .trsv, .tbsv, .tpsv => chk3,
        .ger, .gerc => chk4,
        .syr, .spr => chk5,
        .syr2, .spr2 => chk6,
    };
    try f(T, r, p);
}
const all_types = check.sd_types;

/// Every type `r` exists for, with the `.dat` table natively, ctest's too,
/// and the header's on the rigs.
fn sweep(comptime r: Routine) !void {
    inline for (all_types) |T| {
        if (r == .gerc and !comptime isComplex(T)) continue;
        if (check.rig()) {
            try chk(T, r, withN(dat(T), &rig_n));
        } else {
            try chk(T, r, dat(T));
            try chk(T, r, withN(dat(T), &ctest_n));
        }
    }
}

/// Every routine at `T` on the `N` of `n`, after the error exits, for
/// `check.eachCore`.
pub fn allOf(comptime T: type, n: []const usize) !void {
    try errorExits(T);
    inline for (comptime std.enums.values(Routine)) |r| {
        if (r != .gerc or comptime isComplex(T)) try chk(T, r, withN(dat(T), n));
    }
}

test "blat2: SMVCH self-check on exact data" {
    inline for (all_types) |T| {
        const R = ob.Ref(T);
        const n = 32;
        var a: [n * n]T = undefined;
        var x: [n]T = undefined;
        var y: [n]T = undefined;
        var yy: [n]T = undefined;
        var yt: [n]R = undefined;
        var g: [n]Float(R) = undefined;
        for (0..n) |jj| {
            const j: i64 = @intCast(jj + 1);
            for (0..n) |ii| {
                const i: i64 = @intCast(ii + 1);
                a[ii + jj * n] = ob.val(T, @floatFromInt(@max(i - j + 1, 0)), 0);
            }
            x[jj] = ob.val(T, @floatFromInt(j), 0);
            y[jj] = ob.zero(T);
            yy[jj] = ob.val(T, @floatFromInt(@divTrunc(j * ((j + 1) * j), 2) - @divTrunc((j + 1) * j * (j - 1), 3)), 0);
        }
        for ([_]struct { u8, isize }{ .{ 'N', 1 }, .{ 'T', -1 } }) |c| {
            const res = mvch(T, c[0], n, n, ob.val(T, 1, 0), &a, n, &x, c[1], ob.zero(T), &y, c[1], &yt, &g, &yy);
            for (yt, yy) |t, e| try ob.expect(ob.eql(t, ob.cast(R, e)), "SMVCH: wrong");
            try ob.expect(res.err == 0 and !res.fatal, "SMVCH: ratio not 0");
        }
    }
}

test "blat2: xGEMV" {
    try sweep(.gemv);
}

test "blat2: xGBMV" {
    try sweep(.gbmv);
}

test "blat2: xSYMV, xHEMV" {
    try sweep(.symv);
}

test "blat2: xSBMV, xHBMV" {
    try sweep(.sbmv);
}

test "blat2: xSPMV, xHPMV" {
    try sweep(.spmv);
}

test "blat2: xTRMV" {
    try sweep(.trmv);
}

test "blat2: xTBMV" {
    try sweep(.tbmv);
}

test "blat2: xTPMV" {
    try sweep(.tpmv);
}

test "blat2: xTRSV" {
    try sweep(.trsv);
}

test "blat2: xTBSV" {
    try sweep(.tbsv);
}

test "blat2: xTPSV" {
    try sweep(.tpsv);
}

test "blat2: xGER, xGERU" {
    try sweep(.ger);
}

test "blat2: xGERC" {
    try sweep(.gerc);
}

test "blat2: xSYR, xHER" {
    try sweep(.syr);
}

test "blat2: xSPR, xHPR" {
    try sweep(.spr);
}

test "blat2: xSYR2, xHER2" {
    try sweep(.syr2);
}

test "blat2: xSPR2, xHPR2" {
    try sweep(.spr2);
}

// ---- SCHKE ------------------------------------------------------------------

/// `SCHKE`/`CCHKE`: each call has one bad argument, and must report it.
fn errorExits(comptime T: type) !void {
    const al = testing.allocator;
    const cx = comptime isComplex(T);
    var a: [2]T = @splat(ob.zero(T));
    var x: [2]T = @splat(ob.zero(T));
    var y: [2]T = @splat(ob.zero(T));
    const al0 = ob.zero(T);
    const ra: Float(T) = 0;
    const E = ob.expectError;
    // gemv: a bad letter is one real types do not take
    if (!cx) try E(error.InvalidTrans, l2.gemv(al, ob.io, T, .o, 0, 0, al0, &a, 1, &x, 1, al0, &y, 1), "GEMV INFOT=1");
    try E(error.InvalidM, l2.gemv(al, ob.io, T, .n, -1, 0, al0, &a, 1, &x, 1, al0, &y, 1), "GEMV INFOT=2");
    try E(error.InvalidN, l2.gemv(al, ob.io, T, .n, 0, -1, al0, &a, 1, &x, 1, al0, &y, 1), "GEMV INFOT=3");
    try E(error.InvalidLda, l2.gemv(al, ob.io, T, .n, 2, 0, al0, &a, 1, &x, 1, al0, &y, 1), "GEMV INFOT=6");
    try E(error.InvalidIncx, l2.gemv(al, ob.io, T, .n, 0, 0, al0, &a, 1, &x, 0, al0, &y, 1), "GEMV INFOT=8");
    try E(error.InvalidIncy, l2.gemv(al, ob.io, T, .n, 0, 0, al0, &a, 1, &x, 1, al0, &y, 0), "GEMV INFOT=11");
    // gbmv
    if (!cx) try E(error.InvalidTrans, l2.gbmv(al, ob.io, T, .o, 0, 0, 0, 0, al0, &a, 1, &x, 1, al0, &y, 1), "GBMV INFOT=1");
    try E(error.InvalidM, l2.gbmv(al, ob.io, T, .n, -1, 0, 0, 0, al0, &a, 1, &x, 1, al0, &y, 1), "GBMV INFOT=2");
    try E(error.InvalidN, l2.gbmv(al, ob.io, T, .n, 0, -1, 0, 0, al0, &a, 1, &x, 1, al0, &y, 1), "GBMV INFOT=3");
    try E(error.InvalidK, l2.gbmv(al, ob.io, T, .n, 0, 0, -1, 0, al0, &a, 1, &x, 1, al0, &y, 1), "GBMV INFOT=4");
    try E(error.InvalidK, l2.gbmv(al, ob.io, T, .n, 2, 0, 0, -1, al0, &a, 1, &x, 1, al0, &y, 1), "GBMV INFOT=5");
    try E(error.InvalidLda, l2.gbmv(al, ob.io, T, .n, 0, 0, 1, 0, al0, &a, 1, &x, 1, al0, &y, 1), "GBMV INFOT=8");
    try E(error.InvalidIncx, l2.gbmv(al, ob.io, T, .n, 0, 0, 0, 0, al0, &a, 1, &x, 0, al0, &y, 1), "GBMV INFOT=10");
    try E(error.InvalidIncy, l2.gbmv(al, ob.io, T, .n, 0, 0, 0, 0, al0, &a, 1, &x, 1, al0, &y, 0), "GBMV INFOT=13");
    // symv, hemv
    if (cx) {
        try E(error.InvalidN, l2.hemv(al, ob.io, T, .u, -1, al0, &a, 1, &x, 1, al0, &y, 1), "HEMV INFOT=2");
        try E(error.InvalidLda, l2.hemv(al, ob.io, T, .u, 2, al0, &a, 1, &x, 1, al0, &y, 1), "HEMV INFOT=5");
        try E(error.InvalidIncx, l2.hemv(al, ob.io, T, .u, 0, al0, &a, 1, &x, 0, al0, &y, 1), "HEMV INFOT=7");
        try E(error.InvalidIncy, l2.hemv(al, ob.io, T, .u, 0, al0, &a, 1, &x, 1, al0, &y, 0), "HEMV INFOT=10");
    } else {
        try E(error.InvalidUplo, l2.symv(al, ob.io, T, .v, 0, al0, &a, 1, &x, 1, al0, &y, 1), "SYMV INFOT=1");
        try E(error.InvalidN, l2.symv(al, ob.io, T, .u, -1, al0, &a, 1, &x, 1, al0, &y, 1), "SYMV INFOT=2");
        try E(error.InvalidLda, l2.symv(al, ob.io, T, .u, 2, al0, &a, 1, &x, 1, al0, &y, 1), "SYMV INFOT=5");
        try E(error.InvalidIncx, l2.symv(al, ob.io, T, .u, 0, al0, &a, 1, &x, 0, al0, &y, 1), "SYMV INFOT=7");
        try E(error.InvalidIncy, l2.symv(al, ob.io, T, .u, 0, al0, &a, 1, &x, 1, al0, &y, 0), "SYMV INFOT=10");
    }
    // sbmv, hbmv
    if (!cx) try E(error.InvalidUplo, l2.sbmv(al, ob.io, T, .v, 0, 0, al0, &a, 1, &x, 1, al0, &y, 1), "SBMV INFOT=1");
    if (cx) {
        try E(error.InvalidN, l2.hbmv(al, ob.io, T, .u, -1, 0, al0, &a, 1, &x, 1, al0, &y, 1), "HBMV INFOT=2");
        try E(error.InvalidK, l2.hbmv(al, ob.io, T, .u, 0, -1, al0, &a, 1, &x, 1, al0, &y, 1), "HBMV INFOT=3");
        try E(error.InvalidLda, l2.hbmv(al, ob.io, T, .u, 0, 1, al0, &a, 1, &x, 1, al0, &y, 1), "HBMV INFOT=6");
        try E(error.InvalidIncx, l2.hbmv(al, ob.io, T, .u, 0, 0, al0, &a, 1, &x, 0, al0, &y, 1), "HBMV INFOT=8");
        try E(error.InvalidIncy, l2.hbmv(al, ob.io, T, .u, 0, 0, al0, &a, 1, &x, 1, al0, &y, 0), "HBMV INFOT=11");
        // `HEMVREV` letters reach `hpmv` and the four below only from
        // row-major CBLAS upstream
        try E(error.InvalidUplo, l2.hpmv(al, ob.io, T, .v, 0, al0, &a, &x, 1, al0, &y, 1), "HPMV INFOT=1");
        try E(error.InvalidN, l2.hpmv(al, ob.io, T, .u, -1, al0, &a, &x, 1, al0, &y, 1), "HPMV INFOT=2");
        try E(error.InvalidIncx, l2.hpmv(al, ob.io, T, .u, 0, al0, &a, &x, 0, al0, &y, 1), "HPMV INFOT=6");
        try E(error.InvalidIncy, l2.hpmv(al, ob.io, T, .u, 0, al0, &a, &x, 1, al0, &y, 0), "HPMV INFOT=9");
    } else {
        try E(error.InvalidN, l2.sbmv(al, ob.io, T, .u, -1, 0, al0, &a, 1, &x, 1, al0, &y, 1), "SBMV INFOT=2");
        try E(error.InvalidK, l2.sbmv(al, ob.io, T, .u, 0, -1, al0, &a, 1, &x, 1, al0, &y, 1), "SBMV INFOT=3");
        try E(error.InvalidLda, l2.sbmv(al, ob.io, T, .u, 0, 1, al0, &a, 1, &x, 1, al0, &y, 1), "SBMV INFOT=6");
        try E(error.InvalidIncx, l2.sbmv(al, ob.io, T, .u, 0, 0, al0, &a, 1, &x, 0, al0, &y, 1), "SBMV INFOT=8");
        try E(error.InvalidIncy, l2.sbmv(al, ob.io, T, .u, 0, 0, al0, &a, 1, &x, 1, al0, &y, 0), "SBMV INFOT=11");
        try E(error.InvalidN, l2.spmv(al, ob.io, T, .u, -1, al0, &a, &x, 1, al0, &y, 1), "SPMV INFOT=2");
        try E(error.InvalidIncx, l2.spmv(al, ob.io, T, .u, 0, al0, &a, &x, 0, al0, &y, 1), "SPMV INFOT=6");
        try E(error.InvalidIncy, l2.spmv(al, ob.io, T, .u, 0, al0, &a, &x, 1, al0, &y, 0), "SPMV INFOT=9");
    }
    // the triangular six: a Hermitian `uplo`, an `x`-conjugating `trans`;
    // every `Diag` is valid
    try E(error.InvalidUplo, l2.trmv(al, ob.io, T, .v, .n, .n, 0, &a, 1, &x, 1), "TRMV INFOT=1");
    try E(error.InvalidTrans, l2.trmv(al, ob.io, T, .u, .o, .n, 0, &a, 1, &x, 1), "TRMV INFOT=2");
    inline for (.{ l2.trmv, l2.trsv }, .{ "TRMV", "TRSV" }) |f, nm| {
        try E(error.InvalidN, f(al, ob.io, T, .u, .n, .n, -1, &a, 1, &x, 1), nm ++ " INFOT=4");
        try E(error.InvalidLda, f(al, ob.io, T, .u, .n, .n, 2, &a, 1, &x, 1), nm ++ " INFOT=6");
        try E(error.InvalidIncx, f(al, ob.io, T, .u, .n, .n, 0, &a, 1, &x, 0), nm ++ " INFOT=8");
    }
    inline for (.{ l2.tbmv, l2.tbsv }, .{ "TBMV", "TBSV" }) |f, nm| {
        try E(error.InvalidN, f(al, ob.io, T, .u, .n, .n, -1, 0, &a, 1, &x, 1), nm ++ " INFOT=4");
        try E(error.InvalidK, f(al, ob.io, T, .u, .n, .n, 0, -1, &a, 1, &x, 1), nm ++ " INFOT=5");
        try E(error.InvalidLda, f(al, ob.io, T, .u, .n, .n, 0, 1, &a, 1, &x, 1), nm ++ " INFOT=7");
        try E(error.InvalidIncx, f(al, ob.io, T, .u, .n, .n, 0, 0, &a, 1, &x, 0), nm ++ " INFOT=9");
    }
    inline for (.{ l2.tpmv, l2.tpsv }, .{ "TPMV", "TPSV" }) |f, nm| {
        try E(error.InvalidN, f(al, ob.io, T, .u, .n, .n, -1, &a, &x, 1), nm ++ " INFOT=4");
        try E(error.InvalidIncx, f(al, ob.io, T, .u, .n, .n, 0, &a, &x, 0), nm ++ " INFOT=7");
    }
    // ger, geru, gerc
    inline for (if (cx) .{ l2.ger, l2.gerc } else .{l2.ger}) |f| {
        try E(error.InvalidM, f(al, ob.io, T, -1, 0, al0, &x, 1, &y, 1, &a, 1), "GER INFOT=1");
        try E(error.InvalidN, f(al, ob.io, T, 0, -1, al0, &x, 1, &y, 1, &a, 1), "GER INFOT=2");
        try E(error.InvalidIncx, f(al, ob.io, T, 0, 0, al0, &x, 0, &y, 1, &a, 1), "GER INFOT=5");
        try E(error.InvalidIncy, f(al, ob.io, T, 0, 0, al0, &x, 1, &y, 0, &a, 1), "GER INFOT=7");
        try E(error.InvalidLda, f(al, ob.io, T, 2, 0, al0, &x, 1, &y, 1, &a, 1), "GER INFOT=9");
    }
    // syr, spr, syr2, spr2; her, hpr, her2, hpr2
    if (cx) {
        try E(error.InvalidUplo, l2.her(al, ob.io, T, .m, 0, ra, &x, 1, &a, 1), "HER INFOT=1");
        try E(error.InvalidN, l2.her(al, ob.io, T, .u, -1, ra, &x, 1, &a, 1), "HER INFOT=2");
        try E(error.InvalidIncx, l2.her(al, ob.io, T, .u, 0, ra, &x, 0, &a, 1), "HER INFOT=5");
        try E(error.InvalidLda, l2.her(al, ob.io, T, .u, 2, ra, &x, 1, &a, 1), "HER INFOT=7");
        try E(error.InvalidUplo, l2.hpr(al, ob.io, T, .v, 0, ra, &x, 1, &a), "HPR INFOT=1");
        try E(error.InvalidN, l2.hpr(al, ob.io, T, .u, -1, ra, &x, 1, &a), "HPR INFOT=2");
        try E(error.InvalidIncx, l2.hpr(al, ob.io, T, .u, 0, ra, &x, 0, &a), "HPR INFOT=5");
        try E(error.InvalidUplo, l2.her2(al, ob.io, T, .m, 0, al0, &x, 1, &y, 1, &a, 1), "HER2 INFOT=1");
        try E(error.InvalidN, l2.her2(al, ob.io, T, .u, -1, al0, &x, 1, &y, 1, &a, 1), "HER2 INFOT=2");
        try E(error.InvalidIncx, l2.her2(al, ob.io, T, .u, 0, al0, &x, 0, &y, 1, &a, 1), "HER2 INFOT=5");
        try E(error.InvalidIncy, l2.her2(al, ob.io, T, .u, 0, al0, &x, 1, &y, 0, &a, 1), "HER2 INFOT=7");
        try E(error.InvalidLda, l2.her2(al, ob.io, T, .u, 2, al0, &x, 1, &y, 1, &a, 1), "HER2 INFOT=9");
        try E(error.InvalidUplo, l2.hpr2(al, ob.io, T, .v, 0, al0, &x, 1, &y, 1, &a), "HPR2 INFOT=1");
        try E(error.InvalidN, l2.hpr2(al, ob.io, T, .u, -1, al0, &x, 1, &y, 1, &a), "HPR2 INFOT=2");
        try E(error.InvalidIncx, l2.hpr2(al, ob.io, T, .u, 0, al0, &x, 0, &y, 1, &a), "HPR2 INFOT=5");
        try E(error.InvalidIncy, l2.hpr2(al, ob.io, T, .u, 0, al0, &x, 1, &y, 0, &a), "HPR2 INFOT=7");
    } else {
        try E(error.InvalidN, l2.syr(al, ob.io, T, .u, -1, al0, &x, 1, &a, 1), "SYR INFOT=2");
        try E(error.InvalidIncx, l2.syr(al, ob.io, T, .u, 0, al0, &x, 0, &a, 1), "SYR INFOT=5");
        try E(error.InvalidLda, l2.syr(al, ob.io, T, .u, 2, al0, &x, 1, &a, 1), "SYR INFOT=7");
        try E(error.InvalidN, l2.spr(al, ob.io, T, .u, -1, al0, &x, 1, &a), "SPR INFOT=2");
        try E(error.InvalidIncx, l2.spr(al, ob.io, T, .u, 0, al0, &x, 0, &a), "SPR INFOT=5");
        try E(error.InvalidN, l2.syr2(al, ob.io, T, .u, -1, al0, &x, 1, &y, 1, &a, 1), "SYR2 INFOT=2");
        try E(error.InvalidIncx, l2.syr2(al, ob.io, T, .u, 0, al0, &x, 0, &y, 1, &a, 1), "SYR2 INFOT=5");
        try E(error.InvalidIncy, l2.syr2(al, ob.io, T, .u, 0, al0, &x, 1, &y, 0, &a, 1), "SYR2 INFOT=7");
        try E(error.InvalidLda, l2.syr2(al, ob.io, T, .u, 2, al0, &x, 1, &y, 1, &a, 1), "SYR2 INFOT=9");
        try E(error.InvalidN, l2.spr2(al, ob.io, T, .u, -1, al0, &x, 1, &y, 1, &a), "SPR2 INFOT=2");
        try E(error.InvalidIncx, l2.spr2(al, ob.io, T, .u, 0, al0, &x, 0, &y, 1, &a), "SPR2 INFOT=5");
        try E(error.InvalidIncy, l2.spr2(al, ob.io, T, .u, 0, al0, &x, 1, &y, 0, &a), "SPR2 INFOT=7");
    }
}

test "blat2: error exits" {
    inline for (all_types) |T| errorExits(T) catch |e| {
        ob.say("{c}\n", .{letter(T)});
        return e;
    };
}
