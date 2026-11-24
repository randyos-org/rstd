//! The edge sweep's cases: where a pattern goes, the routines, and their shapes.

const std = @import("std");
const util = @import("../../../util.zig");
const check = @import("../../check.zig");
const guard = @import("../../guard.zig");
const pattern = @import("../../pattern.zig");
const mx = @import("matrix.zig");

const isComplex = check.isComplex;
const Pattern = pattern.Pattern;
const Kind = mx.Kind;
const Layout = mx.Layout;
const Shape = mx.Shape;

/// Where a case's pattern goes.
pub const Where = enum { a, x, y, none };

pub const Case = struct {
    m: usize,
    n: usize,
    /// A band's widths: `kl`, `ku`, or a triangle's `k` in both.
    kl: usize,
    ku: usize,
    inc: [2]isize,
    p: Pattern,
    where: Where,
    /// A second operand's pattern.
    also: ?struct { where: Where, p: Pattern } = null,
    /// The matrix's rows of padding (`Shape.pad`).
    pad: usize = 1,

    /// The pattern operand `w` takes, if any.
    pub fn on(c: Case, w: Where) ?Pattern {
        if (c.where == w) return c.p;
        if (c.also) |a| if (a.where == w) return a.p;
        return null;
    }
};

pub const Routine = enum {
    gemv,
    gbmv,
    symv,
    hemv,
    sbmv,
    hbmv,
    spmv,
    hpmv,
    trmv,
    tbmv,
    tpmv,
    trsv,
    tbsv,
    tpsv,
    ger,
    gerc,
    syr,
    her,
    spr,
    hpr,
    syr2,
    her2,
    spr2,
    hpr2,

    pub fn complexOnly(r: Routine) bool {
        return switch (r) {
            .hemv, .hbmv, .hpmv, .gerc, .her, .hpr, .her2, .hpr2 => true,
            else => false,
        };
    }

    pub const Family = enum { mv, tr, rank };

    pub fn family(r: Routine) Family {
        return switch (r) {
            .gemv, .gbmv, .symv, .hemv, .sbmv, .hbmv, .spmv, .hpmv => .mv,
            .trmv, .tbmv, .tpmv, .trsv, .tbsv, .tpsv => .tr,
            else => .rank,
        };
    }

    pub fn layout(r: Routine) Layout {
        return switch (r) {
            .gbmv, .sbmv, .hbmv, .tbmv, .tbsv => .band,
            .spmv, .hpmv, .tpmv, .tpsv, .spr, .hpr, .spr2, .hpr2 => .packed_,
            else => .full,
        };
    }

    pub fn solve(r: Routine) bool {
        return r == .trsv or r == .tbsv or r == .tpsv;
    }

    /// Variants: `gemv`'s and `gbmv`'s transposes, the triangles (with
    /// `hemv`'s and `hbmv`'s reversed Hermitian ones), and the triangular
    /// routines' triangle, transpose and diagonal (`v / 8`, `v / 2 % 4`,
    /// `v % 2`).
    pub fn variants(r: Routine, complex: bool) usize {
        return switch (r) {
            .gemv, .gbmv => if (complex) 8 else 2,
            .hemv, .hbmv => 4,
            .trmv, .tbmv, .tpmv, .trsv, .tbsv, .tpsv => 16,
            .ger, .gerc => 1,
            else => 2,
        };
    }
};

/// `gemv`'s `TRANS` letters by variant: transposed, `A` conjugated, `x`
/// conjugated.
pub const trans_flags = [8][3]bool{
    .{ false, false, false }, .{ true, false, false }, .{ false, true, false }, .{ true, true, false },
    .{ false, false, true },  .{ true, false, true },  .{ false, true, true },  .{ true, true, true },
};

pub fn shapeOf(r: Routine, v: usize, c: Case) Shape {
    const herm_kind: Kind = if (v >= 2) .hermitian_rev else .hermitian;
    var s: Shape = switch (r) {
        .gemv, .ger, .gerc => .{ .kind = .general, .m = c.m, .n = c.n },
        .gbmv => .{ .kind = .general, .layout = .band, .m = c.m, .n = c.n, .kl = c.kl, .ku = c.ku },
        .symv, .sbmv, .spmv, .syr, .spr, .syr2, .spr2 => .{ .kind = .symmetric, .layout = r.layout(), .m = c.n, .n = c.n, .lower = v == 1, .kl = if (v == 1) c.kl else 0, .ku = if (v == 1) 0 else c.ku },
        .hemv, .hbmv, .hpmv, .her, .hpr, .her2, .hpr2 => .{ .kind = herm_kind, .layout = r.layout(), .m = c.n, .n = c.n, .lower = v % 2 == 1, .kl = if (v % 2 == 1) c.kl else 0, .ku = if (v % 2 == 1) 0 else c.ku },
        .trmv, .tbmv, .tpmv, .trsv, .tbsv, .tpsv => .{ .kind = .triangular, .layout = r.layout(), .m = c.n, .n = c.n, .lower = v / 8 == 1, .unit = v % 2 == 1, .kl = if (v / 8 == 1) c.kl else 0, .ku = if (v / 8 == 1) 0 else c.ku, .solve = r.solve() },
    };
    s.pad = c.pad;
    return s;
}

/// Scratch from the heap holds NaN, so a driver that reads what it did not
/// write shows.
pub var poison: guard.Poison = .{ .child = std.testing.allocator };

pub const val = check.val;

pub fn zeroW(v: anytype) bool {
    return if (comptime isComplex(@TypeOf(v))) v.re == 0 and v.im == 0 else v == 0;
}

pub const isZero = util.isZero;
