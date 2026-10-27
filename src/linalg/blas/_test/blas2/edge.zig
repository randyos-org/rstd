//! Edge cases through every level-2 routine, on every core this machine runs,
//! through the interfaces: the special values of `check.Pattern` in one
//! operand at a time (the matrix, `x` or `y`), against an `f128` reference.
//! What a routine must not read holds NaN: the other triangle, a unit
//! diagonal, a Hermitian diagonal's imaginary part, band corners, padding
//! rows, the gaps of strided vectors, and `y` for a zero `beta` or the matrix
//! and the vectors for a zero `alpha`; where it must not write, the NaN must
//! stay.
//!
//! Not propagation, on purpose, as the reference has it: a zero `beta`
//! overwrites `y`, a zero `alpha` returns once `y` is scaled, a Hermitian
//! update stores a zero imaginary part on the diagonal, entries outside a
//! band or triangle are never multiplied, nor a unit diagonal, and the rank
//! updates pass over a column whose factors are zero (`y` for `ger`; `x`,
//! and `y` for rank 2, otherwise). Complex infinities may meet in another
//! order than the reference's, giving NaN for infinity
//! (`check.expectWideAnyOrder`).
//!
//! A solve divides by the diagonal, whose magnitudes are not its cases: its
//! patterns go in the off-diagonal entries (none huge, which overflow its
//! sums) and the right-hand side, and the diagonal dominates.

const std = @import("std");
const blas = @import("../../__root__.zig");
const Vtable = @import("../../Vtable.zig");
const check = @import("../blas1/check.zig");
const guard = @import("guard.zig");
const portable = @import("../../portable.zig");
const accel = @import("../../apple_accel.zig");

const Allocator = std.mem.Allocator;
const C = std.math.Complex;
const Float = check.FloatOf;
const isComplex = check.isComplex;
const Wide = check.Wide;
const Pattern = check.Pattern;
const Vector = check.Vector;

const Kind = enum { general, symmetric, hermitian, hermitian_rev, triangular };
const Layout = enum { full, band, packed_ };

/// What a matrix operand is and how it is stored.
const Shape = struct {
    kind: Kind,
    layout: Layout = .full,
    m: usize,
    n: usize,
    /// Sub- and superdiagonals of a band; a triangle's other one is 0.
    kl: usize = 0,
    ku: usize = 0,
    lower: bool = false,
    unit: bool = false,
    /// A solve's: the diagonal dominant and left out of the pattern.
    solve: bool = false,

    /// One row of padding below a full or band matrix.
    fn ld(s: Shape) usize {
        return switch (s.layout) {
            .full => s.m + 1,
            .band => s.kl + s.ku + 2,
            .packed_ => 1,
        };
    }

    fn len(s: Shape) usize {
        return if (s.layout == .packed_) s.n * (s.n + 1) / 2 else s.ld() * s.n;
    }

    /// Where `(i, j)` is stored, if it is.
    fn stored(s: Shape, i: usize, j: usize) ?usize {
        if (s.kind == .general) return switch (s.layout) {
            .full => i + j * s.ld(),
            .band => if (i + s.ku >= j and i <= j + s.kl) s.ku + i - j + j * s.ld() else null,
            .packed_ => unreachable,
        };
        if (if (s.lower) i < j else i > j) return null;
        const d = if (s.lower) i - j else j - i;
        return switch (s.layout) {
            .full => i + j * s.ld(),
            .band => if (d > s.kl + s.ku) null else (if (s.lower) d else s.ku - d) + j * s.ld(),
            .packed_ => if (s.lower) j * (2 * s.n - j + 1) / 2 + d else j * (j + 1) / 2 + i,
        };
    }

    /// Whether the routine reads entry `(i, j)`: stored, and not a unit diagonal.
    fn read(s: Shape, i: usize, j: usize) ?usize {
        if (i == j and s.unit) return null;
        return s.stored(i, j);
    }

    /// The entries a pattern goes through, in `init`'s order.
    fn entries(s: Shape) usize {
        var k: usize = 0;
        for (0..s.n) |j| for (0..s.m) |i| {
            if (s.read(i, j) != null and !(s.solve and i == j)) k += 1;
        };
        return k;
    }
};

fn Matrix(comptime T: type) type {
    return struct {
        const Self = @This();
        s: Shape,
        buf: []T,
        was: []T,

        /// Random entries where the routine reads, shaped by `p` in order;
        /// NaN everywhere else, and in a Hermitian diagonal's imaginary part.
        fn init(al: Allocator, s: Shape, seed: u64, p: ?Pattern) !Self {
            const F = Float(T);
            const buf = try al.alloc(T, s.len() + check.tail);
            errdefer al.free(buf);
            @memset(std.mem.sliceAsBytes(buf), 0xff);
            const vals = try al.alloc(T, s.len());
            defer al.free(vals);
            check.fill(T, vals, seed);
            const sc: F = if (s.solve) 1 / @sqrt(@as(F, @floatFromInt(s.n + 1))) else 1;
            var k: usize = 0;
            for (0..s.n) |j| for (0..s.m) |i| {
                const idx = s.read(i, j) orelse continue;
                var v = vals[idx];
                if (s.solve and i == j) {
                    v = if (comptime isComplex(T)) .init(3 + @abs(v.re), v.im) else 3 + @abs(v);
                } else {
                    if (s.solve) v = if (comptime isComplex(T)) .init(v.re * sc, v.im * sc) else v * sc;
                    if (p) |q| v = check.shape(T, v, k, q);
                    k += 1;
                }
                if (comptime isComplex(T)) if (i == j and (s.kind == .hermitian or s.kind == .hermitian_rev)) {
                    buf[idx].re = v.re;
                    continue;
                };
                buf[idx] = v;
            };
            return .{ .s = s, .buf = buf, .was = try al.dupe(T, buf) };
        }

        fn deinit(self: *Self, al: Allocator) void {
            al.free(self.was);
            al.free(self.buf);
        }

        fn ld(self: *const Self) isize {
            return @intCast(self.s.ld());
        }

        /// NaN in every entry: the routine must not read them.
        fn poison(self: *Self) void {
            @memset(std.mem.sliceAsBytes(self.buf), 0xff);
            @memcpy(self.was, self.buf);
        }

        /// Entry `(i, j)` of the matrix the storage stands for, as it was;
        /// null for a zero outside the band or triangle, which the reference
        /// never multiplies.
        fn at(self: *const Self, i: usize, j: usize) ?Wide(T) {
            const s = self.s;
            const z: ?Wide(T) = null;
            switch (s.kind) {
                .general => return if (s.stored(i, j)) |k| check.wide(T, self.was[k]) else z,
                .triangular => {
                    const k = s.stored(i, j) orelse return z;
                    if (i == j and s.unit) return check.wide(T, if (comptime isComplex(T)) T.init(1, 0) else 1);
                    return check.wide(T, self.was[k]);
                },
                .symmetric => {
                    const k = s.stored(i, j) orelse s.stored(j, i) orelse return z;
                    return check.wide(T, self.was[k]);
                },
                .hermitian, .hermitian_rev => {
                    if (comptime !isComplex(T)) unreachable;
                    if (i == j) return .init(self.was[s.stored(i, i).?].re, 0);
                    const rev = s.kind == .hermitian_rev;
                    if (s.stored(i, j)) |k| return if (rev) check.conjW(check.wide(T, self.was[k])) else check.wide(T, self.was[k]);
                    const k = s.stored(j, i) orelse return z;
                    return if (rev) check.wide(T, self.was[k]) else check.conjW(check.wide(T, self.was[k]));
                },
            }
        }

        fn unchanged(self: *const Self) !void {
            for (self.buf, self.was, 0..) |a, b, k| if (!check.same(T, a, b, true)) {
                check.say("matrix storage {d} changed: ", .{k});
                return error.TestUnexpectedResult;
            };
        }
    };
}

/// Where a case's pattern goes.
const Where = enum { a, x, y, none };

const Case = struct {
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

    /// The pattern operand `w` takes, if any.
    fn on(c: Case, w: Where) ?Pattern {
        if (c.where == w) return c.p;
        if (c.also) |a| if (a.where == w) return a.p;
        return null;
    }
};

const Routine = enum {
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

    fn complexOnly(r: Routine) bool {
        return switch (r) {
            .hemv, .hbmv, .hpmv, .gerc, .her, .hpr, .her2, .hpr2 => true,
            else => false,
        };
    }

    const Family = enum { mv, tr, rank };

    fn family(r: Routine) Family {
        return switch (r) {
            .gemv, .gbmv, .symv, .hemv, .sbmv, .hbmv, .spmv, .hpmv => .mv,
            .trmv, .tbmv, .tpmv, .trsv, .tbsv, .tpsv => .tr,
            else => .rank,
        };
    }

    fn layout(r: Routine) Layout {
        return switch (r) {
            .gbmv, .sbmv, .hbmv, .tbmv, .tbsv => .band,
            .spmv, .hpmv, .tpmv, .tpsv, .spr, .hpr, .spr2, .hpr2 => .packed_,
            else => .full,
        };
    }

    fn solve(r: Routine) bool {
        return r == .trsv or r == .tbsv or r == .tpsv;
    }

    /// Variants: `gemv`'s and `gbmv`'s transposes, the triangles (with the
    /// reversed Hermitian ones), and the triangular routines' triangle,
    /// transpose and diagonal (`v / 8`, `v / 2 % 4`, `v % 2`).
    fn variants(r: Routine, complex: bool) usize {
        return switch (r) {
            .gemv, .gbmv => if (complex) 8 else 2,
            .hemv, .hbmv, .hpmv, .her, .hpr, .her2, .hpr2 => 4,
            .trmv, .tbmv, .tpmv, .trsv, .tbsv, .tpsv => 16,
            .ger, .gerc => 1,
            else => 2,
        };
    }
};

/// `gemv`'s `TRANS` letters by variant: transposed, `A` conjugated, `x`
/// conjugated.
const trans_flags = [8][3]bool{
    .{ false, false, false }, .{ true, false, false }, .{ false, true, false }, .{ true, true, false },
    .{ false, false, true },  .{ true, false, true },  .{ false, true, true },  .{ true, true, true },
};

fn shapeOf(r: Routine, v: usize, c: Case) Shape {
    const herm_kind: Kind = if (v >= 2) .hermitian_rev else .hermitian;
    return switch (r) {
        .gemv, .ger, .gerc => .{ .kind = .general, .m = c.m, .n = c.n },
        .gbmv => .{ .kind = .general, .layout = .band, .m = c.m, .n = c.n, .kl = c.kl, .ku = c.ku },
        .symv, .sbmv, .spmv, .syr, .spr, .syr2, .spr2 => .{ .kind = .symmetric, .layout = r.layout(), .m = c.n, .n = c.n, .lower = v == 1, .kl = if (v == 1) c.kl else 0, .ku = if (v == 1) 0 else c.ku },
        .hemv, .hbmv, .hpmv, .her, .hpr, .her2, .hpr2 => .{ .kind = herm_kind, .layout = r.layout(), .m = c.n, .n = c.n, .lower = v % 2 == 1, .kl = if (v % 2 == 1) c.kl else 0, .ku = if (v % 2 == 1) 0 else c.ku },
        .trmv, .tbmv, .tpmv, .trsv, .tbsv, .tpsv => .{ .kind = .triangular, .layout = r.layout(), .m = c.n, .n = c.n, .lower = v / 8 == 1, .unit = v % 2 == 1, .kl = if (v / 8 == 1) c.kl else 0, .ku = if (v / 8 == 1) 0 else c.ku, .solve = r.solve() },
    };
}

/// Scratch from the heap holds NaN, so a driver that reads what it did not
/// write shows.
var poison: guard.Poison = .{ .child = std.testing.allocator };

fn val(comptime T: type, re: f64, im: f64) T {
    return if (comptime isComplex(T)) .init(@floatCast(re), @floatCast(im)) else @floatCast(re);
}

fn zeroW(v: anytype) bool {
    return if (comptime isComplex(@TypeOf(v))) v.re == 0 and v.im == 0 else v == 0;
}

fn isZero(comptime T: type, v: T) bool {
    return if (comptime isComplex(T)) v.re == 0 and v.im == 0 else v == 0;
}

/// `y := alpha*op(A)*x + beta*y` for the matrix-vector products.
fn mvRun(comptime T: type, comptime r: Routine, v: usize, alpha: T, beta: T, c: Case) !void {
    const al = std.testing.allocator;
    const s = shapeOf(r, v, c);
    const t = if (r == .gemv or r == .gbmv) trans_flags[v] else [3]bool{ false, false, false };
    const rows = if (t[0]) s.n else s.m;
    const cols = if (t[0]) s.m else s.n;
    var a: Matrix(T) = try .init(al, s, 7 * c.m + c.n, c.on(.a));
    defer a.deinit(al);
    var x: Vector(T) = try .init(al, cols, c.inc[0], cols + 3, c.on(.x));
    defer x.deinit(al);
    var y: Vector(T) = try .init(al, rows, c.inc[1], rows + 5, c.on(.y));
    defer y.deinit(al);
    if (isZero(T, beta)) y.poison();
    if (isZero(T, alpha)) {
        a.poison();
        x.poison();
    }
    const pa = poison.allocator();
    const m: isize = @intCast(s.m);
    const n: isize = @intCast(s.n);
    const ix = c.inc[0];
    const iy = c.inc[1];
    switch (r) {
        .gemv => try blas.gemv(pa, check.io, T, @fromBackingInt(@intCast(v)), m, n, alpha, a.buf.ptr, a.ld(), x.base(), ix, beta, y.base(), iy),
        .gbmv => try blas.gbmv(pa, check.io, T, @fromBackingInt(@intCast(v)), m, n, @intCast(s.kl), @intCast(s.ku), alpha, a.buf.ptr, a.ld(), x.base(), ix, beta, y.base(), iy),
        .symv => try blas.symv(pa, check.io, T, if (s.lower) .l else .u, n, alpha, a.buf.ptr, a.ld(), x.base(), ix, beta, y.base(), iy),
        .hemv => if (comptime isComplex(T)) try blas.hemv(pa, check.io, T, @fromBackingInt(@intCast(v)), n, alpha, a.buf.ptr, a.ld(), x.base(), ix, beta, y.base(), iy),
        .sbmv => try blas.sbmv(pa, check.io, T, if (s.lower) .l else .u, n, @intCast(s.kl + s.ku), alpha, a.buf.ptr, a.ld(), x.base(), ix, beta, y.base(), iy),
        .hbmv => if (comptime isComplex(T)) try blas.hbmv(pa, check.io, T, @fromBackingInt(@intCast(v)), n, @intCast(s.kl + s.ku), alpha, a.buf.ptr, a.ld(), x.base(), ix, beta, y.base(), iy),
        .spmv => try blas.spmv(pa, check.io, T, if (s.lower) .l else .u, n, alpha, a.buf.ptr, x.base(), ix, beta, y.base(), iy),
        .hpmv => if (comptime isComplex(T)) try blas.hpmv(pa, check.io, T, @fromBackingInt(@intCast(v)), n, alpha, a.buf.ptr, x.base(), ix, beta, y.base(), iy),
        else => unreachable,
    }
    try a.unchanged();
    try x.unchanged();
    const wa = check.wide(T, alpha);
    const wb = check.wide(T, beta);
    for (0..rows) |i| {
        var acc = check.wide(T, val(T, 0, 0));
        var mag: f128 = 0;
        var terms: usize = 0;
        if (!isZero(T, alpha)) for (0..cols) |j| {
            var e = (if (t[0]) a.at(j, i) else a.at(i, j)) orelse continue;
            terms += 1;
            if (t[1]) e = check.conjW(e);
            const xv = if (t[2]) check.conjW(x.w(j)) else x.w(j);
            acc = check.addW(acc, check.mulW(e, xv));
            mag += check.mag(e) * check.mag(xv);
        };
        // a row with no entries (a band's) is only scaled
        var want = if (terms == 0) acc else check.mulW(wa, acc);
        var scale = check.mag(wa) * mag;
        if (!isZero(T, beta)) {
            want = check.addW(want, check.mulW(wb, y.w(i)));
            scale += check.mag(wb) * check.mag(y.w(i));
        }
        check.expectWideAnyOrder(T, want, y.at(i), @floatFromInt(4 * (cols + 2)), scale) catch |e| {
            check.say("y[{d}]: ", .{i});
            return e;
        };
    }
    try y.gaps();
}

/// `x := op(A)*x`, or the solve of `op(A)*z = x`.
fn trRun(comptime T: type, comptime r: Routine, v: usize, c: Case) !void {
    const al = std.testing.allocator;
    const s = shapeOf(r, v, c);
    const op: blas.Trans = @fromBackingInt(@intCast((v / 2) % 4));
    const tr = op == .t or op == .c;
    const cj = op == .r or op == .c;
    var a: Matrix(T) = try .init(al, s, 7 * c.n + v, c.on(.a));
    defer a.deinit(al);
    var x: Vector(T) = try .init(al, s.n, c.inc[0], s.n + 3, c.on(.x));
    defer x.deinit(al);
    const pa = poison.allocator();
    const uplo: blas.Uplo = if (s.lower) .l else .u;
    const diag: blas.Diag = if (s.unit) .u else .n;
    const n: isize = @intCast(s.n);
    const k: isize = @intCast(s.kl + s.ku);
    switch (r) {
        .trmv => try blas.trmv(pa, check.io, T, uplo, op, diag, n, a.buf.ptr, a.ld(), x.base(), c.inc[0]),
        .tbmv => try blas.tbmv(pa, check.io, T, uplo, op, diag, n, k, a.buf.ptr, a.ld(), x.base(), c.inc[0]),
        .tpmv => try blas.tpmv(pa, check.io, T, uplo, op, diag, n, a.buf.ptr, x.base(), c.inc[0]),
        .trsv => try blas.trsv(pa, check.io, T, uplo, op, diag, n, a.buf.ptr, a.ld(), x.base(), c.inc[0]),
        .tbsv => try blas.tbsv(pa, check.io, T, uplo, op, diag, n, k, a.buf.ptr, a.ld(), x.base(), c.inc[0]),
        .tpsv => try blas.tpsv(pa, check.io, T, uplo, op, diag, n, a.buf.ptr, x.base(), c.inc[0]),
        else => unreachable,
    }
    try a.unchanged();
    // op(A)'s (i, j)
    const M = struct {
        fn at(m: *const Matrix(T), t: bool, conj: bool, i: usize, j: usize) ?Wide(T) {
            const e = (if (t) m.at(j, i) else m.at(i, j)) orelse return null;
            return if (conj) check.conjW(e) else e;
        }
    };
    var z: [max_order]Wide(T) = undefined;
    const lower = s.lower != tr;
    for (0..s.n) |q| {
        // the solve goes from the end its triangle starts at
        const i = if (!r.solve() or lower) q else s.n - 1 - q;
        var acc = check.wide(T, val(T, 0, 0));
        var mag: f128 = 0;
        for (0..s.n) |j| {
            if (r.solve() and (j == i or (if (lower) j > i else j < i))) continue;
            const e = M.at(&a, tr, cj, i, j) orelse continue;
            const u = if (r.solve()) z[j] else x.w(j);
            // a unit diagonal is not multiplied by
            acc = check.addW(acc, if (i == j and s.unit) u else check.mulW(e, u));
            mag += check.mag(e) * check.mag(u);
        }
        if (r.solve()) {
            const d = M.at(&a, tr, cj, i, i).?;
            const rhs = check.subW(x.w(i), acc);
            z[i] = if (s.unit) rhs else check.divW(rhs, d);
            mag = (mag + check.mag(x.w(i))) / check.mag(d);
            check.expectWideAnyOrder(T, z[i], x.at(i), @floatFromInt(16 * (s.n + 2)), mag) catch |e| {
                check.say("z[{d}]: ", .{i});
                return e;
            };
        } else {
            check.expectWideAnyOrder(T, acc, x.at(i), @floatFromInt(4 * (s.n + 2)), mag) catch |e| {
                check.say("x[{d}]: ", .{i});
                return e;
            };
        }
    }
    try x.gaps();
}

/// The rank updates: `A += alpha*x*y'` and its symmetric and Hermitian kin.
fn rankRun(comptime T: type, comptime r: Routine, v: usize, alpha: T, c: Case) !void {
    const al = std.testing.allocator;
    const s = shapeOf(r, v, c);
    const sym2 = r == .syr2 or r == .spr2 or r == .her2 or r == .hpr2;
    const two = switch (r) {
        .ger, .gerc, .syr2, .her2, .spr2, .hpr2 => true,
        else => false,
    };
    var a: Matrix(T) = try .init(al, s, 7 * c.m + c.n + v, c.on(.a));
    defer a.deinit(al);
    var x: Vector(T) = try .init(al, s.m, c.inc[0], s.m + 3, c.on(.x));
    defer x.deinit(al);
    var y: Vector(T) = try .init(al, s.n, c.inc[1], s.n + 5, c.on(.y));
    defer y.deinit(al);
    // a zero `y` for `x`'s special values to skip
    if ((r == .ger or r == .gerc) and c.where == .x and s.n > 1) {
        y.was[1] = val(T, 0, 0);
        y.buf[y.idx(1)] = y.was[1];
    }
    const az = isZero(T, alpha);
    if (az) {
        x.poison();
        y.poison();
    }
    const pa = poison.allocator();
    const m: isize = @intCast(s.m);
    const n: isize = @intCast(s.n);
    const ix = c.inc[0];
    const iy = c.inc[1];
    const ar = if (comptime isComplex(T)) alpha.re else alpha;
    const up: blas.Uplo = if (s.lower) .l else .u;
    switch (r) {
        .ger => try blas.ger(pa, check.io, T, m, n, alpha, x.base(), ix, y.base(), iy, a.buf.ptr, a.ld()),
        .gerc => if (comptime isComplex(T)) try blas.gerc(pa, check.io, T, m, n, alpha, x.base(), ix, y.base(), iy, a.buf.ptr, a.ld()),
        .syr => try blas.syr(pa, check.io, T, up, n, alpha, x.base(), ix, a.buf.ptr, a.ld()),
        .spr => try blas.spr(pa, check.io, T, up, n, alpha, x.base(), ix, a.buf.ptr),
        .her => if (comptime isComplex(T)) try blas.her(pa, check.io, T, @fromBackingInt(@intCast(v)), n, ar, x.base(), ix, a.buf.ptr, a.ld()),
        .hpr => if (comptime isComplex(T)) try blas.hpr(pa, check.io, T, @fromBackingInt(@intCast(v)), n, ar, x.base(), ix, a.buf.ptr),
        .syr2 => try blas.syr2(pa, check.io, T, up, n, alpha, x.base(), ix, y.base(), iy, a.buf.ptr, a.ld()),
        .spr2 => try blas.spr2(pa, check.io, T, up, n, alpha, x.base(), ix, y.base(), iy, a.buf.ptr),
        .her2 => if (comptime isComplex(T)) try blas.her2(pa, check.io, T, @fromBackingInt(@intCast(v)), n, alpha, x.base(), ix, y.base(), iy, a.buf.ptr, a.ld()),
        .hpr2 => if (comptime isComplex(T)) try blas.hpr2(pa, check.io, T, @fromBackingInt(@intCast(v)), n, alpha, x.base(), ix, y.base(), iy, a.buf.ptr),
        else => unreachable,
    }
    try x.unchanged();
    try y.unchanged();
    if (az) return a.unchanged();
    const herm = s.kind == .hermitian or s.kind == .hermitian_rev;
    // `her` and `hpr` take the real part
    const wa = check.wide(T, if (r == .her or r == .hpr) val(T, ar, 0) else alpha);
    // everything but the stored entries keeps its NaN
    for (a.buf, 0..) |e, k| {
        var hit = false;
        for (0..s.n) |j| for (0..s.m) |i| {
            hit = hit or s.stored(i, j) == k;
        };
        if (!hit and !check.poisoned(T, e)) {
            check.say("storage {d} written: ", .{k});
            return error.TestUnexpectedResult;
        }
    }
    for (0..s.n) |j| for (0..s.m) |i| {
        const k = s.stored(i, j) orelse continue;
        const xi = x.w(i);
        const yj = if (two) y.w(j) else x.w(j);
        // a column whose factors are zero is passed over, as the reference
        // does
        const skip = if (r == .ger or r == .gerc) zeroW(y.w(j)) else zeroW(x.w(j)) and (!sym2 or zeroW(y.w(j)));
        var u = if (skip) check.wide(T, val(T, 0, 0)) else switch (r) {
            .ger, .syr, .spr => check.mulW(wa, check.mulW(xi, yj)),
            .gerc => check.mulW(wa, check.mulW(xi, check.conjW(yj))),
            // a real alpha, part by part
            .her, .hpr => check.scaleW(ar, check.mulW(xi, check.conjW(yj))),
            .syr2, .spr2 => check.addW(check.mulW(wa, check.mulW(xi, yj)), check.mulW(wa, check.mulW(y.w(i), x.w(j)))),
            .her2, .hpr2 => check.addW(check.mulW(wa, check.mulW(xi, check.conjW(yj))), check.mulW(check.conjW(wa), check.mulW(y.w(i), check.conjW(x.w(j))))),
            else => unreachable,
        };
        if (s.kind == .hermitian_rev) u = check.conjW(u);
        var old = check.wide(T, a.was[k]);
        if (comptime isComplex(T)) if (herm and i == j) {
            old.im = 0;
            u.im = 0;
        };
        const scale = check.mag(old) + check.mag(wa) * (check.mag(xi) * check.mag(yj) + if (sym2) check.mag(y.w(i)) * check.mag(x.w(j)) else 0);
        check.expectWideAnyOrder(T, check.addW(old, u), a.buf[k], 8, scale) catch |e| {
            check.say("A({d},{d}): ", .{ i, j });
            return e;
        };
    };
}

/// The first failure of each routine, so one run names every one.
var failed: ?anyerror = null;

/// Runs `r` over its shapes, variants, operands and patterns, then its
/// special scalars, for `T`; names a failing case and goes on with the next
/// routine.
fn sweep(comptime T: type, comptime r: Routine) void {
    const complex = comptime isComplex(T);
    if (comptime r.complexOnly() and !complex) return;
    const gen = val(T, 1.25, -0.5);
    const beta = val(T, -0.75, 0.25);
    const nan = std.math.nan(Float(T));
    const zero = val(T, 0, 0);
    const wheres: []const Where = switch (comptime r.family()) {
        .mv => &.{ .a, .x, .y },
        .tr => &.{ .a, .x },
        .rank => if (r == .syr or r == .spr or r == .her or r == .hpr) &.{ .a, .x } else &.{ .a, .x, .y },
    };
    var turn: usize = 0;
    for (shapes()) |sh| for (0..r.variants(complex)) |v| {
        if (thinned(sh, v)) continue;
        const base = baseCase(r, sh);
        for (wheres) |where| {
            const s = shapeOf(r, v, base);
            const t = (r == .gemv or r == .gbmv) and trans_flags[v][0];
            const len = switch (where) {
                .a => s.entries(),
                .x => if (r.family() == .rank or t) s.m else s.n,
                .y => if (r.family() == .rank or t) s.n else s.m,
                .none => unreachable,
            };
            if (len == 0) continue;
            var buf: [13]Pattern = undefined;
            for (check.patternsFor(len, &buf)) |p| {
                // huge entries overflow a solve's partial sums, which the
                // reference's do not
                if (r.solve() and where == .a and p == .huge) continue;
                // nor do `x` times itself, where `x` holds huge values
                if ((r == .syr or r == .spr or r == .her or r == .hpr) and where == .x and (p == .huge or p == .mixed)) continue;
                var c = base;
                c.p = p;
                c.where = where;
                turn += 1;
                c.inc = if (turn % 2 == 0) .{ 1, 1 } else .{ -2, 3 };
                run(T, r, v, gen, beta, c) catch |e| {
                    report(e, T, r, v, c);
                    // Accelerate's run reports every failing case
                    if (!accel_run) return;
                };
            }
        }
        // the scalars that are not propagation, and NaN ones that are
        var c = base;
        c.inc = .{ -1, 2 };
        if (r.family() == .mv) for ([_][2]T{ .{ gen, zero }, .{ zero, beta }, .{ val(T, nan, 0), beta }, .{ gen, val(T, nan, 0) } }) |ab| {
            run(T, r, v, ab[0], ab[1], c) catch |e| {
                check.say("alpha {any}, beta {any}: ", .{ ab[0], ab[1] });
                return report(e, T, r, v, c);
            };
        };
        if (r.family() == .rank) for ([_]T{ zero, val(T, nan, 0) }) |a| {
            run(T, r, v, a, beta, c) catch |e| {
                check.say("alpha {any}: ", .{a});
                return report(e, T, r, v, c);
            };
        };
    };
}

/// m, n, kl, ku; the rigs take one, past the unrolls.
const all_shapes = [_][4]usize{ .{ 1, 1, 0, 0 }, .{ 7, 4, 1, 2 }, .{ 4, 17, 3, 0 }, .{ 33, 9, 2, 2 } };

/// Accelerate's: every shape, and one past its small-problem paths.
const accel_shapes = all_shapes ++ [_][4]usize{.{ 130, 67, 4, 6 }};
/// The largest order of any shape.
const max_order = 130;

/// Accelerate's largest shape takes every third variant: its `f128`
/// reference is slow.
fn thinned(sh: [4]usize, v: usize) bool {
    return sh[0] == max_order and v % 3 != 0;
}
var accel_run = false;

fn shapes() []const [4]usize {
    if (accel_run) return &accel_shapes;
    return if (check.rig()) all_shapes[2..3] else &all_shapes;
}

/// `r`'s case of shape `sh`, no pattern yet.
fn baseCase(comptime r: Routine, sh: [4]usize) Case {
    const general = r == .gemv or r == .gbmv or r == .ger or r == .gerc;
    // square routines take `n` from the shape's larger side
    const sq = @max(sh[0], sh[1]);
    return .{ .m = if (general) sh[0] else sq, .n = if (general) sh[1] else sq, .kl = sh[2], .ku = sh[3], .inc = .{ 1, 1 }, .p = .mixed, .where = .none };
}

/// Denormal results: denormal the operands that reach them (`x` tiny where
/// it meets itself, so its products are), the rest and the scalars near 1,
/// so a flushed input or result shows.
fn denormals(comptime T: type, comptime r: Routine) void {
    const complex = comptime isComplex(T);
    if (comptime r.complexOnly() and !complex) return;
    const self = r == .syr or r == .spr or r == .her or r == .hpr;
    for (shapes()) |sh| for (0..r.variants(complex)) |v| for ([_][2]isize{ .{ 1, 1 }, .{ -2, 3 } }) |inc| {
        if (thinned(sh, v)) continue;
        var c = baseCase(r, sh);
        c.inc = inc;
        c.p = .denormal;
        switch (comptime r.family()) {
            .mv => {
                c.where = .a;
                c.also = .{ .where = .y, .p = .denormal };
            },
            .tr => c.where = .x,
            .rank => {
                c.where = .a;
                c.also = .{ .where = .x, .p = if (self) .tiny else .denormal };
            },
        }
        run(T, r, v, val(T, 1.25, -0.5), val(T, -0.75, 0.25), c) catch |e| return report(e, T, r, v, c);
    };
}

fn run(comptime T: type, comptime r: Routine, v: usize, alpha: T, beta: T, c: Case) !void {
    return switch (comptime r.family()) {
        .mv => mvRun(T, r, v, alpha, beta, c),
        .tr => trRun(T, r, v, c),
        .rank => rankRun(T, r, v, alpha, c),
    };
}

fn report(e: anyerror, comptime T: type, comptime r: Routine, v: usize, c: Case) void {
    if (c.where != .none) c.p.describe();
    check.say("in {s}, {s} {s} variant {d} m={d} n={d} k={d},{d} inc={d},{d}\n", .{ @tagName(c.where), @tagName(r), @typeName(T), v, c.m, c.n, c.kl, c.ku, c.inc[0], c.inc[1] });
    failed = e;
}

/// Real and complex apart: one test of both is too big an image for
/// sparc64's debug rig.
const real_types = .{ f32, f64 };
const complex_types = .{ C(f32), C(f64) };

/// `vt`'s kernels, with its own `dtb_entries` and with 3, which takes the
/// triangular drivers through `gemv`; static, as a vtable is large for a
/// rig's stack.
var small_dtb: Vtable = undefined;

fn CheckAll(comptime types: anytype) type {
    return struct {
        fn sweepAll(vt: *const Vtable) !void {
            failed = null;
            {
                const saved = blas.testUse(vt);
                defer blas.testRestore(saved);
                inline for (types) |T| inline for (comptime std.enums.values(Routine)) |r| sweep(T, r);
            }
            small_dtb = vt.*;
            small_dtb.dtb_entries = 3;
            const saved = blas.testUse(&small_dtb);
            defer blas.testRestore(saved);
            inline for (types) |T| inline for (.{ Routine.trmv, Routine.trsv }) |r| sweep(T, r);
            if (failed) |e| return e;
        }

        fn denormalsAll(vt: *const Vtable) !void {
            failed = null;
            const saved = blas.testUse(vt);
            defer blas.testRestore(saved);
            inline for (types) |T| inline for (comptime std.enums.values(Routine)) |r| denormals(T, r);
            if (failed) |e| return e;
        }

        /// Through Apple Accelerate, where the routines go to it.
        fn accelerate() !void {
            if (comptime !accel.covers(f64)) return error.SkipZigTest;
            failed = null;
            const saved = blas.testUse(blas.vtable());
            defer blas.testRestore(saved);
            if (!accel.route(true)) return error.SkipZigTest;
            accel_run = true;
            defer accel_run = false;
            inline for (types) |T| inline for (comptime std.enums.values(Routine)) |r| {
                sweep(T, r);
                denormals(T, r);
            };
            if (failed) |e| return e;
        }

        fn portableKernels() !void {
            var vt: Vtable = .{};
            portable.install(&vt);
            try sweepAll(&vt);
            try denormalsAll(&vt);
        }
    };
}

test "portable level-2 routines, real: special values, and NaN wherever they must not read or write" {
    try CheckAll(real_types).portableKernels();
}

test "portable level-2 routines, complex: special values, and NaN wherever they must not read or write" {
    try CheckAll(complex_types).portableKernels();
}

test "every linked core's level-2 routines, real: special values, and NaN wherever they must not read or write" {
    try check.eachCore(CheckAll(real_types).sweepAll);
}

test "every linked core's level-2 routines, complex: special values, and NaN wherever they must not read or write" {
    try check.eachCore(CheckAll(complex_types).sweepAll);
}

test "the level-2 routines through Accelerate, real: special values, denormals, and NaN wherever they must not read or write" {
    try CheckAll(real_types).accelerate();
}

test "the level-2 routines through Accelerate, complex: special values, denormals, and NaN wherever they must not read or write" {
    try CheckAll(complex_types).accelerate();
}

test "every linked core's level-2 routines, real: denormal results" {
    try check.eachCore(CheckAll(real_types).denormalsAll);
}

test "every linked core's level-2 routines, complex: denormal results" {
    try check.eachCore(CheckAll(complex_types).denormalsAll);
}
