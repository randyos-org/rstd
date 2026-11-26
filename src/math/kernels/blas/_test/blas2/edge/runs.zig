//! One case of a level-2 routine against the `f128` reference: products, solves and rank updates.

const std = @import("std");
const blas = @import("../../../__root__.zig");
const l2 = @import("../../../blas2/core.zig");
const check = @import("../../check.zig");
const bitwise = @import("../../bitwise.zig");
const pattern = @import("../../pattern.zig");
const reference = @import("../../reference.zig");
const ed = @import("../edge.zig");
const mx = @import("matrix.zig");
const cs = @import("case.zig");

const isComplex = check.isComplex;
const Wide = reference.Wide;
const Vector = check.Vector;
const max_order = ed.max_order;
const Matrix = mx.Matrix;
const Case = cs.Case;
const Routine = cs.Routine;
const isZero = cs.isZero;
const shapeOf = cs.shapeOf;
const trans_flags = cs.trans_flags;
const val = cs.val;
const zeroW = cs.zeroW;

/// `y := alpha*op(A)*x + beta*y` for the matrix-vector products.
pub fn mvRun(comptime T: type, comptime r: Routine, v: usize, alpha: T, beta: T, c: Case) !void {
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
    const pa = cs.poison.allocator();
    const m: isize = @intCast(s.m);
    const n: isize = @intCast(s.n);
    const ix = c.inc[0];
    const iy = c.inc[1];
    switch (r) {
        .gemv => try l2.gemv(pa, check.io, T, @fromBackingInt(@intCast(v)), m, n, alpha, a.buf.ptr, a.ld(), x.base(), ix, beta, y.base(), iy),
        .gbmv => try l2.gbmv(pa, check.io, T, @fromBackingInt(@intCast(v)), m, n, @intCast(s.kl), @intCast(s.ku), alpha, a.buf.ptr, a.ld(), x.base(), ix, beta, y.base(), iy),
        .symv => try l2.symv(pa, check.io, T, if (s.lower) .l else .u, n, alpha, a.buf.ptr, a.ld(), x.base(), ix, beta, y.base(), iy),
        .hemv => if (comptime isComplex(T)) try l2.hemv(pa, check.io, T, @fromBackingInt(@intCast(v)), n, alpha, a.buf.ptr, a.ld(), x.base(), ix, beta, y.base(), iy),
        .sbmv => try l2.sbmv(pa, check.io, T, if (s.lower) .l else .u, n, @intCast(s.kl + s.ku), alpha, a.buf.ptr, a.ld(), x.base(), ix, beta, y.base(), iy),
        .hbmv => if (comptime isComplex(T)) try l2.hbmv(pa, check.io, T, @fromBackingInt(@intCast(v)), n, @intCast(s.kl + s.ku), alpha, a.buf.ptr, a.ld(), x.base(), ix, beta, y.base(), iy),
        .spmv => try l2.spmv(pa, check.io, T, if (s.lower) .l else .u, n, alpha, a.buf.ptr, x.base(), ix, beta, y.base(), iy),
        .hpmv => if (comptime isComplex(T)) try l2.hpmv(pa, check.io, T, @fromBackingInt(@intCast(v)), n, alpha, a.buf.ptr, x.base(), ix, beta, y.base(), iy),
        else => unreachable,
    }
    try a.unchanged();
    try x.unchanged();
    const wa = reference.wide(T, alpha);
    const wb = reference.wide(T, beta);
    for (0..rows) |i| {
        var acc = reference.wide(T, val(T, 0, 0));
        var mag: f128 = 0;
        var terms: usize = 0;
        if (!isZero(T, alpha)) for (0..cols) |j| {
            var e = (if (t[0]) a.at(j, i) else a.at(i, j)) orelse continue;
            terms += 1;
            if (t[1]) e = reference.conjW(e);
            const xv = if (t[2]) reference.conjW(x.w(j)) else x.w(j);
            acc = reference.addW(acc, reference.mulW(e, xv));
            mag += reference.mag(e) * reference.mag(xv);
        };
        // a row with no entries (a band's) is only scaled
        var want = if (terms == 0) acc else reference.mulW(wa, acc);
        var scale = reference.mag(wa) * mag;
        if (!isZero(T, beta)) {
            want = reference.addW(want, reference.mulW(wb, y.w(i)));
            scale += reference.mag(wb) * reference.mag(y.w(i));
        }
        pattern.expectWideAnyOrder(T, want, y.at(i), @floatFromInt(4 * (cols + 2)), scale) catch |e| {
            check.say("y[{d}]: ", .{i});
            return e;
        };
    }
    try y.gaps();
}

/// `x := op(A)*x`, or the solve of `op(A)*z = x`.
pub fn trRun(comptime T: type, comptime r: Routine, v: usize, c: Case) !void {
    const al = std.testing.allocator;
    const s = shapeOf(r, v, c);
    const op: blas.Trans = @fromBackingInt(@intCast((v / 2) % 4));
    const tr = op == .t or op == .c;
    const cj = op == .r or op == .c;
    var a: Matrix(T) = try .init(al, s, 7 * c.n + v, c.on(.a));
    defer a.deinit(al);
    var x: Vector(T) = try .init(al, s.n, c.inc[0], s.n + 3, c.on(.x));
    defer x.deinit(al);
    const pa = cs.poison.allocator();
    const uplo: blas.Uplo = if (s.lower) .l else .u;
    const diag: blas.Diag = if (s.unit) .u else .n;
    const n: isize = @intCast(s.n);
    const k: isize = @intCast(s.kl + s.ku);
    switch (r) {
        .trmv => try l2.trmv(pa, check.io, T, uplo, op, diag, n, a.buf.ptr, a.ld(), x.base(), c.inc[0]),
        .tbmv => try l2.tbmv(pa, check.io, T, uplo, op, diag, n, k, a.buf.ptr, a.ld(), x.base(), c.inc[0]),
        .tpmv => try l2.tpmv(pa, check.io, T, uplo, op, diag, n, a.buf.ptr, x.base(), c.inc[0]),
        .trsv => try l2.trsv(pa, check.io, T, uplo, op, diag, n, a.buf.ptr, a.ld(), x.base(), c.inc[0]),
        .tbsv => try l2.tbsv(pa, check.io, T, uplo, op, diag, n, k, a.buf.ptr, a.ld(), x.base(), c.inc[0]),
        .tpsv => try l2.tpsv(pa, check.io, T, uplo, op, diag, n, a.buf.ptr, x.base(), c.inc[0]),
        else => unreachable,
    }
    try a.unchanged();
    // op(A)'s (i, j)
    const M = struct {
        pub fn at(m: *const Matrix(T), t: bool, conj: bool, i: usize, j: usize) ?Wide(T) {
            const e = (if (t) m.at(j, i) else m.at(i, j)) orelse return null;
            return if (conj) reference.conjW(e) else e;
        }
    };
    var z: [max_order]Wide(T) = undefined;
    const lower = s.lower != tr;
    for (0..s.n) |q| {
        // the solve goes from the end its triangle starts at
        const i = if (!r.solve() or lower) q else s.n - 1 - q;
        var acc = reference.wide(T, val(T, 0, 0));
        var mag: f128 = 0;
        for (0..s.n) |j| {
            if (r.solve() and (j == i or (if (lower) j > i else j < i))) continue;
            const e = M.at(&a, tr, cj, i, j) orelse continue;
            const u = if (r.solve()) z[j] else x.w(j);
            // a unit diagonal is not multiplied by
            acc = reference.addW(acc, if (i == j and s.unit) u else reference.mulW(e, u));
            mag += reference.mag(e) * reference.mag(u);
        }
        if (r.solve()) {
            const d = M.at(&a, tr, cj, i, i).?;
            const rhs = reference.subW(x.w(i), acc);
            z[i] = if (s.unit) rhs else reference.divW(rhs, d);
            mag = (mag + reference.mag(x.w(i))) / reference.mag(d);
            pattern.expectWideAnyOrder(T, z[i], x.at(i), reference.solveUlps(s.n), mag) catch |e| {
                check.say("z[{d}]: ", .{i});
                return e;
            };
        } else {
            pattern.expectWideAnyOrder(T, acc, x.at(i), @floatFromInt(4 * (s.n + 2)), mag) catch |e| {
                check.say("x[{d}]: ", .{i});
                return e;
            };
        }
    }
    try x.gaps();
}

pub fn setAt(comptime T: type, v: *Vector(T), i: usize, e: T) void {
    v.was[i] = e;
    v.buf[v.idx(i)] = e;
}

/// Zero factors beside a vector's special values, away from a lone one:
/// at `z1` a column the reference passes over (`y` for `ger`, `x`, both for
/// rank 2), as `-0`; for rank 2, at `z2` only the other vector's zero, a
/// column it does not, so the special value meets that zero. They take the
/// first indices, or the last for an upper triangle, whose columns are full.
pub fn zeroFactors(comptime T: type, comptime r: Routine, c: Case, lower: bool, x: *Vector(T), y: *Vector(T)) void {
    if (c.where != .x and c.where != .y) return;
    const general = r == .ger or r == .gerc;
    if (general and c.where == .y) return;
    const n = if (general) y.n else x.n;
    const lone: ?usize = switch (c.p) {
        .one => |o| o.at,
        else => null,
    };
    var z: [2]usize = undefined;
    var found: usize = 0;
    for (0..n) |k| {
        const i = if (lower or general) k else n - 1 - k;
        if (found < 2 and i != lone) {
            z[found] = i;
            found += 1;
        }
    }
    if (found == 0) return;
    const nz = val(T, -0.0, 0);
    switch (r) {
        .ger, .gerc => setAt(T, y, z[0], nz),
        .syr, .spr, .her, .hpr => setAt(T, x, z[0], nz),
        .syr2, .spr2, .her2, .hpr2 => {
            setAt(T, x, z[0], nz);
            setAt(T, y, z[0], nz);
            if (found == 2) setAt(T, if (c.where == .x) y else x, z[1], val(T, 0, 0));
        },
        else => unreachable,
    }
}

/// The rank updates: `A += alpha*x*y'` and its symmetric and Hermitian kin.
pub fn rankRun(comptime T: type, comptime r: Routine, v: usize, alpha: T, c: Case) !void {
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
    zeroFactors(T, r, c, s.lower, &x, &y);
    const az = isZero(T, alpha);
    if (az) {
        x.poison();
        y.poison();
    }
    const pa = cs.poison.allocator();
    const m: isize = @intCast(s.m);
    const n: isize = @intCast(s.n);
    const ix = c.inc[0];
    const iy = c.inc[1];
    const ar = if (comptime isComplex(T)) alpha.re else alpha;
    const up: blas.Uplo = if (s.lower) .l else .u;
    switch (r) {
        .ger => try l2.ger(pa, check.io, T, m, n, alpha, x.base(), ix, y.base(), iy, a.buf.ptr, a.ld()),
        .gerc => if (comptime isComplex(T)) try l2.gerc(pa, check.io, T, m, n, alpha, x.base(), ix, y.base(), iy, a.buf.ptr, a.ld()),
        .syr => try l2.syr(pa, check.io, T, up, n, alpha, x.base(), ix, a.buf.ptr, a.ld()),
        .spr => try l2.spr(pa, check.io, T, up, n, alpha, x.base(), ix, a.buf.ptr),
        .her => if (comptime isComplex(T)) try l2.her(pa, check.io, T, @fromBackingInt(@intCast(v)), n, ar, x.base(), ix, a.buf.ptr, a.ld()),
        .hpr => if (comptime isComplex(T)) try l2.hpr(pa, check.io, T, @fromBackingInt(@intCast(v)), n, ar, x.base(), ix, a.buf.ptr),
        .syr2 => try l2.syr2(pa, check.io, T, up, n, alpha, x.base(), ix, y.base(), iy, a.buf.ptr, a.ld()),
        .spr2 => try l2.spr2(pa, check.io, T, up, n, alpha, x.base(), ix, y.base(), iy, a.buf.ptr),
        .her2 => if (comptime isComplex(T)) try l2.her2(pa, check.io, T, @fromBackingInt(@intCast(v)), n, alpha, x.base(), ix, y.base(), iy, a.buf.ptr, a.ld()),
        .hpr2 => if (comptime isComplex(T)) try l2.hpr2(pa, check.io, T, @fromBackingInt(@intCast(v)), n, alpha, x.base(), ix, y.base(), iy, a.buf.ptr),
        else => unreachable,
    }
    try x.unchanged();
    try y.unchanged();
    if (az) return a.unchanged();
    const herm = s.kind == .hermitian;
    // `her` and `hpr` take the real part
    const wa = reference.wide(T, if (r == .her or r == .hpr) val(T, ar, 0) else alpha);
    // everything but the stored entries keeps its NaN
    const hit = try al.alloc(bool, a.buf.len);
    defer al.free(hit);
    @memset(hit, false);
    for (0..s.n) |j| for (0..s.m) |i| if (s.stored(i, j)) |k| {
        hit[k] = true;
    };
    for (a.buf, hit, 0..) |e, h, k| {
        if (!h and !bitwise.poisoned(T, e)) {
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
        var u = if (skip) reference.wide(T, val(T, 0, 0)) else switch (r) {
            .ger, .syr, .spr => reference.mulW(wa, reference.mulW(xi, yj)),
            .gerc => reference.mulW(wa, reference.mulW(xi, reference.conjW(yj))),
            // a real alpha, part by part
            .her, .hpr => reference.scaleW(ar, reference.mulW(xi, reference.conjW(yj))),
            .syr2, .spr2 => reference.addW(reference.mulW(wa, reference.mulW(xi, yj)), reference.mulW(wa, reference.mulW(y.w(i), x.w(j)))),
            .her2, .hpr2 => reference.addW(reference.mulW(wa, reference.mulW(xi, reference.conjW(yj))), reference.mulW(reference.conjW(wa), reference.mulW(y.w(i), reference.conjW(x.w(j))))),
            else => unreachable,
        };
        var old = reference.wide(T, a.was[k]);
        if (comptime isComplex(T)) if (herm and i == j) {
            old.im = 0;
            u.im = 0;
        };
        const scale = reference.mag(old) + reference.mag(wa) * (reference.mag(xi) * reference.mag(yj) + if (sym2) reference.mag(y.w(i)) * reference.mag(x.w(j)) else 0);
        pattern.expectWideAnyOrder(T, reference.addW(old, u), a.buf[k], reference.update_ulps, scale) catch |e| {
            check.say("A({d},{d}): ", .{ i, j });
            return e;
        };
    };
}
