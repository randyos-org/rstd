//! HESS (`gehrd`, `orghr`) and SCHUR (`gees`, complex always) over LAPACK
//! (M2), in place of EISPACK's CORTH and COMQR3 (`MATFN4`, `MATFN5`).

const std = @import("std");
const zigstd = @import("../../zigstd.zig");
const scalar = @import("scalar.zig");
const CxMat = @import("CxMat.zig");
const lapack = @import("../lapack/__root__.zig");
const ops = @import("../lapack/ops.zig");
const eig = @import("eig.zig");

const Cx = scalar.Cx;
const Allocator = std.mem.Allocator;
const Io = std.Io;
const ix = scalar.ix;
const EigError = eig.EigError;
const converged = eig.converged;

/// `<P, H>`: `A = P*H*P'`, `P` unitary, `H` Hessenberg (tridiagonal for a
/// Hermitian `A`). Caller frees.
pub const HessResult = struct {
    p_r: []f64,
    p_i: []f64,
    h_r: []f64,
    h_i: []f64,

    pub fn deinit(self: HessResult, allocator: Allocator) void {
        allocator.free(self.p_r);
        allocator.free(self.p_i);
        allocator.free(self.h_r);
        allocator.free(self.h_i);
    }
};

pub const HessError = Allocator.Error || Io.Cancelable;

/// HESS: see `HessResult`. `a` is not mutated (a packed copy is reduced
/// instead).
pub fn hess(allocator: Allocator, io: Io, a: CxMat) HessError!HessResult {
    const n = a.n;
    const h = try a.clone(allocator);
    errdefer h.deinit(allocator);
    const p = try a.clone(allocator);
    errdefer p.deinit(allocator);
    const ni: isize = @intCast(n);
    const ld: isize = @intCast(@max(1, n));
    switch (a.isReal()) {
        inline else => |real| {
            const T = if (real) f64 else Cx;
            const x = try a.toPacked(T, allocator);
            defer allocator.free(x);
            const tau = try allocator.alloc(T, n -| 1);
            defer allocator.free(tau);
            // the arguments come from `a`
            lapack.gehrd(allocator, io, T, ni, 1, ni, x.ptr, ld, tau.ptr) catch |e| return ops.narrow(HessError, e);
            h.fromPacked(T, x);
            lapack.orghr(allocator, io, T, ni, 1, ni, x.ptr, ld, tau.ptr) catch |e| return ops.narrow(HessError, e);
            p.fromPacked(T, x);
        },
    }
    // zero below the first subdiagonal, where `gehrd` keeps `Q`
    h.triu(-1);
    return .{ .p_r = p.r, .p_i = p.i, .h_r = h.r, .h_i = h.i };
}

/// `<U, T>`: `A`'s unitary similarity `U` and (upper triangular) Schur
/// form `T`, with the eigenvalues of `A` on `T`'s diagonal, so
/// `A = U*T*U'`. Caller owns the returned slices (free with `allocator`).
pub const SchurResult = struct {
    u_r: []f64,
    u_i: []f64,
    t_r: []f64,
    t_i: []f64,

    pub fn deinit(self: SchurResult, allocator: Allocator) void {
        allocator.free(self.u_r);
        allocator.free(self.u_i);
        allocator.free(self.t_r);
        allocator.free(self.t_i);
    }
};

pub const SchurError = EigError;

/// SCHUR: see `SchurResult`, by `zgees`. `a` is not mutated (a packed copy
/// is reduced instead). `error.NoConvergence` as `eig`'s.
pub fn schur(allocator: Allocator, io: Io, a: CxMat) SchurError!SchurResult {
    const n = a.n;
    const t = try a.clone(allocator);
    errdefer t.deinit(allocator);
    const u = try a.clone(allocator);
    errdefer u.deinit(allocator);
    const ni: isize = @intCast(n);
    const ld: isize = @intCast(@max(1, n));
    const x = try a.toPacked(Cx, allocator);
    defer allocator.free(x);
    const vs = try allocator.alloc(Cx, n * n);
    defer allocator.free(vs);
    const w = try allocator.alloc(Cx, n);
    defer allocator.free(w);
    var sdim: isize = 0;
    try converged(lapack.gees(allocator, io, Cx, .v, null, ni, x.ptr, ld, &sdim, w.ptr, vs.ptr, ld, null));
    t.fromPacked(Cx, x);
    u.fromPacked(Cx, vs);
    // `zgees` leaves `T` triangular; zero what it does not write
    t.triu(0);
    return .{ .u_r = u.r, .u_i = u.i, .t_r = t.r, .t_i = t.i };
}

// ---- tests ------------------------------------------------------------------

const matmul = @import("matmul.zig").matmul;

fn matSub(alloc: Allocator, n: usize, ar: []const f64, ai: []const f64, br: []const f64, bi: []const f64) ![2][]f64 {
    const dr = try alloc.alloc(f64, n * n);
    const di = try alloc.alloc(f64, n * n);
    for (dr, di, ar, ai, br, bi) |*r, *i, a_r, a_i, b_r, b_i| {
        r.* = a_r - b_r;
        i.* = a_i - b_i;
    }
    return .{ dr, di };
}

test "hess reduces to zero below the first subdiagonal, and P*H*P' reconstructs A" {
    const n = 3;
    var ar = [_]f64{ 1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 10.0 };
    var ai: [n * n]f64 = @splat(0.0);

    var result = try hess(zigstd.testAllocator, std.testing.io, CxMat.initPacked(&ar, &ai, n, n));
    defer result.deinit(zigstd.testAllocator);

    // Zero below the first subdiagonal: only (2,1) and (3,2) survive below
    // the diagonal; (3,1) must be exactly zero.
    try zigstd.testExpectApproxEqAbs(0.0, result.h_r[ix(n, 3, 1)], 1e-9);
    try zigstd.testExpectApproxEqAbs(0.0, result.h_i[ix(n, 3, 1)], 1e-9);

    // A = P*H*P': check by reconstructing P*H*P' and comparing to A.
    var ph_r: [n * n]f64 = undefined;
    var ph_i: [n * n]f64 = undefined;
    const p_mat = CxMat.initPacked(result.p_r, result.p_i, n, n);
    const h_mat = CxMat.initPacked(result.h_r, result.h_i, n, n);
    try matmul(p_mat, h_mat, CxMat.initPacked(&ph_r, &ph_i, n, n), zigstd.testAllocator, std.testing.io);

    // P' = conjugate transpose of P.
    var pt_r: [n * n]f64 = undefined;
    var pt_i: [n * n]f64 = undefined;
    for (0..n) |row| {
        for (0..n) |col| {
            pt_r[ix(n, row + 1, col + 1)] = result.p_r[ix(n, col + 1, row + 1)];
            pt_i[ix(n, row + 1, col + 1)] = -result.p_i[ix(n, col + 1, row + 1)];
        }
    }
    var recon_r: [n * n]f64 = undefined;
    var recon_i: [n * n]f64 = undefined;
    try matmul(CxMat.initPacked(&ph_r, &ph_i, n, n), CxMat.initPacked(&pt_r, &pt_i, n, n), CxMat.initPacked(&recon_r, &recon_i, n, n), zigstd.testAllocator, std.testing.io);

    const diff = try matSub(zigstd.testAllocator, n, &recon_r, &recon_i, &ar, &ai);
    defer zigstd.testAllocator.free(diff[0]);
    defer zigstd.testAllocator.free(diff[1]);
    for (diff[0], diff[1]) |dr, di| {
        try zigstd.testExpectApproxEqAbs(0.0, dr, 1e-13);
        try zigstd.testExpectApproxEqAbs(0.0, di, 1e-13);
    }
}

test "schur produces an upper triangular T with A's eigenvalues on the diagonal" {
    const n = 3;
    var ar = [_]f64{ 1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 10.0 };
    var ai: [n * n]f64 = @splat(0.0);

    var result = try schur(zigstd.testAllocator, std.testing.io, CxMat.initPacked(&ar, &ai, n, n));
    defer result.deinit(zigstd.testAllocator);

    // Upper triangular, exactly.
    for (0..n) |col| {
        for (col + 1..n) |row| {
            try zigstd.testExpectEqual(@as(f64, 0), result.t_r[ix(n, row + 1, col + 1)]);
            try zigstd.testExpectEqual(@as(f64, 0), result.t_i[ix(n, row + 1, col + 1)]);
        }
    }

    // Trace of T (== sum of A's eigenvalues) matches trace of A.
    var trace: f64 = 0.0;
    for (0..n) |k| trace += result.t_r[ix(n, k + 1, k + 1)];
    try zigstd.testExpectApproxEqAbs(16.0, trace, 1e-13);
}

test "hess and schur of empty and 1x1 matrices: empty factors; P = U = 1, H = T = A" {
    const al = zigstd.testAllocator;
    var e: [0]f64 = .{};
    const empty = CxMat.initPacked(&e, &e, 0, 0);
    const h0 = try hess(al, std.testing.io, empty);
    defer h0.deinit(al);
    try zigstd.testExpectEqual(@as(usize, 0), h0.h_r.len + h0.p_r.len);
    const s0 = try schur(al, std.testing.io, empty);
    defer s0.deinit(al);
    try zigstd.testExpectEqual(@as(usize, 0), s0.t_r.len + s0.u_r.len);

    var ar = [_]f64{3.0};
    var ai = [_]f64{-2.0};
    const a = CxMat.initPacked(&ar, &ai, 1, 1);
    const h1 = try hess(al, std.testing.io, a);
    defer h1.deinit(al);
    const s1 = try schur(al, std.testing.io, a);
    defer s1.deinit(al);
    for ([_][4]f64{ .{ h1.p_r[0], h1.p_i[0], h1.h_r[0], h1.h_i[0] }, .{ s1.u_r[0], s1.u_i[0], s1.t_r[0], s1.t_i[0] } }) |got| {
        try zigstd.testExpectEqualSlices(f64, &.{ 1.0, 0.0, 3.0, -2.0 }, &got);
    }
}
