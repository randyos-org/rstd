//! `HESS` and `SCHUR`, ported from `matlab88.f90`'s `MATFN4`/`MATFN5`
//! command bodies, over rstd's LAPACK (M2): `gehrd` and `orghr`, and
//! `gees`, in place of EISPACK's CORTH and COMQR3. As the original's doc
//! text, both apply whether or not `A` is Hermitian; `SCHUR` always takes
//! the complex routine, so `T` is triangular even for a real `A`.

const std = @import("std");
const zigstd = @import("../../zigstd.zig");
const scalar = @import("scalar.zig");
const CxMat = @import("CxMat.zig");
const lapack = @import("../lapack/__root__.zig");
const eig_mod = @import("eig.zig");
const Complex = CxMat.Complex;
const Allocator = std.mem.Allocator;
const Io = std.Io;
const ix = scalar.ix;

/// `<P, H>`: `A`'s unitary similarity `P` and Hessenberg form `H`, so
/// `A = P*H*P'` (`H` zero below the first subdiagonal - tridiagonal if `A`
/// is Hermitian, per the original's own doc text). Caller owns the
/// returned slices (free with `allocator`).
pub const HessResult = struct {
    p_r: []f64,
    p_i: []f64,
    h_r: []f64,
    h_i: []f64,

    pub fn deinit(self: HessResult, allocator: zigstd.Allocator) void {
        allocator.free(self.p_r);
        allocator.free(self.p_i);
        allocator.free(self.h_r);
        allocator.free(self.h_i);
    }
};

pub const HessError = Allocator.Error || Io.Cancelable;

/// Zeros `m` below its first subdiagonal, where `gehrd` keeps `Q`.
fn zeroBelowSubdiagonal(m: CxMat) void {
    for (1..m.n -| 1) |j| {
        const c = m.col(j);
        @memset(c.r[j + 1 ..], 0);
        @memset(c.i[j + 1 ..], 0);
    }
}

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
            const T = if (real) f64 else Complex;
            const x = try a.toPacked(T, allocator);
            defer allocator.free(x);
            const tau = try allocator.alloc(T, n -| 1);
            defer allocator.free(tau);
            lapack.gehrd(allocator, io, T, ni, 1, ni, x.ptr, ld, tau.ptr) catch |e| switch (e) {
                error.OutOfMemory, error.Canceled => |y| return y,
                // the arguments come from `a`
                else => unreachable,
            };
            h.fromPacked(T, x);
            lapack.orghr(allocator, io, T, ni, 1, ni, x.ptr, ld, tau.ptr) catch |e| switch (e) {
                error.OutOfMemory, error.Canceled => |y| return y,
                else => unreachable,
            };
            p.fromPacked(T, x);
        },
    }
    zeroBelowSubdiagonal(h);
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

    pub fn deinit(self: SchurResult, allocator: zigstd.Allocator) void {
        allocator.free(self.u_r);
        allocator.free(self.u_i);
        allocator.free(self.t_r);
        allocator.free(self.t_i);
    }
};

pub const SchurError = eig_mod.EigError;

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
    const x = try a.toPacked(Complex, allocator);
    defer allocator.free(x);
    const vs = try allocator.alloc(Complex, n * n);
    defer allocator.free(vs);
    const w = try allocator.alloc(Complex, n);
    defer allocator.free(w);
    var sdim: isize = 0;
    try eig_mod.converged(lapack.gees(allocator, io, Complex, .v, null, ni, x.ptr, ld, &sdim, w.ptr, vs.ptr, ld, null));
    t.fromPacked(Complex, x);
    u.fromPacked(Complex, vs);
    // `zgees` leaves `T` triangular; zero what it does not write
    for (1..n) |j| {
        const c = t.col(j);
        @memset(c.r[j..], 0);
        @memset(c.i[j..], 0);
    }
    return .{ .u_r = u.r, .u_i = u.i, .t_r = t.r, .t_i = t.i };
}

fn matSub(alloc: zigstd.Allocator, n: usize, ar: []const f64, ai: []const f64, br: []const f64, bi: []const f64) ![2][]f64 {
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
    const matmul_mod = @import("matmul.zig");
    var ph_r: [n * n]f64 = undefined;
    var ph_i: [n * n]f64 = undefined;
    const p_mat = CxMat.initPacked(result.p_r, result.p_i, n, n);
    const h_mat = CxMat.initPacked(result.h_r, result.h_i, n, n);
    try matmul_mod.matmul(p_mat, h_mat, CxMat.initPacked(&ph_r, &ph_i, n, n), zigstd.testAllocator, std.testing.io);

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
    try matmul_mod.matmul(CxMat.initPacked(&ph_r, &ph_i, n, n), CxMat.initPacked(&pt_r, &pt_i, n, n), CxMat.initPacked(&recon_r, &recon_i, n, n), zigstd.testAllocator, std.testing.io);

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
