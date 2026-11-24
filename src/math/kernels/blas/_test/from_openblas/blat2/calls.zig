//! `blat2`'s calls: each routine's arguments, the call, and the sweep's bookkeeping.

const std = @import("std");
const ob = @import("../__root__.zig");
const blas = @import("../../../__root__.zig");
const l2 = @import("../../../blas2/core.zig");
const guard = @import("../../guard.zig");
const b2 = @import("../blat2.zig");
const mk = @import("make.zig");

const Float = ob.Float;
const isComplex = ob.isComplex;
const testing = std.testing;
const Params = b2.Params;
const Routine = b2.Routine;
const incmax = b2.incmax;
const letter = b2.letter;
const name = b2.name;
const nmax = b2.nmax;
const Result = mk.Result;

// ---- the calls ----------------------------------------------------------------

pub fn Work(comptime T: type) type {
    const R = ob.Ref(T);
    return struct {
        a: []T,
        aa: []T,
        as: []T,
        x: []T,
        xx: []T,
        xs: []T,
        y: []T,
        yy: []T,
        ys: []T,
        yt: []R,
        g: []Float(R),
        z: []T,

        pub fn init(al: std.mem.Allocator) !@This() {
            return .{
                .a = try al.alloc(T, nmax * nmax),
                .aa = try al.alloc(T, nmax * nmax),
                .as = try al.alloc(T, nmax * nmax),
                .x = try al.alloc(T, nmax),
                .xx = try al.alloc(T, nmax * incmax),
                .xs = try al.alloc(T, nmax * incmax),
                .y = try al.alloc(T, nmax),
                .yy = try al.alloc(T, nmax * incmax),
                .ys = try al.alloc(T, nmax * incmax),
                .yt = try al.alloc(R, nmax),
                .g = try al.alloc(Float(R), nmax),
                .z = try al.alloc(T, 2 * nmax),
            };
        }
    };
}

/// One call's arguments; `k` is `kl` too.
pub const Call = struct {
    trans: u8 = 'N',
    uplo: u8 = 'U',
    diag: u8 = 'N',
    m: usize = 0,
    n: usize = 0,
    k: isize = 0,
    ku: isize = 0,
    lda: usize = 0,
    incx: isize = 1,
    incy: isize = 1,
    alpha: [2]f64 = .{ 0, 0 },
    beta: [2]f64 = .{ 0, 0 },

    pub fn describe(c: Call, comptime T: type, r: Routine, nc: usize) void {
        ob.say("{c}{s} call {d}: trans={c} uplo={c} diag={c} m={d} n={d} k={d} ku={d} lda={d} incx={d} incy={d} alpha=({d},{d}) beta=({d},{d})\n", .{
            letter(T), name(T, r), nc, c.trans, c.uplo, c.diag, c.m, c.n, c.k, c.ku, c.lda, c.incx, c.incy, c.alpha[0], c.alpha[1], c.beta[0], c.beta[1],
        });
    }
};

pub const transOf = ob.transOf;
pub const uploOf = ob.uploOf;

/// Calls routine `r` on the storage in `w`.
pub fn run(comptime T: type, comptime r: Routine, c: Call, w: anytype) !void {
    var poison: guard.Poison = .{ .child = testing.allocator };
    const al = poison.allocator();
    const cx = comptime isComplex(T);
    const alpha = ob.val(T, c.alpha[0], c.alpha[1]);
    const beta = ob.val(T, c.beta[0], c.beta[1]);
    const m: isize = @intCast(c.m);
    const n: isize = @intCast(c.n);
    const lda: isize = @intCast(c.lda);
    const a = w.aa.ptr;
    const x = w.xx.ptr;
    const y = w.yy.ptr;
    switch (r) {
        .gemv => try l2.gemv(al, ob.io, T, transOf(c.trans), m, n, alpha, a, lda, x, c.incx, beta, y, c.incy),
        .gbmv => try l2.gbmv(al, ob.io, T, transOf(c.trans), m, n, c.k, c.ku, alpha, a, lda, x, c.incx, beta, y, c.incy),
        .symv => if (cx)
            try l2.hemv(al, ob.io, T, uploOf(c.uplo), n, alpha, a, lda, x, c.incx, beta, y, c.incy)
        else
            try l2.symv(al, ob.io, T, uploOf(c.uplo), n, alpha, a, lda, x, c.incx, beta, y, c.incy),
        .sbmv => if (cx)
            try l2.hbmv(al, ob.io, T, uploOf(c.uplo), n, c.k, alpha, a, lda, x, c.incx, beta, y, c.incy)
        else
            try l2.sbmv(al, ob.io, T, uploOf(c.uplo), n, c.k, alpha, a, lda, x, c.incx, beta, y, c.incy),
        .spmv => if (cx)
            try l2.hpmv(al, ob.io, T, uploOf(c.uplo), n, alpha, a, x, c.incx, beta, y, c.incy)
        else
            try l2.spmv(al, ob.io, T, uploOf(c.uplo), n, alpha, a, x, c.incx, beta, y, c.incy),
        .trmv, .tbmv, .tpmv, .trsv, .tbsv, .tpsv => {
            const u = ob.letter(blas.Uplo, c.uplo).?;
            const t = ob.letter(blas.Trans, c.trans).?;
            const d = ob.letter(blas.Diag, c.diag).?;
            switch (r) {
                .trmv => try l2.trmv(al, ob.io, T, u, t, d, n, a, lda, x, c.incx),
                .tbmv => try l2.tbmv(al, ob.io, T, u, t, d, n, c.k, a, lda, x, c.incx),
                .tpmv => try l2.tpmv(al, ob.io, T, u, t, d, n, a, x, c.incx),
                .trsv => try l2.trsv(al, ob.io, T, u, t, d, n, a, lda, x, c.incx),
                .tbsv => try l2.tbsv(al, ob.io, T, u, t, d, n, c.k, a, lda, x, c.incx),
                .tpsv => try l2.tpsv(al, ob.io, T, u, t, d, n, a, x, c.incx),
                else => unreachable,
            }
        },
        .ger => try l2.ger(al, ob.io, T, m, n, alpha, x, c.incx, y, c.incy, a, lda),
        .gerc => if (cx) try l2.gerc(al, ob.io, T, m, n, alpha, x, c.incx, y, c.incy, a, lda) else unreachable,
        .syr => if (cx)
            try l2.her(al, ob.io, T, uploOf(c.uplo), n, alpha.re, x, c.incx, a, lda)
        else
            try l2.syr(al, ob.io, T, uploOf(c.uplo), n, alpha, x, c.incx, a, lda),
        .spr => if (cx)
            try l2.hpr(al, ob.io, T, uploOf(c.uplo), n, alpha.re, x, c.incx, a)
        else
            try l2.spr(al, ob.io, T, uploOf(c.uplo), n, alpha, x, c.incx, a),
        .syr2 => if (cx)
            try l2.her2(al, ob.io, T, uploOf(c.uplo), n, alpha, x, c.incx, y, c.incy, a, lda)
        else
            try l2.syr2(al, ob.io, T, uploOf(c.uplo), n, alpha, x, c.incx, y, c.incy, a, lda),
        .spr2 => if (cx)
            try l2.hpr2(al, ob.io, T, uploOf(c.uplo), n, alpha, x, c.incx, y, c.incy, a)
        else
            try l2.spr2(al, ob.io, T, uploOf(c.uplo), n, alpha, x, c.incx, y, c.incy, a),
    }
}

/// The running state of one routine's sweep: calls made and the largest
/// test ratio.
pub const Sweep = struct {
    nc: usize = 0,
    errmax: f64 = 0,

    /// Makes the call, failing on an error exit (the Fortran's `.NOT.OK`).
    pub fn call(s: *Sweep, comptime T: type, comptime r: Routine, c: Call, w: anytype) !void {
        s.nc += 1;
        run(T, r, c, w) catch |e| {
            ob.say("error exit taken on a valid call: ", .{});
            c.describe(T, r, s.nc);
            return e;
        };
    }

    /// Fails if `ok` is false: argument `what` was changed.
    pub fn same(s: *Sweep, comptime T: type, r: Routine, c: Call, what: []const u8, ok: bool) !void {
        if (ok) return;
        ob.say("{s} was changed incorrectly: ", .{what});
        c.describe(T, r, s.nc);
        return error.TestUnexpectedResult;
    }

    pub fn result(s: *Sweep, comptime T: type, r: Routine, c: Call, res: Result) !void {
        s.errmax = @max(s.errmax, res.err);
        if (!res.fatal) return;
        ob.say("computed result is less than half accurate: ", .{});
        c.describe(T, r, s.nc);
        return error.TestUnexpectedResult;
    }

    pub fn finish(s: *Sweep, comptime T: type, r: Routine, p: Params) !void {
        if (s.errmax < p.thresh) return;
        ob.say("{c}{s} completed {d} calls with maximum test ratio {d}\n", .{ letter(T), name(T, r), s.nc, s.errmax });
        return error.TestUnexpectedResult;
    }
};
