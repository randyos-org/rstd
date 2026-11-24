//! `blat3`'s calls: each routine's arguments, the call, and the sweep's bookkeeping.

const std = @import("std");
const ob = @import("../__root__.zig");
const l3 = @import("../../../blas3/core.zig");
const guard = @import("../../guard.zig");
const b3 = @import("../blat3.zig");
const mk = @import("make.zig");

const Float = ob.Float;
const testing = std.testing;
const Params = b3.Params;
const Routine = b3.Routine;
const letter = b3.letter;
const nmax = b3.nmax;
const Result = mk.Result;
const re = mk.re;

// ---- the calls -----------------------------------------------------------------

pub fn Work(comptime T: type) type {
    const R = ob.Ref(T);
    return struct {
        a: []T,
        aa: []T,
        as: []T,
        b: []T,
        bb: []T,
        bs: []T,
        c: []T,
        cc: []T,
        cs: []T,
        ab: []T,
        w: []T,
        ct: []R,
        g: []Float(R),

        pub fn init(al: std.mem.Allocator) !@This() {
            const nn = nmax * nmax;
            return .{
                .a = try al.alloc(T, nn),
                .aa = try al.alloc(T, nn),
                .as = try al.alloc(T, nn),
                .b = try al.alloc(T, nn),
                .bb = try al.alloc(T, nn),
                .bs = try al.alloc(T, nn),
                .c = try al.alloc(T, nn),
                .cc = try al.alloc(T, nn),
                .cs = try al.alloc(T, nn),
                .ab = try al.alloc(T, 2 * nn),
                .w = try al.alloc(T, 2 * nmax),
                .ct = try al.alloc(R, nmax),
                .g = try al.alloc(Float(R), nmax),
            };
        }
    };
}

pub const Call = struct {
    transa: u8 = 'N',
    transb: u8 = 'N',
    side: u8 = 'L',
    uplo: u8 = 'U',
    diag: u8 = 'N',
    m: usize = 0,
    n: usize = 0,
    k: usize = 0,
    lda: usize = 0,
    ldb: usize = 0,
    ldc: usize = 0,
    alpha: [2]f64 = .{ 0, 0 },
    beta: [2]f64 = .{ 0, 0 },

    pub fn describe(c: Call, comptime T: type, r: Routine, nc: usize) void {
        ob.say("{c}{s} call {d}: transa={c} transb={c} side={c} uplo={c} diag={c} m={d} n={d} k={d} lda={d} ldb={d} ldc={d} alpha=({d},{d}) beta=({d},{d})\n", .{
            letter(T), @tagName(r), nc, c.transa, c.transb, c.side, c.uplo, c.diag, c.m, c.n, c.k, c.lda, c.ldb, c.ldc, c.alpha[0], c.alpha[1], c.beta[0], c.beta[1],
        });
    }
};

pub const asTrans = ob.transOf;
pub const asSide = ob.sideOf;
pub const asUplo = ob.uploOf;
pub const asDiag = ob.diagOf;

pub fn run(comptime T: type, comptime r: Routine, c: Call, w: anytype) !void {
    var poison: guard.Poison = .{ .child = testing.allocator };
    const al = poison.allocator();
    const alpha = ob.val(T, c.alpha[0], c.alpha[1]);
    const beta = ob.val(T, c.beta[0], c.beta[1]);
    const m: isize = @intCast(c.m);
    const n: isize = @intCast(c.n);
    const k: isize = @intCast(c.k);
    const lda: isize = @intCast(c.lda);
    const ldb: isize = @intCast(c.ldb);
    const ldc: isize = @intCast(c.ldc);
    const a = w.aa.ptr;
    const b = w.bb.ptr;
    const cp = w.cc.ptr;
    switch (r) {
        .gemm, .gemm3m => try @field(l3, @tagName(r))(al, ob.io, T, asTrans(c.transa), asTrans(c.transb), m, n, k, alpha, a, lda, b, ldb, beta, cp, ldc),
        .symm, .hemm => try @field(l3, @tagName(r))(al, ob.io, T, asSide(c.side), asUplo(c.uplo), m, n, alpha, a, lda, b, ldb, beta, cp, ldc),
        .trmm, .trsm => try @field(l3, @tagName(r))(al, ob.io, T, asSide(c.side), asUplo(c.uplo), asTrans(c.transa), asDiag(c.diag), m, n, alpha, a, lda, b, ldb),
        .syrk => try l3.syrk(al, ob.io, T, asUplo(c.uplo), asTrans(c.transa), n, k, alpha, a, lda, beta, cp, ldc),
        .herk => try l3.herk(al, ob.io, T, asUplo(c.uplo), asTrans(c.transa), n, k, re(alpha), a, lda, re(beta), cp, ldc),
        .syr2k => try l3.syr2k(al, ob.io, T, asUplo(c.uplo), asTrans(c.transa), n, k, alpha, a, lda, b, ldb, beta, cp, ldc),
        .her2k => try l3.her2k(al, ob.io, T, asUplo(c.uplo), asTrans(c.transa), n, k, alpha, a, lda, b, ldb, re(beta), cp, ldc),
    }
}

/// One routine's sweep: calls made, and the largest test ratio.
pub const Sweep = struct {
    nc: usize = 0,
    errmax: f64 = 0,

    pub fn call(s: *Sweep, comptime T: type, comptime r: Routine, c: Call, w: anytype) !void {
        s.nc += 1;
        run(T, r, c, w) catch |e| {
            if (e == error.SkipZigTest) return e;
            ob.say("error exit taken on a valid call: ", .{});
            c.describe(T, r, s.nc);
            return e;
        };
    }

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
        ob.say("{c}{s} completed {d} calls with maximum test ratio {d}\n", .{ letter(T), @tagName(r), s.nc, s.errmax });
        return error.TestUnexpectedResult;
    }
};

pub fn bump(ld: usize) ?usize {
    return ob.bump(ld, nmax);
}

pub fn val(comptime T: type, v: [2]f64) T {
    return ob.val(T, v[0], v[1]);
}
