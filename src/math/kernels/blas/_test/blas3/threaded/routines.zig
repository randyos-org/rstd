//! Threaded `symm`, `hemm`, the rank updates, `trmm` and `trsm` against the serial loops, bit for bit.

const std = @import("std");
const blas = @import("../../../__root__.zig");
const Vtable = @import("../../../Vtable.zig");
const driver = @import("../../../blas3/driver.zig");
const test_io = @import("../../io.zig");
const types = @import("../../../types.zig");
const check = @import("../../check.zig");
const bitwise = @import("../../bitwise.zig");
const guard = @import("../../guard.zig");
const check3 = @import("../check.zig");
const th = @import("../threaded.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Blocking = driver.Blocking;
const fill = check.fillMixed;
const scalar = check.val;
const sameBits = bitwise.sameBits;
const Depths = th.Depths;
const precisions = th.precisions;

pub const Routine = check3.Routine;

/// The threaded routines besides `gemm`.
pub const routines = [_]Routine{ .symm, .hemm, .syrk, .herk, .syr2k, .her2k, .gemmt, .trmm, .trsm };

/// `r` on `vt` at up to `want` workers: `C` (`trmm`'s and `trsm`'s `B`) is
/// `m x n` (`n x n` for the rank updates), the depth `k`; `opt` is the side
/// or `trans` (`gemmt`'s `transa`), then `upper`, then for `trmm` and
/// `trsm` `trans` (and for complex, conjugating), for `gemmt` `transb`.
pub fn routine(comptime T: type, vt: *const Vtable, al: Allocator, io: Io, r: Routine, opt: [3]bool, m: usize, n: usize, k: usize, alpha: T, beta: T, a: [*]const T, b: [*]const T, ld: usize, c: [*]T, want: u32) !void {
    const tri = r == .trmm or r == .trsm;
    const rank = r == .syrk or r == .herk or r == .syr2k or r == .her2k;
    const case: check3.Case = .{ .r = r, .left = opt[0], .upper = opt[1], .ta = if (tri) .of(opt[2], opt[2]) else .of((rank or r == .gemmt) and opt[0], false), .tb = .of(r == .gemmt and opt[2], false), .m = m, .n = n, .k = k };
    // `herk`'s `alpha` is real
    const al_s = if (comptime types.isComplex(T)) (if (r == .herk) scalar(T, alpha.re, 0) else alpha) else alpha;
    try check3.callAny(T, al, io, vt, r, case, a, ld, b, ld, al_s, beta, c, ld, want);
}

/// A tame triangle for `trsm` in `a` at `ld`: small off the diagonal, at
/// least 1 on it.
pub fn tame(comptime T: type, a: []T, ld: usize) void {
    for (0..ld) |j| for (0..ld) |i| {
        const v = &a[i + j * ld];
        const f = scalar(T, if (i == j) 1.5 else 0.25 / @as(f64, @floatFromInt(ld)), 0);
        v.* = if (comptime types.isComplex(T)) v.mul(f).add(if (i == j) scalar(T, 1, 0) else scalar(T, 0, 0)) else v.* * f + @as(T, if (i == j) 1 else 0);
    };
}

/// `r`'s shapes around every block, panel and triangle edge, all options,
/// against the serial loops bit for bit at 1 to 9 workers on `io`, and with
/// the kernels' depths if `vt` is `Depths.spied`.
pub fn routineShapes(comptime T: type, vt: *const Vtable, io: Io, r: Routine) !void {
    const al = std.testing.allocator;
    const blk = Blocking.of(T, vt);
    const herm = r == .herk or r == .her2k;
    // `gemmt` too has a square `C` and a depth
    const rank = r == .syrk or r == .herk or r == .syr2k or r == .her2k or r == .gemmt;
    const tri = r == .trmm or r == .trsm;
    const ns = [_]usize{ 1, blk.un + 1, 2 * blk.p + blk.um + 1, blk.r + 1, 2 * blk.r + 3 };
    const ks = [_]usize{ 1, blk.q + 1, 2 * blk.q + 3 };
    // past every side: `2 * p + um + 1` outgrows `r` on wide tiles (T142)
    const ld = std.mem.max(usize, &ns) + 2 * blk.q + 8;
    const a = try al.alloc(T, ld * ld);
    defer al.free(a);
    const b = try al.alloc(T, ld * ld);
    defer al.free(b);
    const c0 = try al.alloc(T, ld * ld);
    defer al.free(c0);
    const want = try al.alloc(T, c0.len);
    defer al.free(want);
    const got = try al.alloc(T, c0.len);
    defer al.free(got);
    fill(T, a, 11);
    fill(T, b, 12);
    fill(T, c0, 13);
    if (tri) tame(T, a, ld);
    var poison: guard.Poison = .{ .child = al };
    const pa = poison.allocator();
    const alpha = scalar(T, -0.75, 0.25);
    const beta = if (herm) scalar(T, 0.5, 0) else scalar(T, 0.5, -0.125);
    var x: usize = 0;
    for (0..8) |bits| {
        const opt: [3]bool = .{ bits & 1 != 0, bits & 2 != 0, bits & 4 != 0 };
        if (!tri and r != .gemmt and opt[2]) continue;
        for (ns) |m| for (ns) |n| for (ks) |k| {
            // the rank updates' `C` is square; the others have no depth
            if ((rank and m != n) or (!rank and k != 1)) continue;
            x += 1;
            if (x % 2 != 0) continue;
            @memcpy(want, c0);
            try routine(T, vt, pa, Io.failing, r, opt, m, n, k, alpha, beta, a.ptr, b.ptr, ld, want.ptr, 1);
            const depths = Depths.take();
            for (1..10) |w| {
                @memcpy(got, c0);
                try routine(T, vt, pa, io, r, opt, m, n, k, alpha, beta, a.ptr, b.ptr, ld, got.ptr, @intCast(w));
                if (!sameBits(T, want, got) or !Depths.same(depths)) {
                    std.debug.print("{s} {t} w={d} opt={any} m={d} n={d} k={d}\n", .{ @typeName(T), r, w, opt, m, n, k });
                    return error.TestUnexpectedResult;
                }
            }
        };
    }
}

/// The global vtable with `check.zig`'s shrunk blocking: the core's own
/// copies and kernels, every edge at small sizes.
pub fn shrunkCore(comptime T: type) Vtable {
    var vt = blas.vtable().*;
    const blk = Blocking.of(T, &vt);
    check3.shrinkBlocking(T, &vt, false, blk.um, blk.un);
    return vt;
}

test "threaded symm, hemm, rank-k updates, gemmt, trmm and trsm match the serial loops' bits and kernel depths at 1 to 9 workers" {
    inline for (precisions) |T| {
        const vt = Depths.spied(T, shrunkCore(T));
        for (routines) |r| {
            if (!types.isComplex(T) and (r == .hemm or r == .herk or r == .her2k)) continue;
            const Run = struct {
                vt: *const Vtable,
                r: Routine,
                fn f(run: @This(), io: Io) !void {
                    try routineShapes(T, run.vt, io, run.r);
                }
            };
            try test_io.eachTestIo(Run{ .vt = &vt, .r = r }, Run.f);
        }
    }
}
