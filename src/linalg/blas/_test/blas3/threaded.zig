//! Threaded level 3 (`blas3/threaded.zig`) against the serial driver, bit for
//! bit, at 1 to 9 workers, on every `Io`.

const std = @import("std");
const zigstd = @import("../../../../zigstd.zig");
const blas = @import("../../__root__.zig");
const Vtable = @import("../../Vtable.zig");
const driver = @import("../../blas3/driver.zig");
const threading = @import("../../threading.zig");
const test_io = @import("../io.zig");
const types = @import("../../types.zig");
const apple_accel = @import("../../apple_accel.zig");
const portable = @import("../../portable.zig");
const level3_kernels = @import("../../blas3/kernels.zig");
const check = @import("../check.zig");
const bitwise = @import("../bitwise.zig");
const guard = @import("../guard.zig");
const check3 = @import("check.zig");
const interface = @import("../../blas3/interface.zig");
const FpState = @import("../../FpState.zig");
const workers = @import("../../blas3/threaded.zig").workers;

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Blocking = driver.Blocking;
const Error = driver.Error;

const fill = check.fillMixed;

const scalar = check.val;

const sameBits = bitwise.sameBits;

/// The portable kernels, blocking shrunk as `check.zig` shrinks it, so every
/// block and panel edge comes at small sizes; `r` leaves column tails inside.
fn shrunk(comptime T: type) Vtable {
    var vt: Vtable = .{};
    level3_kernels.installGeneric(&vt);
    const blk = Blocking.of(T, &vt);
    check3.shrinkBlocking(T, &vt, false, blk.um, blk.un);
    return vt;
}

/// The `op` pairs a type has: transposes, and for complex conjugates too.
fn transPairs(comptime T: type) []const [2]types.Trans {
    return comptime blk: {
        const all: []const types.Trans = if (types.isComplex(T)) &.{ .n, .t, .r, .c } else &.{ .n, .t };
        var out: [all.len * all.len][2]types.Trans = undefined;
        for (all, 0..) |x, i| for (all, 0..) |y, j| {
            out[i * all.len + j] = .{ x, y };
        };
        const final = out;
        break :blk &final;
    };
}

/// The depths the level 3 kernels see, a bit each (past 255 in the last):
/// `Q`'s slices, which must not follow the worker count (threading.md,
/// Guardrails).
const Depths = struct {
    // 32-bit words: 32-bit targets may lack 64-bit atomics (powerpc-power9)
    var seen: [8]std.atomic.Value(u32) = @splat(.init(0));

    fn note(k: types.BlasLong) void {
        const d: usize = @intCast(std.math.clamp(k, 0, 255));
        _ = seen[d / 32].fetchOr(@as(u32, 1) << @intCast(d % 32), .monotonic);
    }

    /// Whether the depths seen since the last `take` are `want`; clears them.
    fn same(want: [8]u32) bool {
        return std.mem.eql(u32, &want, &take());
    }

    /// The depths seen since the last `take`.
    fn take() [8]u32 {
        var out: [8]u32 = undefined;
        for (&seen, &out) |*s, *o| o.* = s.swap(0, .monotonic);
        return out;
    }

    /// Kernel `name` noting its depth, its third argument, then calling `real`.
    fn Spy(comptime name: []const u8) type {
        const F = @typeInfo(@FieldType(Vtable, name)).optional.child;
        const Fn = @typeInfo(F).pointer.child;
        const A = std.meta.ArgsTuple(Fn);
        const R = @typeInfo(Fn).@"fn".return_type.?;
        const cc = types.kernelCc(f64);
        // Zig 0.17 names it `fields`
        const info = @typeInfo(A).@"struct";
        const n = if (@hasField(@TypeOf(info), "field_names")) info.field_names.len else info.fields.len;
        return struct {
            var real: F = undefined;
            fn t(comptime i: usize) type {
                return @FieldType(A, std.fmt.comptimePrint("{d}", .{i}));
            }
            const wrap: F = switch (n) {
                8 => &struct {
                    fn f(a0: t(0), a1: t(1), k: t(2), a3: t(3), a4: t(4), a5: t(5), a6: t(6), a7: t(7)) callconv(cc) R {
                        note(k);
                        return real(a0, a1, k, a3, a4, a5, a6, a7);
                    }
                }.f,
                9 => &struct {
                    fn f(a0: t(0), a1: t(1), k: t(2), a3: t(3), a4: t(4), a5: t(5), a6: t(6), a7: t(7), a8: t(8)) callconv(cc) R {
                        note(k);
                        return real(a0, a1, k, a3, a4, a5, a6, a7, a8);
                    }
                }.f,
                10 => &struct {
                    fn f(a0: t(0), a1: t(1), k: t(2), a3: t(3), a4: t(4), a5: t(5), a6: t(6), a7: t(7), a8: t(8), a9: t(9)) callconv(cc) R {
                        note(k);
                        return real(a0, a1, k, a3, a4, a5, a6, a7, a8, a9);
                    }
                }.f,
                else => @compileError(name),
            };
        };
    }

    /// `vt` with `T`'s `gemm`, `trsm` and `trmm` kernels noting their depths.
    fn spied(comptime T: type, vt: Vtable) Vtable {
        const l = comptime portable.letter(T).?;
        const complex = comptime types.isComplex(T);
        const sides = .{ "LN", "LT", "RN", "RT" } ++ if (complex) .{ "LR", "LC", "RR", "RC" } else .{};
        const gemms = if (complex) .{ "gemm_kernel_n", "gemm_kernel_l", "gemm_kernel_r", "gemm_kernel_b" } else .{"gemm_kernel"};
        var out = vt;
        inline for (gemms) |s| patch(&out, l ++ s);
        inline for (sides) |s| {
            patch(&out, l ++ "trsm_kernel_" ++ s);
            patch(&out, l ++ "trmm_kernel_" ++ s);
        }
        return out;
    }

    fn patch(vt: *Vtable, comptime name: []const u8) void {
        if (@field(vt, name)) |real| {
            Spy(name).real = real;
            @field(vt, name) = Spy(name).wrap;
        }
    }
};

/// One shape against `driver.gemm`, bit for bit, at 1 to 9 workers on `io`,
/// and with the kernels' depths if `vt` is `Depths.spied`. Operand padding
/// holds NaN; `C`'s, a sentinel that must survive; scratch comes poisoned.
fn shape(comptime T: type, vt: *const Vtable, io: Io, tr: [2]types.Trans, m: usize, n: usize, k: usize, alpha: T, beta: T) !void {
    const al = std.testing.allocator;
    const ta = tr[0].transposed();
    const tb = tr[1].transposed();
    const conj: driver.Conj = .of(types.isComplex(T) and tr[0].conjA(), types.isComplex(T) and tr[1].conjA());
    const pad = 3;
    const ra = if (ta) k else m;
    const ca = if (ta) m else k;
    const rb = if (tb) n else k;
    const cb = if (tb) k else n;
    const a = try al.alloc(T, (ra + pad) * ca);
    defer al.free(a);
    const b = try al.alloc(T, (rb + pad) * cb);
    defer al.free(b);
    const c0 = try al.alloc(T, (m + pad) * n + pad);
    defer al.free(c0);
    const want = try al.alloc(T, c0.len);
    defer al.free(want);
    const got = try al.alloc(T, c0.len);
    defer al.free(got);
    fill(T, a, m * 31 + k);
    fill(T, b, n * 17 + k);
    fill(T, c0, m + n);
    const nan = scalar(T, std.math.nan(f64), std.math.nan(f64));
    for (0..ca) |j| @memset(a[j * (ra + pad) + ra ..][0..pad], nan);
    for (0..cb) |j| @memset(b[j * (rb + pad) + rb ..][0..pad], nan);
    const sentinel = scalar(T, 1234.5, -6.5);
    for (0..n) |j| @memset(c0[j * (m + pad) + m ..][0..pad], sentinel);
    @memset(c0[n * (m + pad) ..], sentinel);

    var poison: guard.Poison = .{ .child = al };
    const pa = poison.allocator();
    @memcpy(want, c0);
    try driver.gemm(pa, Io.failing, T, vt, ta, tb, conj, m, n, k, alpha, a.ptr, ra + pad, b.ptr, rb + pad, beta, want.ptr, m + pad, 1);
    const depths = Depths.take();
    for (1..10) |w| {
        @memcpy(got, c0);
        try driver.gemm(pa, io, T, vt, ta, tb, conj, m, n, k, alpha, a.ptr, ra + pad, b.ptr, rb + pad, beta, got.ptr, m + pad, @intCast(w));
        if (!sameBits(T, want, got) or !Depths.same(depths)) {
            std.debug.print("{s} w={d} {t}{t} m={d} n={d} k={d}\n", .{ @typeName(T), w, tr[0], tr[1], m, n, k });
            return error.TestUnexpectedResult;
        }
    }
}

/// Every `op` pair over shapes around every block and panel edge (complex,
/// fewer shapes over its sixteen pairs).
fn shapes(comptime T: type, io: Io) !void {
    const vt = Depths.spied(T, shrunk(T));
    const blk = Blocking.of(T, &vt);
    const um = blk.um;
    const un = blk.un;
    const rr = blk.r;
    const ms = [_]usize{ 1, um + 1, 2 * blk.p + um + 1, 4 * blk.p + 3, 9 * um + 5 };
    const ns = [_]usize{ 1, un + 1, 3 * un + 1, rr + 1, 2 * rr + 3, 5 * rr };
    const ks = [_]usize{ 1, blk.q + 1, 2 * blk.q + 3 };
    const step: usize = if (comptime types.isComplex(T)) 5 else 2;
    const alpha = scalar(T, -0.75, 0.25);
    const beta = scalar(T, 0.5, -0.125);
    for (transPairs(T)) |tr| {
        var x: usize = 0;
        for (ms) |m| for (ns) |n| for (ks) |k| {
            x += 1;
            if (x % step != 0) continue;
            try shape(T, &vt, io, tr, m, n, k, alpha, beta);
        };
        // beta 0 and 1, which the first slice folds differently
        try shape(T, &vt, io, tr, 9 * um + 5, 2 * rr + 3, 2 * blk.q + 3, scalar(T, 1.5, 0), scalar(T, 0, 0));
        try shape(T, &vt, io, tr, 9 * um + 5, 2 * rr + 3, 2 * blk.q + 3, scalar(T, 1, 0), scalar(T, 1, 0));
    }
}

const precisions = check.sd_types;

test "threaded gemm matches the serial driver's bits and kernel depths at 1 to 9 workers, on every Io" {
    inline for (precisions) |T| try test_io.eachTestIo(T, shapes);
}

const Routine = check3.Routine;

/// The threaded routines besides `gemm`.
const routines = [_]Routine{ .symm, .hemm, .syrk, .herk, .syr2k, .her2k, .trmm, .trsm };

/// `r` on `vt` at up to `want` workers: `C` (`trmm`'s and `trsm`'s `B`) is
/// `m x n` (`n x n` for the rank updates), the depth `k`; `opt` is the side
/// or `trans`, then `upper`, then for `trmm` and `trsm` `trans` (and for
/// complex, conjugating).
fn routine(comptime T: type, vt: *const Vtable, al: Allocator, io: Io, r: Routine, opt: [3]bool, m: usize, n: usize, k: usize, alpha: T, beta: T, a: [*]const T, b: [*]const T, ld: usize, c: [*]T, want: u32) !void {
    const tri = r == .trmm or r == .trsm;
    const rank = r == .syrk or r == .herk or r == .syr2k or r == .her2k;
    const case: check3.Case = .{ .r = r, .left = opt[0], .upper = opt[1], .ta = if (tri) .of(opt[2], opt[2]) else .of(rank and opt[0], false), .m = m, .n = n, .k = k };
    // `herk`'s `alpha` is real
    const al_s = if (comptime types.isComplex(T)) (if (r == .herk) scalar(T, alpha.re, 0) else alpha) else alpha;
    try check3.callAny(T, al, io, vt, r, case, a, ld, b, ld, al_s, beta, c, ld, want);
}

/// A tame triangle for `trsm` in `a` at `ld`: small off the diagonal, at
/// least 1 on it.
fn tame(comptime T: type, a: []T, ld: usize) void {
    for (0..ld) |j| for (0..ld) |i| {
        const v = &a[i + j * ld];
        const f = scalar(T, if (i == j) 1.5 else 0.25 / @as(f64, @floatFromInt(ld)), 0);
        v.* = if (comptime types.isComplex(T)) v.mul(f).add(if (i == j) scalar(T, 1, 0) else scalar(T, 0, 0)) else v.* * f + @as(T, if (i == j) 1 else 0);
    };
}

/// `r`'s shapes around every block, panel and triangle edge, all options,
/// against the serial loops bit for bit at 1 to 9 workers on `io`, and with
/// the kernels' depths if `vt` is `Depths.spied`.
fn routineShapes(comptime T: type, vt: *const Vtable, io: Io, r: Routine) !void {
    const al = std.testing.allocator;
    const blk = Blocking.of(T, vt);
    const herm = r == .herk or r == .her2k;
    const rank = r == .syrk or r == .herk or r == .syr2k or r == .her2k;
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
        if (!tri and opt[2]) continue;
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
fn shrunkCore(comptime T: type) Vtable {
    var vt = blas.vtable().*;
    const blk = Blocking.of(T, &vt);
    check3.shrinkBlocking(T, &vt, false, blk.um, blk.un);
    return vt;
}

test "threaded symm, hemm, rank-k updates, trmm and trsm match the serial loops' bits and kernel depths at 1 to 9 workers" {
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

test "threaded gemm keeps the serial bits under the caller's rounding toward zero and flush to zero" {
    if (comptime !test_io.threaded_tests or !FpState.governs) return error.SkipZigTest;
    const own: FpState = .capture();
    defer own.apply();
    var threaded: test_io.TestIo = .init(std.testing.allocator, .{});
    defer threaded.deinit();
    try fpShapes(own, threaded.io());
    if (comptime test_io.evented_workers) {
        var ev: threading.Evented = undefined;
        try zigstd.testExpectTrue(test_io.initEvented(&ev, 3));
        defer ev.deinit();
        try fpShapes(own, ev.io());
    }
}

/// `gemm` on `io` under `own`, then rounding toward zero.
fn fpShapes(own: FpState, io: Io) !void {
    inline for (precisions) |T| {
        const vt = shrunk(T);
        const blk = Blocking.of(T, &vt);
        // the workers start under the default state first
        for ([_]FpState{ own, own.towardZero() }) |s| {
            s.apply();
            try shape(T, &vt, io, .{ .n, .t }, 9 * blk.um + 5, 2 * blk.r + 3, 2 * blk.q + 3, scalar(T, -0.75, 0.25), scalar(T, 0.5, -0.125));
            try zigstd.testExpectTrue(FpState.capture().eql(s));
        }
    }
}

test "threaded gemm keeps the serial bits on a fiber Io with every peer on the caller's thread" {
    const Evented = threading.Evented;
    if (comptime !test_io.hosted or Evented == void) return error.SkipZigTest;
    // one thread, so a waiter that never yielded would hang; with worker
    // threads std's `Uring` and `Dispatch` crash in `Group` (issue #208)
    var ev: Evented = undefined;
    if (!test_io.initEvented(&ev, 0)) return error.SkipZigTest;
    defer ev.deinit();
    const vt = shrunk(f64);
    const blk = Blocking.of(f64, &vt);
    // `ic` and `jc` both past 1, several slices
    try shape(f64, &vt, ev.io(), .{ .n, .t }, 9 * blk.um + 5, 2 * blk.r + 3, 2 * blk.q + 3, -0.75, 0.5);
}

test "threaded gemm matches the serial driver on the selected core's own blocking" {
    if (comptime !test_io.threaded_tests) return error.SkipZigTest;
    var threaded: test_io.TestIo = .init(std.testing.allocator, .{});
    defer threaded.deinit();
    const vt = blas.vtable();
    inline for (precisions) |T| {
        const blk = Blocking.of(T, vt);
        try shape(T, vt, threaded.io(), .{ .n, .t }, 2 * blk.p + blk.um + 3, 3 * blk.un * 8 + 5, blk.q + 7, scalar(T, -0.75, 0.5), scalar(T, 0.5, 0));
        try shape(T, vt, threaded.io(), .{ .t, .n }, 3 * blk.p + 1, 97, 2 * blk.q + 1, scalar(T, 1.25, -1), scalar(T, 1, 0));
    }
}

/// Takes the public routines off Accelerate (macOS), so they run rstd's
/// threading; returns whether it was routing, for `apple_accel.route`.
fn rstdOnly() bool {
    // after init, which turns routing on
    _ = blas.vtable();
    const was = apple_accel.routing();
    _ = apple_accel.route(false);
    return was;
}

test "the public gemm threads above the threshold and keeps the serial bits" {
    if (comptime !test_io.threaded_tests) return error.SkipZigTest;
    _ = blas.coreName();
    const saved = blas.config.cpu_count;
    defer blas.config.cpu_count = saved;
    const accel = rstdOnly();
    defer _ = apple_accel.route(accel);
    var threaded: test_io.TestIo = .init(std.testing.allocator, .{});
    defer threaded.deinit();
    const al = std.testing.allocator;
    inline for (precisions) |T| {
        // past `m * n * k` of 8 serial units for every type's rule
        const m = 211;
        const n = 197;
        const k = 173;
        const a = try al.alloc(T, m * k);
        defer al.free(a);
        const b = try al.alloc(T, k * n);
        defer al.free(b);
        const want = try al.alloc(T, m * n);
        defer al.free(want);
        const got = try al.alloc(T, m * n);
        defer al.free(got);
        fill(T, a, 5);
        fill(T, b, 6);
        fill(T, want, 7);
        @memcpy(got, want);
        const alpha = scalar(T, 0.5, 2);
        const beta = scalar(T, -1.5, 0.25);
        try driver.gemm(al, Io.failing, T, blas.vtable(), true, false, .n, m, n, k, alpha, a.ptr, k, b.ptr, k, beta, want.ptr, m, 1);
        for ([_]u32{ 1, 3, 8 }) |cpus| {
            blas.config.cpu_count = cpus;
            try zigstd.testExpectEqual(cpus, workers(.gemm, T, m * n * k));
            const c = try al.dupe(T, got);
            defer al.free(c);
            try interface.gemm(al, threaded.io(), T, .t, .n, m, n, k, alpha, a.ptr, k, b.ptr, k, beta, c.ptr, m);
            try zigstd.testExpectTrue(sameBits(T, want, c));
        }
    }
}

test "the public symm, herk, syr2k and trmm thread above their thresholds and keep the serial bits" {
    if (comptime !test_io.threaded_tests) return error.SkipZigTest;
    // selecting the core fills `blas.config` once; override after
    _ = blas.coreName();
    const saved = blas.config.cpu_count;
    defer blas.config.cpu_count = saved;
    blas.config.cpu_count = 5;
    const accel = rstdOnly();
    defer _ = apple_accel.route(accel);
    var threaded: test_io.TestIo = .init(std.testing.allocator, .{});
    defer threaded.deinit();
    const io = threaded.io();
    const al = std.testing.allocator;
    const T = types.Complex64;
    const n = 157;
    const k = 141;
    const a = try al.alloc(T, n * n);
    defer al.free(a);
    const b = try al.alloc(T, n * n);
    defer al.free(b);
    const want = try al.alloc(T, n * n);
    defer al.free(want);
    const got = try al.alloc(T, n * n);
    defer al.free(got);
    fill(T, a, 21);
    fill(T, b, 22);
    fill(T, want, 23);
    @memcpy(got, want);
    const vt = blas.vtable();
    const alpha = scalar(T, 0.5, -1);
    try zigstd.testExpectEqual(@as(u32, 5), workers(.symm, T, 2 * n * n * n));
    try driver.symm(al, Io.failing, T, vt, false, false, true, n, n, alpha, a.ptr, n, b.ptr, n, alpha, want.ptr, n, 1);
    try interface.symm(al, io, T, .r, .u, n, n, alpha, a.ptr, n, b.ptr, n, alpha, got.ptr, n);
    try zigstd.testExpectTrue(sameBits(T, want, got));
    try zigstd.testExpectEqual(@as(u32, 5), workers(.syrk, T, (n + 1) * n * k));
    try driver.syrk(al, Io.failing, T, vt, true, false, true, n, k, scalar(T, 0.25, 0), a.ptr, k, null, 0, scalar(T, 2, 0), want.ptr, n, 1);
    try interface.herk(al, io, T, .l, .c, n, k, 0.25, a.ptr, k, 2, got.ptr, n);
    try zigstd.testExpectTrue(sameBits(T, want, got));
    try zigstd.testExpectEqual(@as(u32, 5), workers(.syr2k, T, n * k));
    try driver.syrk(al, Io.failing, T, vt, false, true, false, n, k, alpha, a.ptr, n, b.ptr, n, alpha, want.ptr, n, 1);
    try interface.syr2k(al, io, T, .u, .n, n, k, alpha, a.ptr, n, b.ptr, n, alpha, got.ptr, n);
    try zigstd.testExpectTrue(sameBits(T, want, got));
    try zigstd.testExpectEqual(@as(u32, 5), workers(.trsm, T, n * k));
    try driver.trmm(al, Io.failing, T, vt, false, false, true, true, true, k, n, alpha, a.ptr, n, want.ptr, n, 1);
    try interface.trmm(al, io, T, .r, .l, .c, .u, k, n, alpha, a.ptr, n, got.ptr, n);
    try zigstd.testExpectTrue(sameBits(T, want, got));
}

/// `r` (`gemm` if null) on `io` at 2 to 8 workers, round after round,
/// against one worker's bits on `Io.failing`. `C` (`B`) is `m x n` at
/// `ld`; every element past it is a guard band that must keep its value.
/// Scratch comes NaN-poisoned, so a read before its packing shows.
fn stress(comptime T: type, io: Io, r: ?Routine, opt: [3]bool, m: usize, n: usize, k: usize, rounds: usize) !void {
    const al = std.testing.allocator;
    const vt = shrunkCore(T);
    const ld = @max(m, n, k) + 3;
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
    fill(T, a, 41);
    fill(T, b, 42);
    fill(T, c0, 43);
    if (r == .trsm) tame(T, a, ld);
    var poison: guard.Poison = .{ .child = al };
    const pa = poison.allocator();
    const alpha = scalar(T, -0.75, 0.25);
    const beta = scalar(T, 0.5, 0);
    for (0..rounds + 1) |round| {
        const w: u32 = if (round == 0) 1 else @intCast(2 + round % 7);
        const out = if (round == 0) want else got;
        @memcpy(out, c0);
        const rio = if (round == 0) Io.failing else io;
        if (r) |rr| {
            try routine(T, &vt, pa, rio, rr, opt, m, n, k, alpha, beta, a.ptr, b.ptr, ld, out.ptr, w);
        } else {
            try driver.gemm(pa, rio, T, &vt, opt[0], opt[1], .n, m, n, k, alpha, a.ptr, ld, b.ptr, ld, beta, out.ptr, ld, w);
        }
        if (round > 0 and !sameBits(T, want, got)) {
            std.debug.print("{s} {?t} round {d} w={d}\n", .{ @typeName(T), r, round, w });
            return error.TestUnexpectedResult;
        }
    }
}

test "race stress: threaded level 3 keeps one worker's bits and its guard bands, round after round" {
    if (comptime !test_io.threaded_tests) return error.SkipZigTest;
    var threaded: test_io.TestIo = .init(std.testing.allocator, .{});
    defer threaded.deinit();
    try stresses(threaded.io());
    if (comptime test_io.evented_workers) {
        var ev: threading.Evented = undefined;
        try zigstd.testExpectTrue(test_io.initEvented(&ev, 3));
        defer ev.deinit();
        try stresses(ev.io());
    }
}

/// The race stress test's routines on `io`.
fn stresses(io: Io) !void {
    const rounds = test_io.stress_rounds;
    // shared packed `B` behind its barriers: `ic` and `jc` past 1, several slices
    {
        const blk = Blocking.of(f64, &shrunkCore(f64));
        try stress(f64, io, null, .{ false, true, false }, 9 * blk.um + 5, 2 * blk.r + 3, 2 * blk.q + 3, rounds);
    }
    // a triangle's balanced groups, both passes
    {
        const blk = Blocking.of(types.Complex32, &shrunkCore(types.Complex32));
        try stress(types.Complex32, io, .syr2k, .{ false, true, false }, 2 * blk.r + 3, 2 * blk.r + 3, 2 * blk.q + 3, rounds / 2);
    }
    // independent panels through `forEach`
    {
        const blk = Blocking.of(f32, &shrunkCore(f32));
        try stress(f32, io, .trsm, .{ true, false, false }, 2 * blk.r + 3, 2 * blk.r + 3, 1, rounds / 2);
    }
}

/// A kernel that holds every call until released, counting them and the
/// threads that make them.
const Hold = struct {
    const Kernel = @typeInfo(@FieldType(Vtable, "dgemm_kernel")).optional.child;
    const cc = types.kernelCc(f64);
    const BlasLong = types.BlasLong;
    const Atomic = std.atomic.Value(u32);

    var real: Kernel = undefined;
    var calls: Atomic = .init(0);
    var inside: Atomic = .init(0);
    var held: std.atomic.Value(bool) = .init(true);
    var first: std.atomic.Value(std.Thread.Id) = .init(0);
    var others: Atomic = .init(0);
    var settled: std.atomic.Value(bool) = .init(false);

    fn kernel(m: BlasLong, n: BlasLong, k: BlasLong, alpha: f64, a: [*]f64, b: [*]f64, c: [*]f64, ldc: BlasLong) callconv(cc) i32 {
        _ = calls.fetchAdd(1, .acq_rel);
        _ = inside.fetchAdd(1, .acq_rel);
        const id = std.Thread.getCurrentId();
        if (first.cmpxchgStrong(0, id, .acq_rel, .acquire)) |f| {
            if (f != id) _ = others.fetchAdd(1, .acq_rel);
        }
        defer _ = inside.fetchSub(1, .acq_rel);
        while (held.load(.acquire)) std.Thread.yield() catch {};
        return real(m, n, k, alpha, a, b, c, ldc);
    }

    /// Settles once a second thread is in the kernel (or after ten seconds
    /// without one), then releases the calls 50 ms later.
    fn release(io: Io) void {
        while (calls.load(.acquire) == 0) std.Thread.yield() catch {};
        // workers start well after the caller on a slow runner (issue #300)
        for (0..10_000) |_| {
            if (others.load(.acquire) > 0) break;
            io.sleep(.fromMilliseconds(1), .awake) catch {};
        }
        settled.store(true, .release);
        io.sleep(.fromMilliseconds(50), .awake) catch {};
        held.store(false, .release);
    }
};

fn canceledGemm(io: Io, vt: *const Vtable, a: []f64, b: []f64, c: []f64, n: usize, k: usize) Error!void {
    return driver.gemm(std.testing.allocator, io, f64, vt, false, false, .n, n, n, k, 1, a.ptr, n, b.ptr, k, 1, c.ptr, n, 4);
}

test "a canceled threaded gemm returns error.Canceled with no worker left running" {
    if (comptime !test_io.hosted) return error.SkipZigTest;
    var vt = shrunk(f64);
    Hold.real = vt.dgemm_kernel.?;
    vt.dgemm_kernel = &Hold.kernel;
    const n = 64;
    const k = 2000;
    const al = std.testing.allocator;
    const a = try al.alloc(f64, n * k);
    defer al.free(a);
    const b = try al.alloc(f64, k * n);
    defer al.free(b);
    const c = try al.alloc(f64, n * n);
    defer al.free(c);
    fill(f64, a, 1);
    fill(f64, b, 2);
    for ([_]Io.Limit{ .limited(3), .limited(8) }) |limit| {
        var threaded: Io.Threaded = .init(al, .{ .async_limit = limit, .concurrent_limit = limit });
        defer threaded.deinit();
        try cancelGemm(threaded.io(), &vt, a, b, c, n, k);
    }
    if (comptime test_io.evented_workers) {
        var ev: threading.Evented = undefined;
        try zigstd.testExpectTrue(test_io.initEvented(&ev, 3));
        defer ev.deinit();
        try cancelGemm(ev.io(), &vt, a, b, c, n, k);
    }
}

/// Cancels a held `gemm` on `io` once workers are in the kernel.
fn cancelGemm(io: Io, vt: *const Vtable, a: []f64, b: []f64, c: []f64, n: usize, k: usize) !void {
    Hold.calls.store(0, .release);
    Hold.held.store(true, .release);
    Hold.first.store(0, .release);
    Hold.others.store(0, .release);
    Hold.settled.store(false, .release);
    // a thread of its own sleeps on a blocking `Io`; a fiber `Io` needs a fiber
    const releaser = try std.Thread.spawn(.{}, Hold.release, .{Io.Threaded.global_single_threaded.io()});
    var future = io.concurrent(canceledGemm, .{ io, vt, a, b, c, n, k }) catch return error.SkipZigTest;
    while (!Hold.settled.load(.acquire)) std.Thread.yield() catch {};
    const result = future.cancel(io);
    releaser.join();
    try zigstd.testExpectError(error.Canceled, result);
    // workers did run beside the caller
    try zigstd.testExpectTrue(Hold.others.load(.acquire) > 0);
    // nothing calls the kernel once the call has returned
    const after = Hold.calls.load(.acquire);
    try zigstd.testExpectEqual(@as(u32, 0), Hold.inside.load(.acquire));
    io.sleep(.fromMilliseconds(20), .awake) catch {};
    try zigstd.testExpectEqual(after, Hold.calls.load(.acquire));
}
