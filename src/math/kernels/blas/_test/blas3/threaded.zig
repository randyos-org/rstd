//! Threaded level 3 (`blas3/threaded.zig`) against the serial driver, bit for
//! bit, at 1 to 9 workers, on every `Io`.

const std = @import("std");
const zigstd = @import("../../../../../zigstd.zig");
const blas = @import("../../__root__.zig");
const Vtable = @import("../../Vtable.zig");
const driver = @import("../../blas3/driver.zig");
const threading = @import("../../threading.zig");
const test_io = @import("../io.zig");
const types = @import("../../types.zig");
const portable = @import("../../portable.zig");
const level3_kernels = @import("../../blas3/kernels.zig");
const check = @import("../check.zig");
const bitwise = @import("../bitwise.zig");
const guard = @import("../guard.zig");
const check3 = @import("check.zig");
const FpState = @import("../../FpState.zig");

const Io = std.Io;
const Blocking = driver.Blocking;
const fill = check.fillMixed;
const scalar = check.val;
const sameBits = bitwise.sameBits;

/// The portable kernels, blocking shrunk as `check.zig` shrinks it, so every
/// block and panel edge comes at small sizes; `r` leaves column tails inside.
pub fn shrunk(comptime T: type) Vtable {
    var vt: Vtable = .{};
    level3_kernels.installGeneric(&vt);
    const blk = Blocking.of(T, &vt);
    check3.shrinkBlocking(T, &vt, false, blk.um, blk.un);
    return vt;
}

/// The `op` pairs a type has: transposes, and for complex conjugates too.
pub fn transPairs(comptime T: type) []const [2]types.Trans {
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
pub const Depths = struct {
    // 32-bit words: 32-bit targets may lack 64-bit atomics (powerpc-power9)
    pub var seen: [8]std.atomic.Value(u32) = @splat(.init(0));

    pub fn note(k: types.BlasLong) void {
        const d: usize = @intCast(std.math.clamp(k, 0, 255));
        _ = seen[d / 32].fetchOr(@as(u32, 1) << @intCast(d % 32), .monotonic);
    }

    /// Whether the depths seen since the last `take` are `want`; clears them.
    pub fn same(want: [8]u32) bool {
        return std.mem.eql(u32, &want, &take());
    }

    /// The depths seen since the last `take`.
    pub fn take() [8]u32 {
        var out: [8]u32 = undefined;
        for (&seen, &out) |*s, *o| o.* = s.swap(0, .monotonic);
        return out;
    }

    /// Kernel `name` noting its depth, its third argument, then calling `real`.
    pub fn Spy(comptime name: []const u8) type {
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
    pub fn spied(comptime T: type, vt: Vtable) Vtable {
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

    pub fn patch(vt: *Vtable, comptime name: []const u8) void {
        if (@field(vt, name)) |real| {
            Spy(name).real = real;
            @field(vt, name) = Spy(name).wrap;
        }
    }
};

/// One shape against `driver.gemm`, bit for bit, at 1 to 9 workers on `io`,
/// and with the kernels' depths if `vt` is `Depths.spied`. Operand padding
/// holds NaN; `C`'s, a sentinel that must survive; scratch comes poisoned.
pub fn shape(comptime T: type, vt: *const Vtable, io: Io, tr: [2]types.Trans, m: usize, n: usize, k: usize, alpha: T, beta: T) !void {
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
pub fn shapes(comptime T: type, io: Io) !void {
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

pub const precisions = check.sd_types;

test "threaded gemm matches the serial driver's bits and kernel depths at 1 to 9 workers, on every Io" {
    inline for (precisions) |T| try test_io.eachTestIo(T, shapes);
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
pub fn fpShapes(own: FpState, io: Io) !void {
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

test {
    _ = @import("threaded/cancel.zig");
    _ = @import("threaded/public.zig");
    _ = @import("threaded/routines.zig");
    _ = @import("threaded/stress.zig");
    _ = @import("threaded/team.zig");
}
