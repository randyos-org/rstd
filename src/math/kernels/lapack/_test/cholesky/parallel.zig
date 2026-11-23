//! Tests of `potrf`'s parallel algorithm: its factor, its failing minor, one worker's bits at any worker count, and `crew`.

const std = @import("std");
const blas = @import("../../../blas/__root__.zig");
const types = @import("../../../blas/types.zig");
const driver = @import("../../../blas/blas3/driver.zig");
const accel = @import("../../../blas/apple_accel.zig");
const Vtable = @import("../../../blas/Vtable.zig");
const interface = @import("../../interface.zig");
const cholesky = @import("../../cholesky.zig");
const check = @import("../../../blas/_test/check.zig");
const guard = @import("../../../blas/_test/guard.zig");
const lcheck = @import("../check.zig");
const fixtures = @import("../fixtures.zig");
const ch = @import("../cholesky.zig");

const Io = std.Io;
const factorParallel = cholesky.factorParallel;
const potrf = cholesky.potrf;
const potri = cholesky.potri;
const testing = std.testing;
const test_types = lcheck.test_types;
const fromF64 = lcheck.fromF64;
const expectSameBits = lcheck.expectSameBits;
const own_sizes = ch.own_sizes;
const shrunk_sizes = ch.shrunk_sizes;

pub fn parallelCases(vt: *const Vtable) anyerror!void {
    return lcheck.eachBlocking(vt, lcheck.shrink, test_types, parallelOn);
}

pub fn parallelOn(comptime T: type, v: *const Vtable, vi: usize) anyerror!void {
    const al = testing.allocator;
    for ([_]bool{ true, false }) |upper| {
        const blk = driver.Blocking.of(T, v);
        const sizes: []const usize = if (vi == 0) &shrunk_sizes else &own_sizes;
        for (sizes) |n| {
            if (check.rig() and n > 70) continue;
            var p = try fixtures.CholProblem(T).init(al, n, upper, n * 29 + @as(u64, @intFromBool(upper)));
            defer p.deinit(al);
            var poison: guard.Poison = .{ .child = al };
            const k = try factorParallel(poison.allocator(), check.io, T, interface.table(T), upper, blk, n, p.a.ptr, p.lda, 1, 1, 1);
            try testing.expectEqual(@as(?usize, null), k);
            p.expectFactor() catch |e| {
                check.say("{s}, blocking {s}, parallel algorithm\n", .{ @typeName(T), if (vi == 0) "shrunk" else "own" });
                return e;
            };
            if (n < 13) continue;
            // a failure past the first block
            p.store();
            p.a[(n - 2) * (p.lda + 1)] = fromF64(T, -100);
            try testing.expectEqual(@as(?usize, n - 2), try factorParallel(al, check.io, T, interface.table(T), upper, blk, n, p.a.ptr, p.lda, 1, 1, 1));
        }
    }
}

test "potrf's parallel algorithm: A = L Lᴴ, and the failing minor, on every linked core" {
    try check.eachCore(parallelCases);
}

/// Inputs of a threaded-bits case: `A`, copied into each run's `out`.
pub fn Threaded(comptime T: type) type {
    return struct {
        p: *const fixtures.CholProblem(T),

        pub fn run(c: @This(), io: Io, want: u32, out: []T) anyerror!void {
            const vt = interface.table(T);
            @memcpy(out, c.p.a);
            _ = try factorParallel(testing.allocator, io, T, vt, c.p.upper, driver.Blocking.of(T, vt), c.p.n, out.ptr, c.p.lda, want, 1, 1);
        }

        /// The public `potrf`, its workers `blas.config`'s and its teams by
        /// `interface.team_grain`.
        pub fn grained(c: @This(), io: Io, want: u32, out: []T) anyerror!void {
            const saved = blas.config.cpu_count;
            defer blas.config.cpu_count = saved;
            blas.config.cpu_count = want;
            @memcpy(out, c.p.a);
            try potrf(testing.allocator, io, T, if (c.p.upper) .u else .l, c.p.n, out, c.p.lda);
        }

        /// The public `potri`, threading by `blas.config`'s workers.
        pub fn inverse(c: @This(), io: Io, want: u32, out: []T) anyerror!void {
            const saved = blas.config.cpu_count;
            defer blas.config.cpu_count = saved;
            blas.config.cpu_count = want;
            @memcpy(out, c.p.a);
            try potri(testing.allocator, io, T, if (c.p.upper) .u else .l, c.p.n, out, c.p.lda);
        }
    };
}

test "potrf's parallel algorithm and potri give one worker's bits at any worker count" {
    _ = blas.coreName();
    const was = accel.routing();
    _ = accel.route(false);
    defer _ = accel.route(was);
    const al = testing.allocator;
    inline for (.{ f32, f64, types.Complex32, types.Complex64 }) |T| for ([_]bool{ true, false }) |upper| {
        // past every type's `potri` threshold
        const n: usize = if (check.rig()) 70 else 210;
        var p = try fixtures.CholProblem(T).init(al, n, upper, 41);
        defer p.deinit(al);
        const c: Threaded(T) = .{ .p = &p };
        try expectSameBits(T, c, Threaded(T).run, p.a.len, "potrf");
        _ = try p.factorIt(al);
        try expectSameBits(T, c, Threaded(T).inverse, p.a.len, "potri");
    };
}

test "potrf gives one worker's bits at any worker count where its grain allows teams" {
    _ = blas.coreName();
    const was = accel.routing();
    _ = accel.route(false);
    defer _ = accel.route(was);
    const al = testing.allocator;
    inline for (.{ f32, f64, types.Complex32, types.Complex64 }) |T| for ([_]bool{ true, false }) |upper| {
        // the first `trsm` and `herk` take teams for every type, on Windows' grain too
        const n: usize = if (check.rig()) 70 else 420;
        var p = try fixtures.CholProblem(T).init(al, n, upper, 43);
        defer p.deinit(al);
        const c: Threaded(T) = .{ .p = &p };
        try expectSameBits(T, c, Threaded(T).grained, p.a.len, "potrf, grained");
    };
}

test "crew: one worker per grain, a complex product four, never more than asked" {
    const grain = interface.block_grain;
    try testing.expectEqual(@as(u32, 1), interface.crew(f64, 12, grain, 64, 64, 64));
    try testing.expectEqual(@as(u32, 1), interface.crew(f64, 12, grain, 100, 100, 199));
    try testing.expectEqual(@as(u32, 2), interface.crew(f64, 12, grain, 100, 100, 200));
    try testing.expectEqual(@as(u32, 4), interface.crew(types.Complex64, 12, grain, 100, 100, 100));
    try testing.expectEqual(@as(u32, 12), interface.crew(f32, 12, grain, 512, 512, 512));
    try testing.expectEqual(@as(u32, 1), interface.crew(f32, 1, grain, 512, 512, 512));
    try testing.expectEqual(@as(u32, 3), interface.crew(f32, 3, 0, 2, 2, 2));
    try testing.expectEqual(@as(u32, 5), interface.crew(f64, 5, grain, std.math.maxInt(u64), 2, 2));
}
