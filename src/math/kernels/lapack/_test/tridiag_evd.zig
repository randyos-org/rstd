//! Tests of `../tridiag_evd.zig`.

const std = @import("std");
const ops = @import("../ops.zig");
const troot = @import("../tridiag_evd.zig");

const Allocator = std.mem.Allocator;
const divideAndConquer = troot.divideAndConquer;
const qrAlgorithm = troot.qrAlgorithm;

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../../blas/_test/check.zig");

const Solver = enum { qr, dc4, dc8 };

/// `U diag(s) Uᵀ` against the tridiagonal per element, faer's
/// `ApproxEq::eps()` (`128 eps`) times `max(1, ‖T‖)`; `U` orthogonal
/// alike, `s` nondecreasing.
fn expectTridiag(comptime F: type, al: Allocator, diag: []const F, off: []const F, solver: Solver) !void {
    const n = diag.len;
    const d = try al.dupe(F, diag);
    defer al.free(d);
    const e = try al.alloc(F, n);
    defer al.free(e);
    @memcpy(e[0 .. n - 1], off);
    e[n - 1] = 0;
    const u = try ops.alloc(F, al, n, n);
    defer ops.free(F, al, u);
    u.fill(std.math.nan(F));
    switch (solver) {
        .qr => try qrAlgorithm(F, d, e, u),
        .dc4 => try divideAndConquer(al, check.io, F, d, e, u, 4),
        .dc8 => try divideAndConquer(al, check.io, F, d, e, u, 8),
    }
    var tmax: F = 0;
    for (diag) |x| tmax = @max(tmax, @abs(x));
    for (off) |x| tmax = @max(tmax, @abs(x));
    const tol = 128 * std.math.floatEps(F) * @max(tmax, 1);
    for (1..n) |i| try testing.expect(d[i - 1] <= d[i]);
    for (0..n) |j| for (0..n) |i| {
        var s: F = 0;
        var g: F = 0;
        for (0..n) |k| {
            s += u.get(i, k) * d[k] * u.get(j, k);
            g += u.get(k, i) * u.get(k, j);
        }
        const want: F = if (i == j) diag[j] else if (i == j + 1) off[j] else if (j == i + 1) off[i] else 0;
        const id: F = if (i == j) 1 else 0;
        if (!(@abs(s - want) <= tol) or !(@abs(g - id) <= tol * @sqrt(@as(F, @floatFromInt(n))))) {
            check.say("tridiag evd {s} {t} n={d}: ({d},{d}) U S Ut {e} want {e}, Ut U {e}\n", .{ @typeName(F), solver, n, i, j, s, want, g });
            return error.TestUnexpectedResult;
        }
    };
}

test "faer tridiag evd: faer's evd_qr_tests through QR iteration and divide and conquer" {
    const al = testing.allocator;
    const x = 1e-6;
    const cases = [_]struct { d: []const f64, e: []const f64 }{
        .{ .d = &.{ 1, 1 }, .e = &.{0} },
        .{ .d = &.{ 1, 1 }, .e = &.{0.5213289} },
        .{ .d = &.{ 1.79069356, 1.20930644, 1 }, .e = &.{ -4.06813537e-01, 0 } },
        .{ .d = &.{ 1.95069537, 2.44845332, 2.56957029, 3.03128102, 1 }, .e = &.{ -7.02200909e-01, -1.11661820e+00, -6.81418803e-01, 0 } },
        .{ .d = &.{ 3, 2, 1, 0, 1, 2, 3 }, .e = &.{ 1, 1, 1, 1, 1, 1 } },
        .{ .d = &.{ 3, 2, 1, 0, 1, 2, 3, 3, 2, 1, 0, 1, 2, 3 }, .e = &.{ 1, 1, 1, 1, 1, 1, x, 1, 1, 1, 1, 1, 1 } },
        .{ .d = &.{ 0, 0, 1.0769230769230773, -0.4290761869709236, -0.8278050499098524, 0.07994922044020283, -0.35579623371016944, 0.6487378508167678, -0.9347442346214521, -0.08624745233962683, -0.4999243909534632, 1.3708277457481026, -0.2592167303689501, -0.5929351972647323, -0.5863220906879729, 0.15069873027683844, 0.2449309426221532, 0.5599151389028441, 0.440084861097156, 9.811634162559901e-17 }, .e = &.{ 1.7320508075688772, 2.081665999466133, 2.0303418353670932, 1.2463948607107287, 1.5895840148470526, 1.3810057029812097, 1.265168346300635, 0.8941431038915991, 1.007512301091709, 0.5877505835309086, 1.0370970338888965, 0.8628932798233644, 1.1935059937001073, 1.1614143449715744, 0.41040224297074174, 0.561318309959268, 3.1807090401145072e-15, 0.4963971959331084, 1.942890293094024e-16 } },
        .{ .d = &.{ 1.0000000000000002, 1.0000000000000002 }, .e = &.{7.216449660063518e-16} },
    };
    for (cases) |c| for ([_]Solver{ .qr, .dc4, .dc8 }) |s| try expectTridiag(f64, al, c.d, c.e, s);
}

test "faer tridiag evd: random, clustered, graded and zero tridiagonals, f32 and f64" {
    const al = testing.allocator;
    inline for (.{ f32, f64 }) |F| {
        for ([_]usize{ 3, 9, 33, 100 }) |n| for (0..5) |kind| {
            const d = try al.alloc(F, n);
            defer al.free(d);
            const e = try al.alloc(F, n - 1);
            defer al.free(e);
            check.fill(F, d, n + kind);
            check.fill(F, e, 2 * n + kind);
            switch (kind) {
                0 => {},
                // clustered: equal diagonal, tiny coupling
                1 => {
                    @memset(d, 1);
                    for (e) |*x| x.* *= std.math.floatEps(F);
                },
                // graded over the type's range of significant digits
                2 => for (d, 0..) |*x, i| {
                    x.* = std.math.pow(F, 10, -@as(F, @floatFromInt(i % 12)));
                },
                // a zero block in the middle
                3 => {
                    @memset(d[n / 3 .. n / 2], 0);
                    @memset(e[n / 3 .. n / 2], 0);
                },
                // all zero
                else => {
                    @memset(d, 0);
                    @memset(e, 0);
                },
            }
            for ([_]Solver{ .qr, .dc4, .dc8 }) |s| try expectTridiag(F, al, d, e, s);
        };
    }
}
