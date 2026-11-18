//! The `inline` section (M4; only when named, any host): `rstd.linalg`'s
//! inline `@Vector` forms (`vector.dot` and `length`, `matrix.any_size`'s
//! `mul`, `det` and `solve`) against the BLAS and LAPACK calls that serve
//! past the inline limit, at sizes 4 to 16, `f32` and `f64`, one worker:
//! where the inline code stops winning, to confirm the limit of 8 or move it
//! (`src/linalg/design.md`, Dispatch).

const std = @import("std");
const Io = std.Io;
const rstd = @import("rstd");
const bench = @import("bench.zig");
const bb = @import("blas_bench.zig");
const libs = @import("libs.zig");

const linalg = rstd.linalg;
const blas = linalg.blas;
const lapack = linalg.lapack;
const fail = libs.fail;

const timing: bench.Timing = .{ .min_ns = 5 * std.time.ns_per_ms, .batches = 9 };

/// The races: `gemv` and `gemm` are `mul` by a vector and by a matrix; `det`
/// races `getrf`, `solve` `gesv` with one right-hand side.
const Op = enum { dot, length, gemv, gemm, det, solve };

pub fn run(io: Io, gpa: std.mem.Allocator) !void {
    try bb.clock(io, "inline start");
    std.debug.print(
        \\
        \\ns per call, best of {d} batches of at least {d} ms; one worker; the inline forms against the BLAS and LAPACK at the same size
        \\
        \\| op | T | n | inline | library | library/inline |
        \\| -- | - | -: | -----: | ------: | -------------: |
        \\
    , .{ timing.batches, @divTrunc(timing.min_ns, std.time.ns_per_ms) });
    inline for (.{ f32, f64 }) |T| {
        inline for (comptime std.meta.tags(Op)) |op| {
            inline for (4..17) |n| try race(io, gpa, T, n, op);
        }
    }
}

fn race(io: Io, gpa: std.mem.Allocator, comptime T: type, comptime n: usize, comptime op: Op) !void {
    var in: Inline(T, n, op) = .{};
    var lib: Library(T, n, op) = .{ .gpa = gpa, .io = io };
    fill(T, n, &lib.a, &lib.b, &lib.x);
    inline for (0..n) |j| {
        in.a.data[j] = lib.a[j * n ..][0..n].*;
        in.b.data[j] = lib.b[j * n ..][0..n].*;
    }
    in.x = lib.x;
    in.y = lib.x;
    lib.y = lib.x;
    // the library's errors surface here, before the race
    try lib.call();
    const ctxs = .{ &in, &lib };
    const calls = [_]bench.Call{ .of(ctxs[0]), .of(ctxs[1]) };
    var secs: [2]f64 = undefined;
    try bench.race(io, timing, &calls, &secs);
    if (libs.lib_err) |err| return err;
    std.debug.print("| {t} | {s} | {d} | {d:.1} | {d:.1} | {d:.2} |\n", .{ op, @typeName(T), n, secs[0] * 1e9, secs[1] * 1e9, secs[1] / secs[0] });
}

/// A diagonally dominant `a` (well conditioned for `det` and `solve`), `b`
/// and `x` from small integers, the same for every size.
fn fill(comptime T: type, comptime n: usize, a: *[n * n]T, b: *[n * n]T, x: *[n]T) void {
    for (0..n) |j| {
        for (0..n) |i| {
            const v: T = @floatFromInt(@as(i32, @intCast((i * 7 + j * 3) % 11)) - 5);
            a[i + j * n] = if (i == j) v + 4 * @as(T, @floatFromInt(n)) else v / 8;
            b[i + j * n] = @floatFromInt(@as(i32, @intCast((i * 5 + j * 11) % 13)) - 6);
        }
        x[j] = @floatFromInt(@as(i32, @intCast(j % 5)) - 2);
    }
}

fn Inline(comptime T: type, comptime n: usize, comptime op: Op) type {
    return struct {
        a: linalg.StaticMat(n, n, T) = undefined,
        b: linalg.StaticMat(n, n, T) = undefined,
        x: linalg.Vec(n, T) = undefined,
        y: linalg.Vec(n, T) = undefined,

        pub fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            const any = linalg.matrix.any_size;
            switch (op) {
                .dot => std.mem.doNotOptimizeAway(linalg.vector.dot(s.x, s.y)),
                .length => std.mem.doNotOptimizeAway(linalg.vector.length(s.x)),
                .gemv => std.mem.doNotOptimizeAway(any.mul(s.a, s.x)),
                .gemm => std.mem.doNotOptimizeAway(any.mul(s.a, s.b)),
                .det => std.mem.doNotOptimizeAway(any.det(s.a)),
                .solve => if (any.solve(s.a, s.x)) |r| std.mem.doNotOptimizeAway(r) else |e| fail(e),
            }
        }
    };
}

fn Library(comptime T: type, comptime n: usize, comptime op: Op) type {
    return struct {
        gpa: std.mem.Allocator,
        io: Io,
        a: [n * n]T = undefined,
        b: [n * n]T = undefined,
        c: [n * n]T = undefined,
        lu: [n * n]T = undefined,
        x: [n]T = undefined,
        y: [n]T = undefined,
        ipiv: [n]usize = undefined,

        fn call(s: *@This()) !void {
            const one = Io.failing;
            switch (op) {
                .dot => std.mem.doNotOptimizeAway(try blas.dot(s.io, T, .of(&s.x), .of(&s.y))),
                .length => std.mem.doNotOptimizeAway(blas.nrm2(s.io, T, .of(&s.x))),
                .gemv => try blas.gemv(s.gpa, one, T, .n, 1, &s.a, n, .of(&s.x), 0, .of(&s.y)),
                .gemm => try blas.gemm(s.gpa, one, T, .n, .n, n, n, n, 1, &s.a, n, &s.b, n, 0, &s.c, n),
                .det => {
                    s.lu = s.a;
                    _ = try lapack.getrf(s.gpa, one, T, n, n, &s.lu, n, &s.ipiv);
                    var d: T = 1;
                    for (0..n) |i| d *= s.lu[i * n + i];
                    std.mem.doNotOptimizeAway(d);
                },
                .solve => {
                    s.lu = s.a;
                    s.y = s.x;
                    try lapack.gesv(s.gpa, one, T, n, 1, &s.lu, n, &s.ipiv, &s.y, n);
                },
            }
        }

        pub fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            s.call() catch |e| fail(e);
        }
    };
}
