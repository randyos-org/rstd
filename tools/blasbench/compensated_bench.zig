//! The `compensated` section (T138, T198; only when named, any host): `sum`,
//! `asum` and `dot` in each `blas.Summation` and, when loaded, OpenBLAS's,
//! `f32`, `f64` and `Complex(f32)` at n = 100 to 10^7, on rstd's `Io`
//! (`--threads`): the price of each mode, for choosing the defaults.

const std = @import("std");
const Io = std.Io;
const rstd = @import("rstd");
const bench = @import("bench.zig");
const bb = @import("blas_bench.zig");
const libs = @import("libs.zig");

const blas = rstd.linalg.blas;
const Summation = blas.Summation;
const fail = libs.fail;

const timing: bench.Timing = .{ .min_ns = 20 * std.time.ns_per_ms, .batches = 7 };
const sizes = [_]usize{ 100, 1000, 10_000, 100_000, 1_000_000, 10_000_000 };
const Op = enum { sum, asum, dot };
const modes = [_]Summation{ .plain, .compensated, .fast };

pub fn run(io: Io, gpa: std.mem.Allocator) !void {
    _ = io;
    try bb.clock(libs.rstd_io, "compensated start");
    std.debug.print(
        \\
        \\ns per element, best of {d} batches of at least {d} ms; rstd's `Io` ({s}); `openblas` is its `cblas_*` routine (`dotu_sub` for complex `dot`); ratios are over plain
        \\
        \\| type | n | op | plain | comp. | fast | openblas | comp./plain | fast/plain | ob/plain |
        \\| ---- | -: | -- | ----: | ----: | ---: | -------: | ----------: | ---------: | -------: |
        \\
    , .{ timing.batches, @divTrunc(timing.min_ns, std.time.ns_per_ms), libs.io_note });
    inline for (.{ f32, f64, std.math.Complex(f32) }) |T| {
        for (sizes) |n| {
            if (n > libs.n_max) continue;
            try race(T, gpa, n);
        }
    }
}

fn race(comptime T: type, gpa: std.mem.Allocator, n: usize) !void {
    const xy = try gpa.alloc(T, 2 * n);
    defer gpa.free(xy);
    var rng = std.Random.DefaultPrng.init(n);
    const r = rng.random();
    for (xy) |*v| v.* = if (comptime isComplex(T)) .init(r.float(f32) * 2 - 1, r.float(f32) * 2 - 1) else r.float(T) * 2 - 1;
    inline for (comptime std.meta.tags(Op)) |op| {
        var ctx: [modes.len]Call(T) = undefined;
        for (&ctx, modes) |*c, s| c.* = .{ .x = xy[0..n], .y = xy[n..], .op = op, .s = s };
        const O = Ob(T, op);
        O.f = libs.lookup(.openblas, *const O.Sig, O.name);
        var ob: O = .{ .x = xy[0..n], .y = xy[n..] };
        var calls: [modes.len + 1]bench.Call = undefined;
        for (calls[0..modes.len], &ctx) |*c, *x| c.* = .of(x);
        calls[modes.len] = .of(&ob);
        const count: usize = if (O.f != null) calls.len else modes.len;
        var secs: [calls.len]f64 = @splat(0);
        try bench.race(libs.rstd_io, timing, calls[0..count], &secs);
        if (libs.lib_err) |err| return err;
        var ns: [calls.len]f64 = undefined;
        for (&ns, secs) |*o, s| o.* = s * 1e9 / @as(f64, @floatFromInt(n));
        const name = if (comptime isComplex(T)) "c32" else @typeName(T);
        std.debug.print("| {s} | {d} | {t} | {d:.3} | {d:.3} | {d:.3} | ", .{ name, n, op, ns[0], ns[1], ns[2] });
        if (count > modes.len) std.debug.print("{d:.3} | ", .{ns[3]}) else std.debug.print("- | ", .{});
        std.debug.print("{d:.2} | {d:.2} | ", .{ ns[1] / ns[0], ns[2] / ns[0] });
        if (count > modes.len) std.debug.print("{d:.2} |\n", .{ns[3] / ns[0]}) else std.debug.print("- |\n", .{});
    }
}

fn isComplex(comptime T: type) bool {
    return @typeInfo(T) == .@"struct";
}

fn Call(comptime T: type) type {
    return struct {
        x: []const T,
        y: []const T,
        op: Op,
        s: Summation,

        pub fn run(o: *anyopaque) void {
            const c: *@This() = @ptrCast(@alignCast(o));
            const io = libs.rstd_io;
            const r = switch (c.op) {
                .sum => blas.sumWith(io, T, .of(c.x), c.s),
                .asum => blas.asumWith(io, T, .of(c.x), c.s),
                .dot => if (comptime isComplex(T)) blk: {
                    const d = blas.dotWith(io, T, .of(c.x), .of(c.y), c.s) catch |e| break :blk e;
                    break :blk d.re + d.im;
                } else blas.dotWith(io, T, .of(c.x), .of(c.y), c.s),
            };
            // kept, so the call is not dropped
            std.mem.doNotOptimizeAway(r catch |e| return fail(e));
        }
    };
}

/// OpenBLAS's `cblas_*` routine for `op` over `T`.
fn Ob(comptime T: type, comptime op: Op) type {
    const cx = isComplex(T);
    const F = if (cx) @FieldType(T, "re") else T;
    const re = if (F == f32) "s" else "d";
    const s = if (cx) (if (F == f32) "c" else "z") else re;
    return struct {
        x: []const T,
        y: []const T,

        const name = switch (op) {
            .sum, .asum => "cblas_" ++ (if (cx) re else "") ++ s ++ @tagName(op),
            .dot => "cblas_" ++ s ++ (if (cx) "dotu_sub" else "dot"),
        };
        const Sig = switch (op) {
            .sum, .asum => fn (c_int, *const anyopaque, c_int) callconv(.c) F,
            .dot => if (cx) fn (c_int, *const anyopaque, c_int, *const anyopaque, c_int, *anyopaque) callconv(.c) void else fn (c_int, *const anyopaque, c_int, *const anyopaque, c_int) callconv(.c) F,
        };
        var f: ?*const Sig = null;

        pub fn run(o: *anyopaque) void {
            const c: *@This() = @ptrCast(@alignCast(o));
            const n: c_int = @intCast(c.x.len);
            const g = f.?;
            const r = switch (op) {
                .sum, .asum => g(n, c.x.ptr, 1),
                .dot => if (cx) blk: {
                    var out: T = undefined;
                    g(n, c.x.ptr, 1, c.y.ptr, 1, &out);
                    break :blk out.re + out.im;
                } else g(n, c.x.ptr, 1, c.y.ptr, 1),
            };
            std.mem.doNotOptimizeAway(r);
        }
    };
}
