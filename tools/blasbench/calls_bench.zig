//! The `calls` section (T210, issue #594; only when named, any host): the
//! fixed cost of a level-1 reduction's call, ns a call at n = 16 to 1000:
//! its kernel through `blas.kernel`, the public routine (`blas.sum`,
//! `sumWith(.fast)`, `blas.dot`) and, when loaded, OpenBLAS's `cblas_*`.

const std = @import("std");
const Io = std.Io;
const rstd = @import("rstd");
const bench = @import("bench.zig");
const bb = @import("blas_bench.zig");
const libs = @import("libs.zig");

const blas = rstd.math.kernels.blas;
const fail = libs.fail;

const timing: bench.Timing = .{ .min_ns = 20 * std.time.ns_per_ms, .batches = 7 };
const sizes = [_]usize{ 16, 100, 1000 };

/// A reduction raced: its slot, the public call and OpenBLAS's routine.
const Case = struct {
    name: []const u8,
    T: type,
    slot: blas.Slot,
    fast: bool,
    ob: [:0]const u8,
};

const cases = [_]Case{
    .{ .name = "f64 sum", .T = f64, .slot = .dsum_k, .fast = false, .ob = "cblas_dsum" },
    .{ .name = "f64 dot", .T = f64, .slot = .ddot_k, .fast = false, .ob = "cblas_ddot" },
    .{ .name = "f32 sum", .T = f32, .slot = .ssum_wide_k, .fast = false, .ob = "cblas_ssum" },
    .{ .name = "f32 sum, fast", .T = f32, .slot = .ssum_k, .fast = true, .ob = "cblas_ssum" },
};

pub fn run(io: Io, gpa: std.mem.Allocator) !void {
    _ = io;
    try bb.clock(libs.rstd_io, "calls start");
    std.debug.print(
        \\
        \\ns a call, best of {d} batches of at least {d} ms; rstd's `Io` ({s}); `kernel` is the slot's kernel alone, `call` the public routine
        \\
        \\| case | n | kernel | call | openblas | call - kernel |
        \\| ---- | -: | -----: | ---: | -------: | ------------: |
        \\
    , .{ timing.batches, @divTrunc(timing.min_ns, std.time.ns_per_ms), libs.io_note });
    inline for (cases) |c| {
        for (sizes) |n| try race(c, gpa, n);
    }
}

fn race(comptime c: Case, gpa: std.mem.Allocator, n: usize) !void {
    const T = c.T;
    const xy = try gpa.alloc(T, 2 * n);
    defer gpa.free(xy);
    var rng = std.Random.DefaultPrng.init(n);
    const r = rng.random();
    for (xy) |*v| v.* = r.float(T) * 2 - 1;
    var k: Kernel(c) = .{ .x = xy[0..n], .y = xy[n..] };
    var p: Public(c) = .{ .x = xy[0..n], .y = xy[n..] };
    const O = Ob(c);
    O.f = libs.lookup(.openblas, *const O.Sig, c.ob);
    var ob: O = .{ .x = xy[0..n], .y = xy[n..] };
    const calls = [_]bench.Call{ .of(&k), .of(&p), .of(&ob) };
    const count: usize = if (O.f != null) 3 else 2;
    var secs: [3]f64 = @splat(0);
    try bench.race(libs.rstd_io, timing, calls[0..count], &secs);
    if (libs.lib_err) |err| return err;
    const ns = [3]f64{ secs[0] * 1e9, secs[1] * 1e9, secs[2] * 1e9 };
    std.debug.print("| {s} | {d} | {d:.1} | {d:.1} | ", .{ c.name, n, ns[0], ns[1] });
    if (count == 3) std.debug.print("{d:.1} | ", .{ns[2]}) else std.debug.print("- | ", .{});
    std.debug.print("{d:.1} |\n", .{ns[1] - ns[0]});
}

fn isDot(comptime c: Case) bool {
    return std.mem.endsWith(u8, @tagName(c.slot), "dot_k");
}

/// The slot's kernel, called as the routines call it.
fn Kernel(comptime c: Case) type {
    return struct {
        x: []const c.T,
        y: []const c.T,

        pub fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            const kern = blas.kernel(c.slot);
            const x: [*]c.T = @constCast(s.x.ptr);
            const n: isize = @intCast(s.x.len);
            const r = if (comptime isDot(c)) kern(n, x, 1, @constCast(s.y.ptr), 1) else kern(n, x, 1);
            std.mem.doNotOptimizeAway(r);
        }
    };
}

/// The public routine.
fn Public(comptime c: Case) type {
    return struct {
        x: []const c.T,
        y: []const c.T,

        pub fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            const io = libs.rstd_io;
            const r = if (comptime isDot(c))
                blas.dot(io, c.T, .of(s.x), .of(s.y)) catch |e| return fail(e)
            else if (c.fast)
                blas.sumWith(io, c.T, .of(s.x), .fast) catch |e| return fail(e)
            else
                blas.sum(io, c.T, .of(s.x)) catch |e| return fail(e);
            std.mem.doNotOptimizeAway(r);
        }
    };
}

/// OpenBLAS's routine.
fn Ob(comptime c: Case) type {
    return struct {
        x: []const c.T,
        y: []const c.T,

        const Sig = if (isDot(c)) fn (c_int, *const anyopaque, c_int, *const anyopaque, c_int) callconv(.c) c.T else fn (c_int, *const anyopaque, c_int) callconv(.c) c.T;
        var f: ?*const Sig = null;

        pub fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            const n: c_int = @intCast(s.x.len);
            const g = f.?;
            const r = if (comptime isDot(c)) g(n, s.x.ptr, 1, s.y.ptr, 1) else g(n, s.x.ptr, 1);
            std.mem.doNotOptimizeAway(r);
        }
    };
}
