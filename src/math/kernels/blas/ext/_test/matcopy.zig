//! `omatcopy`'s tiles (`transpose.zig`) and bands (`interface.zig`)
//! against the plain loops (`matcopy.plain`), bit for bit.

const std = @import("std");
const sysinfo = @import("builtin");
const blas = @import("../../__root__.zig");
const threading = @import("../../threading.zig");
const types = @import("../../types.zig");
const check = @import("../../_test/check.zig");
const bitwise = @import("../../_test/bitwise.zig");
const pattern = @import("../../_test/pattern.zig");
const test_io = @import("../../_test/io.zig");
const matcopy = @import("../matcopy.zig");
const transpose = @import("../transpose.zig");
const interface = @import("../interface.zig");

const Io = std.Io;
const Trans = types.Trans;
const Copy = matcopy.Copy;
const isComplex = types.isComplex;
const tile = transpose.tile;

/// The same bits, or a NaN for a NaN: a compiler may spell one product two ways.
fn same(comptime T: type, w: T, g: T) bool {
    if (comptime isComplex(T)) return same(types.Float(T), w.re, g.re) and same(types.Float(T), w.im, g.im);
    return bitwise.same(T, w, g, true) or (std.math.isNan(w) and std.math.isNan(g));
}

fn expectSame(comptime T: type, want: []const T, got: []const T) !void {
    for (want, got, 0..) |w, g, i| if (!same(T, w, g)) {
        check.say("element {d}: ", .{i});
        return error.TestExpectedEqual;
    };
}

fn sentinel(comptime T: type) T {
    return if (comptime isComplex(T)) .init(7, 7) else 7;
}

/// `A` with special values, the padding past its rows too.
fn fillA(comptime T: type, a: []T, seed: u64) void {
    check.fill(T, a, seed);
    for (a, 0..) |*v, i| v.* = pattern.shape(T, v.*, i, .mixed);
}

/// Sizes around the block and the tile's edges; fewer on the rigs.
const dims = [_]usize{ 1, 3, 4, 5, tile - 3, tile - 1, tile, tile + 1, tile + 3, tile + 4, tile + 5, 2 * tile - 1, 2 * tile, 2 * tile + 1, 2 * tile + 3 };
const rig_dims = [_]usize{ 3, tile + 1, 2 * tile + 3 };

/// How a case lays its matrices out: leading dimensions rounded to the
/// block (which peels to a vector boundary) or not, and where in their
/// buffers `A` and `B` start.
const Layout = struct { round: bool, at_a: usize, at_b: usize };
const layouts = [_]Layout{ .{ .round = false, .at_a = 0, .at_b = 0 }, .{ .round = true, .at_a = 0, .at_b = 0 }, .{ .round = true, .at_a = 1, .at_b = 2 }, .{ .round = true, .at_a = 3, .at_b = 1 } };

fn ld(n: usize, round: bool) usize {
    return if (round) (n + 4) / 4 * 4 else n + 1;
}

fn tiledCase(comptime F: type, al: std.mem.Allocator, rows: usize, cols: usize, alpha: F, l: Layout) !void {
    const lda = ld(rows, l.round);
    const ldb = ld(cols, l.round);
    const a = try al.alloc(F, l.at_a + lda * cols);
    defer al.free(a);
    const b = try al.alloc(F, 2 * (l.at_b + ldb * rows));
    defer al.free(b);
    const want = b[0 .. b.len / 2];
    const got = b[b.len / 2 ..];
    fillA(F, a, rows * 131 + cols);
    @memset(b, sentinel(F));
    matcopy.plain(F, .ct, rows, cols, alpha, a.ptr + l.at_a, lda, want.ptr + l.at_b, ldb);
    _ = matcopy.Portable(F).omatcopy_k_ct(@intCast(rows), @intCast(cols), alpha, a.ptr + l.at_a, @intCast(lda), got.ptr + l.at_b, @intCast(ldb));
    // the padding and the ends too: nothing else is written
    try expectSame(F, want, got);
}

test "the tiled transposed copy gives the plain loops' bits at every size around its blocks and tiles, scaled or not, aligned or not" {
    const al = std.testing.allocator;
    const ds: []const usize = if (check.rig()) &rig_dims else &dims;
    inline for (.{ f32, f64 }) |F| {
        for (ds) |rows| for (ds) |cols| for ([_]F{ 1.5, 1, -1, 0 }) |alpha| for (layouts) |l| {
            tiledCase(F, al, rows, cols, alpha, l) catch |e| {
                check.say("{s} {d}x{d} alpha {d} layout {any}: ", .{ @typeName(F), rows, cols, alpha, l });
                return e;
            };
        };
    }
}

/// The kernel `trans` takes for `T`, as `interface.zig` picks it.
fn copyOf(comptime T: type, trans: Trans) Copy {
    const c = isComplex(T) and trans.conjA();
    return if (trans.transposed()) (if (c) .ctc else .ct) else (if (c) .cnc else .cn);
}

fn bandedCase(comptime T: type, al: std.mem.Allocator, tio: ?Io, trans: Trans, rows: usize, cols: usize, alpha: T) !void {
    const t = trans.transposed();
    const out_rows = if (t) cols else rows;
    const out_cols = if (t) rows else cols;
    const lda = rows + 3;
    const ldb = out_rows + 1;
    const a = try al.alloc(T, lda * cols);
    defer al.free(a);
    const b = try al.alloc(T, 2 * ldb * out_cols);
    defer al.free(b);
    const want = b[0 .. b.len / 2];
    const got = b[b.len / 2 ..];
    fillA(T, a, rows * 17 + cols);
    @memset(want, sentinel(T));
    switch (copyOf(T, trans)) {
        inline else => |op| matcopy.plain(T, op, rows, cols, alpha, a.ptr, lda, want.ptr, ldb),
    }
    for (1..10) |workers| {
        if (workers > 1 and tio == null) break;
        blas.config.cpu_count = @intCast(workers);
        @memset(got, sentinel(T));
        try interface.omatcopy(if (workers > 1) tio.? else Io.failing, T, trans, @intCast(rows), @intCast(cols), alpha, a.ptr, @intCast(lda), got.ptr, @intCast(ldb));
        expectSame(T, want, got) catch |e| {
            check.say("{d} workers: ", .{workers});
            return e;
        };
    }
}

test "omatcopy in bands gives the plain loops' bits for every letter, alpha and odd size, at 1 to 9 workers" {
    _ = blas.coreName();
    const cpus = blas.config.cpu_count;
    defer blas.config.cpu_count = cpus;
    // bands from a few hundred elements
    threading.test_serial_max = 512;
    defer threading.test_serial_max = null;
    const al = std.testing.allocator;
    var threaded: if (test_io.threaded_tests) test_io.TestIo else void = if (comptime test_io.threaded_tests) .init(al, .{}) else {};
    defer if (comptime test_io.threaded_tests) threaded.deinit();
    const tio: ?Io = if (comptime test_io.threaded_tests) threaded.io() else null;
    const sizes = [_][2]usize{ .{ 2 * tile + 3, 4 * tile + 3 }, .{ 4 * tile + 2, 2 * tile + 1 }, .{ tile + 1, 6 * tile + 9 }, .{ 200, tile - 1 }, .{ 3 * tile + 1, 2 * tile }, .{ 1, 300 }, .{ 300, 1 } };
    inline for (.{ f32, f64, types.Complex32 }) |T| {
        const alphas = if (comptime isComplex(T)) [_]T{ .init(1.5, -0.75), .init(1, 0), .init(0, 0) } else [_]T{ 1.5, 1, 0 };
        for ([_]Trans{ .n, .t, .r, .c }) |trans| for (sizes) |s| for (alphas) |alpha| {
            bandedCase(T, al, tio, trans, s[0], s[1], alpha) catch |e| {
                check.say("{s} {t} {d}x{d} alpha {any}: ", .{ @typeName(T), trans, s[0], s[1], alpha });
                return e;
            };
        };
    }
}

test "omatcopy's bands: the reference's threshold, real types only, transposed bands whole tiles at least 32 wide" {
    _ = blas.coreName();
    const cpus = blas.config.cpu_count;
    defer blas.config.cpu_count = cpus;
    try std.testing.expect(interface.threaded(f32) and interface.threaded(f64));
    inline for (.{ types.Complex32, types.Complex64, types.XDouble, f16 }) |T| {
        if (T != f64) try std.testing.expect(!interface.threaded(T));
    }
    if (comptime sysinfo.single_threaded) return;
    blas.config.cpu_count = 8;
    if (blas.config.threadCeiling()) |c| if (c < 2) return;
    inline for (.{ f32, f64 }) |T| {
        // `rows * cols < 65536 * 4` is serial
        try std.testing.expect(interface.bands(T, false, 511, 512) == null);
        try std.testing.expect(interface.bands(T, true, 512, 511) == null);
        try std.testing.expect(interface.bands(T, false, 512, 512) != null);
        try std.testing.expect(interface.bands(T, true, 512, 512) != null);
        // too few columns for two bands, whatever the size
        try std.testing.expect(interface.bands(T, true, 100_000, 2 * interface.band_min - 1) == null);
        try std.testing.expect(interface.bands(T, false, 100_000, 2 * interface.band_min - 1) != null);
        for ([_][2]usize{ .{ 512, 512 }, .{ 100_000, 64 }, .{ 100_000, 100 }, .{ 4097, 4099 }, .{ 70_000, 333 } }) |s| {
            const p = interface.bands(T, true, s[0], s[1]).?;
            try std.testing.expect(p.chunks.count <= s[1] / interface.band_min);
            for (0..p.chunks.count) |i| {
                const ch = p.chunks.get(i);
                try std.testing.expect(ch.start % tile == 0);
                try std.testing.expect(ch.end - ch.start >= interface.band_min);
            }
            try std.testing.expectEqual(s[1], p.chunks.get(p.chunks.count - 1).end);
        }
    }
    // one worker: no bands
    blas.config.cpu_count = 1;
    try std.testing.expect(interface.bands(f64, true, 4096, 4096) == null);
}
