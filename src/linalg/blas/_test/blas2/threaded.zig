//! Threaded level 2 against one worker, bit for bit, on every linked core:
//! each routine at one worker, then at several on a threaded `Io` and on an
//! inline one, over chunks made small (`threading.test_serial_max`). The
//! routines that fold partial vectors are also checked against one chunk.

const std = @import("std");
const blas = @import("../../__root__.zig");
const threading = @import("../../threading.zig");
const types = @import("../../types.zig");
const Vtable = @import("../../Vtable.zig");
const check = @import("../blas1/check.zig");
const interface = @import("../../blas2/interface.zig");
const band = @import("../../blas2/banded.zig");
const tri = @import("../../blas2/triangular.zig");

const Io = std.Io;
const Allocator = std.mem.Allocator;
const Trans = types.Trans;
const isComplex = types.isComplex;
const Float = types.Float;

const all = .{ f32, f64, types.Complex32, types.Complex64 } ++ check.x_real ++ check.x_complex;

/// The long side of a test matrix: several 64-row chunks and a tail.
const long = 5 * 64 + 37;
/// The short side.
const short = 97;

fn val(comptime T: type, re: f64, im: f64) T {
    if (comptime isComplex(T)) return .init(@floatCast(re), @floatCast(im));
    return @floatCast(re);
}

const Body = fn (Allocator, Io, *std.ArrayList(u8)) anyerror!void;

/// `body`'s output at one worker against 3 and 8, on `tio` and inline.
fn sameBits(al: Allocator, tio: Io, comptime body: Body) !void {
    var want: std.ArrayList(u8) = .empty;
    defer want.deinit(al);
    blas.config.cpu_count = 1;
    try body(al, Io.failing, &want);
    for ([_]u32{ 3, 8 }) |cpus| for ([_]Io{ tio, Io.failing }) |io| {
        blas.config.cpu_count = cpus;
        var got: std.ArrayList(u8) = .empty;
        defer got.deinit(al);
        try body(al, io, &got);
        try std.testing.expectEqualSlices(u8, want.items, got.items);
    };
}

/// `body`'s output over chunks of `most` against one chunk, within rounding:
/// partial vectors change the order of the sums only.
fn nearSerial(comptime T: type, al: Allocator, comptime body: Body, most: u64) !void {
    blas.config.cpu_count = 1;
    var want: std.ArrayList(u8) = .empty;
    defer want.deinit(al);
    threading.test_serial_max = std.math.maxInt(u64);
    try body(al, Io.failing, &want);
    var got: std.ArrayList(u8) = .empty;
    defer got.deinit(al);
    threading.test_serial_max = most;
    try body(al, Io.failing, &got);
    try std.testing.expectEqual(want.items.len, got.items.len);
    const z = check.bitsLen(T);
    for (0..want.items.len / z) |i| {
        const w = check.fromBits(T, want.items[i * z ..][0..z]);
        const g = check.fromBits(T, got.items[i * z ..][0..z]);
        check.expectClose(T, w, g, check.tolerance(Float(T), 8 * long)) catch |e| {
            check.say("element {d}: ", .{i});
            return e;
        };
    }
}

/// A filled buffer of `len`.
fn filled(comptime T: type, al: Allocator, len: usize, seed: u64) ![]T {
    const v = try al.alloc(T, len);
    check.fill(T, v, seed);
    return v;
}

fn Cases(comptime T: type) type {
    return struct {
        const alpha = val(T, 1.25, -0.5);
        const beta = val(T, 0.5, 0.25);
        const incs = [_][2]isize{ .{ 1, 1 }, .{ -2, 3 } };
        /// Past `gemv_rows` rows twice over, and past 64 columns many times.
        const shapes = [_][2]usize{ .{ 2 * 4096 + 300, 9 }, .{ 9, 2 * 4096 + 300 }, .{ long, short }, .{ short, long } };
        const most = (2 * 4096 + 300) * 3;
        const transes = if (isComplex(T)) std.enums.values(Trans) else &[_]Trans{ .n, .t };

        fn gemv(al: Allocator, io: Io, out: *std.ArrayList(u8)) anyerror!void {
            const a = try filled(T, al, (long + 3) * long, 1);
            defer al.free(a);
            const x = try filled(T, al, most, 2);
            defer al.free(x);
            const y0 = try filled(T, al, most, 3);
            defer al.free(y0);
            const y = try al.alloc(T, y0.len);
            defer al.free(y);
            for (transes) |t| for (shapes) |s| for (incs) |inc| {
                @memcpy(y, y0);
                try interface.gemv(al, io, T, t, @intCast(s[0]), @intCast(s[1]), alpha, a.ptr, @intCast(s[0] + 3), x.ptr, inc[0], beta, y.ptr, inc[1]);
                try check.appendDigest(al, out, y);
            };
        }

        fn ger(al: Allocator, io: Io, out: *std.ArrayList(u8)) anyerror!void {
            const a0 = try filled(T, al, (long + 3) * long, 4);
            defer al.free(a0);
            const a = try al.alloc(T, a0.len);
            defer al.free(a);
            const x = try filled(T, al, most, 5);
            defer al.free(x);
            const y = try filled(T, al, most, 6);
            defer al.free(y);
            for (shapes) |s| for (incs) |inc| for ([_]bool{ false, true }) |holes| {
                // zeros in `y` split it into runs that cross the chunks
                if (holes) for ([_]usize{ 0, 50, 51, 130, 200, 201, 202, 300 }) |j| {
                    if (j < s[1]) y[j * @abs(inc[1])] = val(T, 0, 0);
                };
                @memcpy(a, a0);
                try interface.ger(al, io, T, @intCast(s[0]), @intCast(s[1]), alpha, x.ptr, inc[0], y.ptr, inc[1], a.ptr, @intCast(s[0] + 3));
                try check.appendDigest(al, out, a[0 .. (s[0] + 3) * s[1]]);
                if (comptime isComplex(T)) {
                    @memcpy(a, a0);
                    try interface.gerc(al, io, T, @intCast(s[0]), @intCast(s[1]), alpha, x.ptr, inc[0], y.ptr, inc[1], a.ptr, @intCast(s[0] + 3));
                    try check.appendDigest(al, out, a[0 .. (s[0] + 3) * s[1]]);
                }
            };
        }

        fn gbmv(al: Allocator, io: Io, out: *std.ArrayList(u8)) anyerror!void {
            const kl = 9;
            const ku = 8;
            const lda = kl + ku + 3;
            const a = try filled(T, al, lda * long, 10);
            defer al.free(a);
            const x = try filled(T, al, 3 * long, 11);
            defer al.free(x);
            const y0 = try filled(T, al, 3 * long, 12);
            defer al.free(y0);
            const y = try al.alloc(T, y0.len);
            defer al.free(y);
            for (transes) |t| for ([_][2]usize{ .{ long, short }, .{ short, long } }) |s| for (incs) |inc| {
                @memcpy(y, y0);
                try band.gbmv(al, io, T, t, @intCast(s[0]), @intCast(s[1]), kl, ku, alpha, a.ptr, lda, x.ptr, inc[0], beta, y.ptr, inc[1]);
                try check.appendBits(al, out, y);
            };
        }

        /// `symv`, and `hemv` over every triangle.
        fn symv(al: Allocator, io: Io, out: *std.ArrayList(u8)) anyerror!void {
            const n = long;
            const lda = n + 3;
            const a = try filled(T, al, lda * n, 13);
            defer al.free(a);
            const x = try filled(T, al, 3 * n, 14);
            defer al.free(x);
            const y0 = try filled(T, al, 3 * n, 15);
            defer al.free(y0);
            const y = try al.alloc(T, y0.len);
            defer al.free(y);
            for (std.enums.values(types.Uplo)) |u| for (incs) |inc| {
                if (!u.isRev()) {
                    @memcpy(y, y0);
                    try interface.symv(al, io, T, u, n, alpha, a.ptr, lda, x.ptr, inc[0], beta, y.ptr, inc[1]);
                    try check.appendBits(al, out, y);
                }
                if (comptime isComplex(T)) {
                    @memcpy(y, y0);
                    try interface.hemv(al, io, T, u, n, alpha, a.ptr, lda, x.ptr, inc[0], beta, y.ptr, inc[1]);
                    try check.appendBits(al, out, y);
                }
            };
        }

        /// `sbmv`, and `hbmv` over every triangle; `2k` past `gbmv`'s 15.
        fn sbmv(al: Allocator, io: Io, out: *std.ArrayList(u8)) anyerror!void {
            const n = long;
            const k = 9;
            const lda = k + 2;
            const a = try filled(T, al, lda * n, 16);
            defer al.free(a);
            const x = try filled(T, al, 3 * n, 17);
            defer al.free(x);
            const y0 = try filled(T, al, 3 * n, 18);
            defer al.free(y0);
            const y = try al.alloc(T, y0.len);
            defer al.free(y);
            for (std.enums.values(types.Uplo)) |u| for (incs) |inc| {
                if (!u.isRev()) {
                    @memcpy(y, y0);
                    try band.sbmv(al, io, T, u, n, k, alpha, a.ptr, lda, x.ptr, inc[0], beta, y.ptr, inc[1]);
                    try check.appendBits(al, out, y);
                }
                if (comptime isComplex(T)) {
                    @memcpy(y, y0);
                    try band.hbmv(al, io, T, u, n, k, alpha, a.ptr, lda, x.ptr, inc[0], beta, y.ptr, inc[1]);
                    try check.appendBits(al, out, y);
                }
            };
        }

        /// `spmv`, and `hpmv`.
        fn spmv(al: Allocator, io: Io, out: *std.ArrayList(u8)) anyerror!void {
            const n = long;
            const a = try filled(T, al, n * (n + 1) / 2, 19);
            defer al.free(a);
            const x = try filled(T, al, 3 * n, 20);
            defer al.free(x);
            const y0 = try filled(T, al, 3 * n, 21);
            defer al.free(y0);
            const y = try al.alloc(T, y0.len);
            defer al.free(y);
            for ([_]types.Uplo{ .u, .l }) |u| for (incs) |inc| {
                @memcpy(y, y0);
                try band.spmv(al, io, T, u, n, alpha, a.ptr, x.ptr, inc[0], beta, y.ptr, inc[1]);
                try check.appendBits(al, out, y);
                if (comptime isComplex(T)) {
                    @memcpy(y, y0);
                    try band.hpmv(al, io, T, u, n, alpha, a.ptr, x.ptr, inc[0], beta, y.ptr, inc[1]);
                    try check.appendBits(al, out, y);
                }
            };
        }

        /// `trmv`, `tbmv` and `tpmv` over every triangle, `trans` and
        /// diagonal.
        fn trmv(al: Allocator, io: Io, out: *std.ArrayList(u8)) anyerror!void {
            const n = long;
            const k = 9;
            const lda = n + 3;
            const a = try filled(T, al, lda * n, 22);
            defer al.free(a);
            const x0 = try filled(T, al, 3 * n, 23);
            defer al.free(x0);
            const x = try al.alloc(T, x0.len);
            defer al.free(x);
            const ops: []const Trans = if (isComplex(T)) &.{ .n, .t, .r, .c } else &.{ .n, .t };
            for ([_]types.Uplo{ .u, .l }) |u| for (ops) |t| for (std.enums.values(types.Diag)) |d| for (incs) |inc| {
                inline for (0..3) |r| {
                    @memcpy(x, x0);
                    switch (r) {
                        0 => try tri.trmv(al, io, T, u, t, d, n, a.ptr, lda, x.ptr, inc[0]),
                        1 => try tri.tbmv(al, io, T, u, t, d, n, k, a.ptr, k + 1, x.ptr, inc[0]),
                        else => try tri.tpmv(al, io, T, u, t, d, n, a.ptr, x.ptr, inc[0]),
                    }
                    try check.appendBits(al, out, x);
                }
            };
        }

        /// The rank updates, full and packed, rank 1 and 2, with zeros in
        /// `x` and `y` for the columns they skip.
        fn rank(al: Allocator, io: Io, out: *std.ArrayList(u8)) anyerror!void {
            const n = long;
            const lda = n + 3;
            const a0 = try filled(T, al, lda * n, 7);
            defer al.free(a0);
            const a = try al.alloc(T, a0.len);
            defer al.free(a);
            const x = try filled(T, al, 3 * n, 8);
            defer al.free(x);
            const y = try filled(T, al, 3 * n, 9);
            defer al.free(y);
            for ([_]usize{ 5, 77, 78, 200 }) |j| x[j] = val(T, 0, 0);
            for ([_]usize{ 77, 150 }) |j| y[j] = val(T, 0, 0);
            const ra = Float(T);
            for ([_]types.Uplo{ .u, .l }) |u| for (incs) |inc| {
                const ni: isize = n;
                const ld: isize = lda;
                inline for (0..4) |r| {
                    @memcpy(a, a0);
                    switch (r) {
                        0 => try band.syr(al, io, T, u, ni, alpha, x.ptr, inc[0], a.ptr, ld),
                        1 => try band.spr(al, io, T, u, ni, alpha, x.ptr, inc[0], a.ptr),
                        2 => try band.syr2(al, io, T, u, ni, alpha, x.ptr, inc[0], y.ptr, inc[1], a.ptr, ld),
                        else => try band.spr2(al, io, T, u, ni, alpha, x.ptr, inc[0], y.ptr, inc[1], a.ptr),
                    }
                    try check.appendDigest(al, out, a);
                    if (comptime isComplex(T)) {
                        @memcpy(a, a0);
                        switch (r) {
                            0 => try band.her(al, io, T, u, ni, @as(ra, 0.75), x.ptr, inc[0], a.ptr, ld),
                            1 => try band.hpr(al, io, T, u, ni, @as(ra, 0.75), x.ptr, inc[0], a.ptr),
                            2 => try band.her2(al, io, T, u, ni, alpha, x.ptr, inc[0], y.ptr, inc[1], a.ptr, ld),
                            else => try band.hpr2(al, io, T, u, ni, alpha, x.ptr, inc[0], y.ptr, inc[1], a.ptr),
                        }
                        try check.appendDigest(al, out, a);
                    }
                }
            };
        }
    };
}

fn onCore(vt: *const Vtable) !void {
    const saved = blas.testUse(vt);
    defer blas.testRestore(saved);
    const al = std.testing.allocator;
    var threaded: threading.TestIo = .init(al, .{});
    defer threaded.deinit();
    inline for (all) |T| {
        const C = Cases(T);
        // a dozen chunks or so of each; the last fold partial vectors
        const bodies = .{ C.gemv, C.ger, C.rank, C.gbmv, C.symv, C.sbmv, C.spmv, C.trmv };
        const names = .{ "gemv", "ger", "rank updates", "gbmv", "symv and hemv", "sbmv and hbmv", "spmv and hpmv", "trmv, tbmv and tpmv" };
        inline for (bodies, names, .{ 3000, 3000, 40, 3000, 40, 12000, 40, 40 }, 0..) |body, name, most, i| {
            threading.test_serial_max = most;
            sameBits(al, threaded.io(), body) catch |e| {
                check.say("{s} {s}: ", .{ name, @typeName(T) });
                return e;
            };
            if (i >= 3) nearSerial(T, al, body, most) catch |e| {
                check.say("{s} {s} against one chunk: ", .{ name, @typeName(T) });
                return e;
            };
        }
    }
}

test "threaded level 2 gives one worker's bits at any worker count, Io and core" {
    if (comptime !threading.threaded_tests) return error.SkipZigTest;
    _ = blas.coreName();
    const cpus = blas.config.cpu_count;
    defer blas.config.cpu_count = cpus;
    defer threading.test_serial_max = null;
    try check.eachCore(onCore);
}

/// One round of the race stress test: `gemv` by columns, `ger`, `syr2` and
/// `symv`, whose partial vectors fold, with the inputs filled once.
fn Stress(comptime T: type) type {
    return struct {
        const Self = @This();
        const alpha = val(T, 1.25, -0.5);
        const beta = val(T, 0.5, 0.25);
        /// `gemv`'s and `ger`'s rows and columns: 64-column chunks and a tail.
        const rows = 13;
        const cols = 3 * 64 + 37;
        /// `syr2`'s and `symv`'s order.
        const order = 96;
        const lda = order + 3;

        a: []T,
        x: []T,
        y: []T,
        /// Outputs; `lda` leaves guard bands below each column of `A`.
        out_a: []T,
        out_y: []T,

        fn init(al: Allocator) !Self {
            return .{
                .a = try filled(T, al, lda * cols, 31),
                .x = try filled(T, al, 2 * cols, 32),
                .y = try filled(T, al, 2 * cols, 33),
                .out_a = try al.alloc(T, lda * cols),
                .out_y = try al.alloc(T, 2 * cols),
            };
        }

        fn deinit(s: *Self, al: Allocator) void {
            inline for (.{ s.a, s.x, s.y, s.out_a, s.out_y }) |v| al.free(v);
        }

        fn round(s: *Self, al: Allocator, io: Io, out: *std.ArrayList(u8)) !void {
            // four chunks of `gemv` and `ger`, eight of `syr2` and `symv`
            threading.test_serial_max = 32 * rows;
            @memcpy(s.out_y, s.y);
            try interface.gemv(al, io, T, .t, rows, cols, alpha, s.a.ptr, lda, s.x.ptr, 1, beta, s.out_y.ptr, 2);
            try check.appendBits(al, out, s.out_y);
            @memcpy(s.out_a, s.a);
            try interface.ger(al, io, T, rows, cols, alpha, s.x.ptr, 2, s.y.ptr, 1, s.out_a.ptr, lda);
            try check.appendBits(al, out, s.out_a);
            threading.test_serial_max = order / 8 + 1;
            @memcpy(s.out_a, s.a);
            try band.syr2(al, io, T, .l, order, alpha, s.x.ptr, 1, s.y.ptr, 2, s.out_a.ptr, lda);
            try check.appendBits(al, out, s.out_a);
            @memcpy(s.out_y, s.y);
            try interface.symv(al, io, T, .u, order, alpha, s.a.ptr, lda, s.x.ptr, 1, beta, s.out_y.ptr, 2);
            try check.appendBits(al, out, s.out_y);
        }
    };
}

test "race stress: threaded level 2 keeps one worker's bits and its guard bands, round after round" {
    if (comptime !threading.threaded_tests) return error.SkipZigTest;
    _ = blas.coreName();
    const cpus = blas.config.cpu_count;
    defer blas.config.cpu_count = cpus;
    defer threading.test_serial_max = null;
    const al = std.testing.allocator;
    var threaded: threading.TestIo = .init(al, .{ .async_limit = .limited(8) });
    defer threaded.deinit();
    inline for (.{ f64, types.Complex32 }) |T| {
        var s: Stress(T) = try .init(al);
        defer s.deinit(al);
        var want: std.ArrayList(u8) = .empty;
        defer want.deinit(al);
        blas.config.cpu_count = 1;
        try s.round(al, Io.failing, &want);
        var got: std.ArrayList(u8) = .empty;
        defer got.deinit(al);
        for (0..threading.stress_rounds) |i| {
            blas.config.cpu_count = @intCast(2 + i % 7);
            got.clearRetainingCapacity();
            try s.round(al, threaded.io(), &got);
            if (!check.sameBytes(want.items, got.items)) std.testing.expectEqualSlices(u8, want.items, got.items) catch |e| {
                check.say("{s} round {d}: ", .{ @typeName(T), i });
                return e;
            };
        }
    }
}
