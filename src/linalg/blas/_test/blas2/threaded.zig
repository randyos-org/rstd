//! Threaded level 2 against one worker, bit for bit, on every linked core:
//! each routine at one worker, then at several on a threaded `Io` and on an
//! inline one, over chunks made small (`threading.test_serial_max`).

const std = @import("std");
const blas = @import("../../__root__.zig");
const threading = @import("../../threading.zig");
const types = @import("../../types.zig");
const Vtable = @import("../../Vtable.zig");
const check = @import("../blas1/check.zig");
const interface = @import("../../blas2/interface.zig");
const band = @import("../../blas2/banded.zig");

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
                try out.appendSlice(al, std.mem.sliceAsBytes(y));
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
                try out.appendSlice(al, std.mem.sliceAsBytes(a[0 .. (s[0] + 3) * s[1]]));
                if (comptime isComplex(T)) {
                    @memcpy(a, a0);
                    try interface.gerc(al, io, T, @intCast(s[0]), @intCast(s[1]), alpha, x.ptr, inc[0], y.ptr, inc[1], a.ptr, @intCast(s[0] + 3));
                    try out.appendSlice(al, std.mem.sliceAsBytes(a[0 .. (s[0] + 3) * s[1]]));
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
                try out.appendSlice(al, std.mem.sliceAsBytes(y));
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
                    try out.appendSlice(al, std.mem.sliceAsBytes(a));
                    if (comptime isComplex(T)) {
                        @memcpy(a, a0);
                        switch (r) {
                            0 => try band.her(al, io, T, u, ni, @as(ra, 0.75), x.ptr, inc[0], a.ptr, ld),
                            1 => try band.hpr(al, io, T, u, ni, @as(ra, 0.75), x.ptr, inc[0], a.ptr),
                            2 => try band.her2(al, io, T, u, ni, alpha, x.ptr, inc[0], y.ptr, inc[1], a.ptr, ld),
                            else => try band.hpr2(al, io, T, u, ni, alpha, x.ptr, inc[0], y.ptr, inc[1], a.ptr),
                        }
                        try out.appendSlice(al, std.mem.sliceAsBytes(a));
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
    var threaded: Io.Threaded = .init(al, .{});
    defer threaded.deinit();
    inline for (all) |T| {
        const C = Cases(T);
        // a dozen chunks or so of each
        inline for (.{ C.gemv, C.ger, C.gbmv, C.rank }, .{ "gemv", "ger", "gbmv", "rank updates" }, .{ 3000, 3000, 3000, 40 }) |body, name, most| {
            threading.test_serial_max = most;
            sameBits(al, threaded.io(), body) catch |e| {
                check.say("{s} {s}: ", .{ name, @typeName(T) });
                return e;
            };
        }
    }
}

test "threaded level 2 gives one worker's bits at any worker count, Io and core" {
    if (comptime !threading.hosted) return error.SkipZigTest;
    _ = blas.coreName();
    const cpus = blas.config.cpu_count;
    defer blas.config.cpu_count = cpus;
    defer threading.test_serial_max = null;
    try check.eachCore(onCore);
}
