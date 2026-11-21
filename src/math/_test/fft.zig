//! `rstd.math.fft` against a direct DFT, round trips, and its errors.

const std = @import("std");
const fft = @import("../fft.zig");

const testing = std.testing;
const Complex = std.math.Complex;

fn tol(comptime T: type) T {
    return if (T == f32) 1e-4 else 1e-11;
}

/// `x`'s DFT in a plain double loop.
fn dft(comptime T: type, out: []Complex(T), x: []const Complex(T), sign: T) void {
    const n: T = @floatFromInt(x.len);
    for (out, 0..) |*o, k| {
        var s: Complex(T) = .{ .re = 0, .im = 0 };
        for (x, 0..) |v, j| {
            const a = sign * 2 * std.math.pi * @as(T, @floatFromInt((j * k) % x.len)) / n;
            s = s.add(v.mul(.{ .re = @cos(a), .im = @sin(a) }));
        }
        o.* = s;
    }
}

fn fill(comptime T: type, x: []Complex(T), seed: u64) void {
    var prng: std.Random.DefaultPrng = .init(seed);
    const r = prng.random();
    for (x) |*v| v.* = .{ .re = r.float(T) * 2 - 1, .im = r.float(T) * 2 - 1 };
}

fn expectClose(comptime T: type, want: []const Complex(T), got: []const Complex(T)) !void {
    try testing.expectEqual(want.len, got.len);
    for (want, got) |w, g| {
        try testing.expectApproxEqAbs(w.re, g.re, tol(T));
        try testing.expectApproxEqAbs(w.im, g.im, tol(T));
    }
}

test "fft and ifft: the DFT, and back" {
    const al = testing.allocator;
    inline for (.{ f32, f64 }) |T| {
        for ([_]usize{ 1, 2, 7, 16, 60, 97 }) |n| {
            const x = try al.alloc(Complex(T), n);
            defer al.free(x);
            fill(T, x, n);
            const want = try al.alloc(Complex(T), n);
            defer al.free(want);
            dft(T, want, x, -1);
            const y = try fft.fft(al, testing.io, T, x);
            defer al.free(y);
            try expectClose(T, want, y);
            const back = try fft.ifft(al, testing.io, T, y);
            defer al.free(back);
            try expectClose(T, x, back);
            // ortho: both sides scaled by 1/√n, still a round trip
            const yo = try fft.fftWith(al, testing.io, T, x, .{ .normalization = .ortho });
            defer al.free(yo);
            const s = 1 / @sqrt(@as(T, @floatFromInt(n)));
            for (want, yo) |w, g| {
                try testing.expectApproxEqAbs(w.re * s, g.re, tol(T));
                try testing.expectApproxEqAbs(w.im * s, g.im, tol(T));
            }
            const bo = try fft.ifftWith(al, testing.io, T, yo, .{ .normalization = .ortho });
            defer al.free(bo);
            try expectClose(T, x, bo);
        }
    }
}

test "rfft and irfft: the half spectrum, and back" {
    const al = testing.allocator;
    inline for (.{ f32, f64 }) |T| {
        for ([_]usize{ 1, 2, 9, 32, 50, 101 }) |n| {
            const cx = try al.alloc(Complex(T), n);
            defer al.free(cx);
            fill(T, cx, n + 3);
            const x = try al.alloc(T, n);
            defer al.free(x);
            for (x, cx) |*r, *c| {
                r.* = c.re;
                c.im = 0;
            }
            const full = try al.alloc(Complex(T), n);
            defer al.free(full);
            dft(T, full, cx, -1);
            const half = try fft.rfft(al, testing.io, T, x);
            defer al.free(half);
            try expectClose(T, full[0 .. n / 2 + 1], half);
            const back = try fft.irfft(al, testing.io, T, n, half);
            defer al.free(back);
            for (x, back) |w, g| try testing.expectApproxEqAbs(w, g, tol(T));
            // forward: the forward side divides by n, the inverse not
            const hf = try fft.rfftWith(al, testing.io, T, x, .{ .normalization = .forward });
            defer al.free(hf);
            const bf = try fft.irfftWith(al, testing.io, T, n, hf, .{ .normalization = .forward });
            defer al.free(bf);
            for (x, bf) |w, g| try testing.expectApproxEqAbs(w, g, tol(T));
        }
    }
}

test "errors: empty input, a half spectrum of the wrong length" {
    const al = testing.allocator;
    try testing.expectError(error.ZeroLength, fft.fft(al, testing.io, f64, &.{}));
    try testing.expectError(error.ZeroLength, fft.rfft(al, testing.io, f64, &.{}));
    try testing.expectError(error.ZeroLength, fft.irfft(al, testing.io, f64, 0, &.{}));
    const half: [3]Complex(f64) = @splat(.{ .re = 1, .im = 0 });
    try testing.expectError(error.LengthMismatch, fft.irfft(al, testing.io, f64, 6, &half));
    const back = try fft.irfft(al, testing.io, f64, 5, &half);
    defer al.free(back);
    try testing.expectEqual(@as(usize, 5), back.len);
}
