//! PocketFFT's `fftblue` (`pocketfft_hdronly.h`, cpp branch, commit
//! 4d865cba, BSD-3, LICENSE): Bluestein's algorithm, a length-`n` DFT as a
//! convolution through `cfftp` of the 2, 3, 5, 7, 11-smooth `n2 >= 2n - 1`.
//! `exec`'s scratch is `2 n2` (PocketFFT allocates it per call).

const std = @import("std");
const cmplx = @import("cmplx.zig");
const twiddle = @import("twiddle.zig");
const util = @import("util.zig");
const cfftp = @import("cfftp/__root__.zig");

const Allocator = std.mem.Allocator;
const C = cmplx.C;

pub fn Bluestein(comptime T: type) type {
    return struct {
        n: usize,
        n2: usize,
        plan: cfftp.Cfftp(T),
        /// `bk` (`n`), then `bkf` (`n2 / 2 + 1`).
        mem: []C(T),

        const Self = @This();

        pub const InitError = Allocator.Error || error{TooLarge};

        pub fn init(al: Allocator, n: usize) InitError!Self {
            const n2 = util.goodSizeComplex(2 * n - 1) orelse return error.TooLarge;
            var plan = try cfftp.Cfftp(T).init(al, n2);
            errdefer plan.deinit(al);
            const mem = try al.alloc(C(T), n + n2 / 2 + 1);
            errdefer al.free(mem);
            const bk = mem[0..n];
            const bkf = mem[n..];
            // b_k
            var tw = try twiddle.Twiddles(T).init(al, 2 * n);
            defer tw.deinit(al);
            bk[0] = .{ .re = 1, .im = 0 };
            var coeff: usize = 0;
            for (1..n) |m| {
                coeff += 2 * m - 1;
                if (coeff >= 2 * n) coeff -= 2 * n;
                bk[m] = tw.at(coeff);
            }
            // the zero-padded, Fourier transformed b_k, with the normalization
            const tbkf = try al.alloc(C(T), 2 * n2);
            defer al.free(tbkf);
            const xn2 = 1 / @as(T, @floatFromInt(n2));
            tbkf[0] = cmplx.scale(T, bk[0], xn2);
            for (1..n) |m| {
                tbkf[m] = cmplx.scale(T, bk[m], xn2);
                tbkf[n2 - m] = tbkf[m];
            }
            for (n..n2 - n + 1) |m| tbkf[m] = .{ .re = 0, .im = 0 };
            plan.exec(tbkf[0..n2], tbkf[n2..], 1, true);
            @memcpy(bkf, tbkf[0 .. n2 / 2 + 1]);
            return .{ .n = n, .n2 = n2, .plan = plan, .mem = mem };
        }

        pub fn deinit(self: *Self, al: Allocator) void {
            self.plan.deinit(al);
            al.free(self.mem);
            self.* = undefined;
        }

        pub fn scratchLen(self: Self) usize {
            return 2 * self.n2;
        }

        pub fn exec(self: *const Self, c: []C(T), scratch: []C(T), fct: T, fwd: bool) void {
            if (fwd) self.fft(true, c, scratch, fct) else self.fft(false, c, scratch, fct);
        }

        /// Scratch for `execReal`: `n` complex values, then `exec`'s.
        pub fn realScratchLen(self: Self) usize {
            return self.n + self.scratchLen();
        }

        /// `exec_r`: `n` reals to the halfcomplex order (`fwd`) or back.
        pub fn execReal(self: *const Self, c: []T, scratch: []C(T), fct: T, fwd: bool) void {
            const n = self.n;
            const tmp = scratch[0..n];
            const rest = scratch[n..];
            const zero = 0 * c[0];
            if (fwd) {
                for (tmp, c[0..n]) |*t, x| t.* = .{ .re = x, .im = zero };
                self.fft(true, tmp, rest, fct);
                c[0] = tmp[0].re;
                // `tmp[1..]`'s parts, interleaved (the struct's layout is not)
                for (1..n) |j| c[j] = if (j & 1 == 1) tmp[(j + 1) / 2].re else tmp[j / 2].im;
            } else {
                tmp[0] = .{ .re = c[0], .im = c[0] * 0 };
                for (1..n) |j| {
                    if (j & 1 == 1) tmp[(j + 1) / 2].re = c[j] else tmp[j / 2].im = c[j];
                }
                if (n & 1 == 0) tmp[n / 2].im = 0 * c[0];
                var m: usize = 1;
                while (2 * m < n) : (m += 1) tmp[n - m] = .{ .re = tmp[m].re, .im = -tmp[m].im };
                self.fft(false, tmp, rest, fct);
                for (c[0..n], tmp) |*x, t| x.* = t.re;
            }
        }

        fn fft(self: *const Self, comptime fwd: bool, c: []C(T), scratch: []C(T), fct: T) void {
            const n = self.n;
            const n2 = self.n2;
            const bk = self.mem[0..n];
            const bkf = self.mem[n..];
            const akf = scratch[0..n2];
            const ch = scratch[n2 .. 2 * n2];
            // a_k, then its FFT
            for (0..n) |m| akf[m] = cmplx.twiddle(T, fwd, c[m], bk[m]);
            const zero = cmplx.scale(T, akf[0], 0);
            for (akf[n..]) |*a| a.* = zero;
            self.plan.exec(akf, ch, 1, true);
            // the convolution
            akf[0] = cmplx.twiddle(T, !fwd, akf[0], bkf[0]);
            for (1..(n2 + 1) / 2) |m| {
                akf[m] = cmplx.twiddle(T, !fwd, akf[m], bkf[m]);
                akf[n2 - m] = cmplx.twiddle(T, !fwd, akf[n2 - m], bkf[m]);
            }
            if (n2 & 1 == 0) akf[n2 / 2] = cmplx.twiddle(T, !fwd, akf[n2 / 2], bkf[n2 / 2]);
            // the inverse FFT, then times b_k
            self.plan.exec(akf, ch, 1, false);
            for (0..n) |m| c[m] = cmplx.scale(T, cmplx.twiddle(T, fwd, akf[m], bk[m]), fct);
        }
    };
}
