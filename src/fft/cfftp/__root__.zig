//! PocketFFT's `cfftp` (`pocketfft_hdronly.h`, cpp branch, commit 4d865cba,
//! BSD-3, LICENSE): the complex mixed-radix FFT (FFTPACK's): the length's
//! factors, their twiddles, and a pass per factor. `exec` takes its
//! ping-pong buffer from the caller, where PocketFFT allocates one.

const std = @import("std");
const mem = @import("../../mem.zig");
const cmplx = @import("../cmplx.zig");
const twiddle = @import("../twiddle.zig");
const passes = @import("passes.zig");
const passg_ = @import("passg.zig");

const Allocator = std.mem.Allocator;
const C = cmplx.C;

pub fn Cfftp(comptime T: type) type {
    return struct {
        length: usize,
        nfact: usize = 0,
        fact: [max_factors]Factor = undefined,
        /// Every factor's twiddles, `tw` and (past 11) `tws`.
        mem: []C(T) = &.{},

        const Self = @This();
        const P = passes.Passes(T);
        /// A length's factors past the first are at least 2.
        const max_factors = @bitSizeOf(usize);
        const Factor = struct { fct: usize, tw: usize = 0, tws: usize = 0 };

        pub fn init(al: Allocator, length: usize) Allocator.Error!Self {
            var self: Self = .{ .length = length };
            if (length <= 1) return self;
            self.factorize();
            self.mem = try al.alloc(C(T), self.twsize());
            errdefer al.free(self.mem);
            try self.compTwiddle(al);
            return self;
        }

        pub fn deinit(self: *Self, al: Allocator) void {
            al.free(self.mem);
            self.* = undefined;
        }

        /// `c := fct * FFT(c)`, `ch` (at least `length`) the ping-pong buffer.
        pub fn exec(self: *const Self, c: []C(T), ch: []C(T), fct: T, fwd: bool) void {
            if (fwd) self.passAll(true, c, ch, fct) else self.passAll(false, c, ch, fct);
        }

        fn passAll(self: *const Self, comptime fwd: bool, c: []C(T), ch: []C(T), fct: T) void {
            const len = self.length;
            if (len == 1) {
                c[0] = cmplx.scale(T, c[0], fct);
                return;
            }
            var l1: usize = 1;
            var p1: [*]C(T) = c.ptr;
            var p2: [*]C(T) = ch.ptr;
            for (self.fact[0..self.nfact]) |f| {
                const ip = f.fct;
                const l2 = ip * l1;
                const ido = len / l2;
                const tw = self.mem.ptr + f.tw;
                switch (ip) {
                    4 => P.pass4(fwd, ido, l1, p1, p2, tw),
                    8 => P.pass8(fwd, ido, l1, p1, p2, tw),
                    2 => P.pass2(fwd, ido, l1, p1, p2, tw),
                    3 => P.passOdd(3, fwd, ido, l1, p1, p2, tw),
                    5 => P.passOdd(5, fwd, ido, l1, p1, p2, tw),
                    7 => P.passOdd(7, fwd, ido, l1, p1, p2, tw),
                    11 => P.passOdd(11, fwd, ido, l1, p1, p2, tw),
                    else => {
                        passg_.passg(T, fwd, ido, ip, l1, p1, p2, tw, self.mem.ptr + f.tws);
                        mem.swap([*]C(T), &p1, &p2);
                    },
                }
                mem.swap([*]C(T), &p1, &p2);
                l1 = l2;
            }
            if (p1 != c.ptr) {
                if (fct != 1) {
                    for (c, ch[0..len]) |*x, y| x.* = cmplx.scale(T, y, fct);
                } else @memcpy(c, p1[0..len]);
            } else if (fct != 1) {
                for (c) |*x| x.* = cmplx.scale(T, x.*, fct);
            }
        }

        fn addFactor(self: *Self, f: usize) void {
            self.fact[self.nfact] = .{ .fct = f };
            self.nfact += 1;
        }

        fn factorize(self: *Self) void {
            var len = self.length;
            while (len & 7 == 0) {
                self.addFactor(8);
                len >>= 3;
            }
            while (len & 3 == 0) {
                self.addFactor(4);
                len >>= 2;
            }
            if (len & 1 == 0) {
                len >>= 1;
                // factor 2 should be at the front of the factor list
                self.addFactor(2);
                mem.swap(usize, &self.fact[0].fct, &self.fact[self.nfact - 1].fct);
            }
            var divisor: usize = 3;
            while (divisor * divisor <= len) : (divisor += 2) {
                while (len % divisor == 0) {
                    self.addFactor(divisor);
                    len /= divisor;
                }
            }
            if (len > 1) self.addFactor(len);
        }

        fn twsize(self: *const Self) usize {
            var size: usize = 0;
            var l1: usize = 1;
            for (self.fact[0..self.nfact]) |f| {
                const ido = self.length / (l1 * f.fct);
                size += (f.fct - 1) * (ido - 1);
                if (f.fct > 11) size += f.fct;
                l1 *= f.fct;
            }
            return size;
        }

        fn compTwiddle(self: *Self, al: Allocator) Allocator.Error!void {
            var tw = try twiddle.Twiddles(T).init(al, self.length);
            defer tw.deinit(al);
            var l1: usize = 1;
            var ofs: usize = 0;
            for (self.fact[0..self.nfact]) |*f| {
                const ip = f.fct;
                const ido = self.length / (l1 * ip);
                f.tw = ofs;
                ofs += (ip - 1) * (ido - 1);
                for (1..ip) |j| for (1..ido) |i| {
                    self.mem[f.tw + (j - 1) * (ido - 1) + i - 1] = tw.at(j * l1 * i);
                };
                if (ip > 11) {
                    f.tws = ofs;
                    ofs += ip;
                    for (0..ip) |j| self.mem[f.tws + j] = tw.at(j * l1 * ido);
                }
                l1 *= ip;
            }
        }
    };
}
