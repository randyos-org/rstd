//! PocketFFT's `rfftp` (`pocketfft_hdronly.h`, cpp branch, commit 4d865cba,
//! BSD-3, LICENSE): FFTPACK's real FFT, `n` reals to the halfcomplex order
//! `r0, r1, i1, r2, i2, ...` (`r(n/2)` last for even `n`) and back. `exec`
//! takes its ping-pong buffer from the caller, where PocketFFT allocates one.

const std = @import("std");
const mem = @import("../../../../mem.zig");
const twiddle = @import("../twiddle.zig");
const forward = @import("forward.zig");
const backward = @import("backward.zig");
const general = @import("general.zig");

const Allocator = std.mem.Allocator;

pub fn Rfftp(comptime T: type) type {
    return struct {
        length: usize,
        nfact: usize = 0,
        fact: [max_factors]Factor = undefined,
        /// Every factor's twiddles, `tw` and (past 5) `tws`.
        mem: []T = &.{},

        const Self = @This();
        const F = forward.Forward(T);
        const B = backward.Backward(T);
        const max_factors = @bitSizeOf(usize);
        const Factor = struct { fct: usize, tw: usize = 0, tws: usize = 0 };

        pub fn init(al: Allocator, length: usize) Allocator.Error!Self {
            var self: Self = .{ .length = length };
            if (length <= 1) return self;
            self.factorize();
            self.mem = try al.alloc(T, self.twsize());
            errdefer al.free(self.mem);
            try self.compTwiddle(al);
            return self;
        }

        pub fn deinit(self: *Self, al: Allocator) void {
            al.free(self.mem);
            self.* = undefined;
        }

        /// `c := fct * FFT(c)`: real to halfcomplex (`r2hc`) or back; `ch` (at
        /// least `length`) the ping-pong buffer.
        pub fn exec(self: *const Self, c: []T, ch: []T, fct: T, r2hc: bool) void {
            const len = self.length;
            if (len == 1) {
                c[0] *= fct;
                return;
            }
            const nf = self.nfact;
            var p1: [*]T = c.ptr;
            var p2: [*]T = ch.ptr;
            if (r2hc) {
                var l1 = len;
                for (0..nf) |k1| {
                    const f = self.fact[nf - k1 - 1];
                    const ip = f.fct;
                    const ido = len / l1;
                    l1 /= ip;
                    const tw = self.mem.ptr + f.tw;
                    switch (ip) {
                        4 => F.radf4(ido, l1, p1, p2, tw),
                        2 => F.radf2(ido, l1, p1, p2, tw),
                        3 => F.radf3(ido, l1, p1, p2, tw),
                        5 => F.radf5(ido, l1, p1, p2, tw),
                        else => {
                            general.radfg(T, ido, ip, l1, p1, p2, tw, self.mem.ptr + f.tws);
                            mem.swap([*]T, &p1, &p2);
                        },
                    }
                    mem.swap([*]T, &p1, &p2);
                }
            } else {
                var l1: usize = 1;
                for (self.fact[0..nf]) |f| {
                    const ip = f.fct;
                    const ido = len / (ip * l1);
                    const tw = self.mem.ptr + f.tw;
                    switch (ip) {
                        4 => B.radb4(ido, l1, p1, p2, tw),
                        2 => B.radb2(ido, l1, p1, p2, tw),
                        3 => B.radb3(ido, l1, p1, p2, tw),
                        5 => B.radb5(ido, l1, p1, p2, tw),
                        else => general.radbg(T, ido, ip, l1, p1, p2, tw, self.mem.ptr + f.tws),
                    }
                    mem.swap([*]T, &p1, &p2);
                    l1 *= ip;
                }
            }
            // `copy_and_norm`
            if (p1 != c.ptr) {
                if (fct != 1) {
                    for (c, p1[0..len]) |*x, y| x.* = fct * y;
                } else @memcpy(c, p1[0..len]);
            } else if (fct != 1) {
                for (c) |*x| x.* *= fct;
            }
        }

        fn addFactor(self: *Self, f: usize) void {
            self.fact[self.nfact] = .{ .fct = f };
            self.nfact += 1;
        }

        fn factorize(self: *Self) void {
            var len = self.length;
            while (len % 4 == 0) {
                self.addFactor(4);
                len >>= 2;
            }
            if (len % 2 == 0) {
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
                if (f.fct > 5) size += 2 * f.fct;
                l1 *= f.fct;
            }
            return size;
        }

        fn compTwiddle(self: *Self, al: Allocator) Allocator.Error!void {
            var tw = try twiddle.Twiddles(T).init(al, self.length);
            defer tw.deinit(al);
            var l1: usize = 1;
            var ofs: usize = 0;
            for (self.fact[0..self.nfact], 0..) |*f, k| {
                const ip = f.fct;
                const ido = self.length / (l1 * ip);
                // the last factor doesn't need twiddles
                if (k < self.nfact - 1) {
                    f.tw = ofs;
                    ofs += (ip - 1) * (ido - 1);
                    for (1..ip) |j| for (1..(ido - 1) / 2 + 1) |i| {
                        const w = tw.at(j * l1 * i);
                        self.mem[f.tw + (j - 1) * (ido - 1) + 2 * i - 2] = w.re;
                        self.mem[f.tw + (j - 1) * (ido - 1) + 2 * i - 1] = w.im;
                    };
                }
                // special factors required by the general passes
                if (ip > 5) {
                    f.tws = ofs;
                    ofs += 2 * ip;
                    const s = self.mem[f.tws..][0 .. 2 * ip];
                    s[0] = 1;
                    s[1] = 0;
                    var i: usize = 2;
                    var ic: usize = 2 * ip - 2;
                    while (i <= ic) : ({
                        i += 2;
                        ic -= 2;
                    }) {
                        const w = tw.at(i / 2 * (self.length / ip));
                        s[i] = w.re;
                        s[i + 1] = w.im;
                        s[ic] = w.re;
                        s[ic + 1] = -w.im;
                    }
                }
                l1 *= ip;
            }
        }
    };
}
