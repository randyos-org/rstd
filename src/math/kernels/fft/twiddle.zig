//! PocketFFT's `sincos_2pibyn` (`pocketfft_hdronly.h`, cpp branch, commit
//! 4d865cba, BSD-3, LICENSE): `exp(2πi k/n)` from two tables of `√n` entries
//! each, octant by octant, in f64 (wider for wider types).

const std = @import("std");
const cmplx = @import("cmplx.zig");

const Allocator = std.mem.Allocator;
const C = cmplx.C;

pub fn Twiddles(comptime T: type) type {
    return struct {
        n: usize,
        mask: usize,
        shift: std.math.Log2Int(usize),
        v1: []C(H),
        v2: []C(H),

        const Self = @This();
        /// `Thigh`: double, or `T` when wider.
        const H = if (@bitSizeOf(T) > 64) T else f64;

        pub fn init(al: Allocator, n: usize) Allocator.Error!Self {
            const ang: H = 0.25 * std.math.pi / @as(H, @floatFromInt(n));
            const nval = (n + 2) / 2;
            var shift: std.math.Log2Int(usize) = 1;
            while ((@as(usize, 1) << shift) * (@as(usize, 1) << shift) < nval) shift += 1;
            const mask = (@as(usize, 1) << shift) - 1;
            const v1 = try al.alloc(C(H), mask + 1);
            errdefer al.free(v1);
            v1[0] = .{ .re = 1, .im = 0 };
            for (v1[1..], 1..) |*v, i| v.* = calc(i, n, ang);
            const v2 = try al.alloc(C(H), (nval + mask) / (mask + 1));
            v2[0] = .{ .re = 1, .im = 0 };
            for (v2[1..], 1..) |*v, i| v.* = calc(i * (mask + 1), n, ang);
            return .{ .n = n, .mask = mask, .shift = shift, .v1 = v1, .v2 = v2 };
        }

        pub fn deinit(self: *Self, al: Allocator) void {
            al.free(self.v1);
            al.free(self.v2);
            self.* = undefined;
        }

        /// `exp(2πi idx/n)`, `idx < n`.
        pub fn at(self: Self, idx: usize) C(T) {
            const flip = 2 * idx > self.n;
            const j = if (flip) self.n - idx else idx;
            const x1 = self.v1[j & self.mask];
            const x2 = self.v2[j >> self.shift];
            const re: T = @floatCast(x1.re * x2.re - x1.im * x2.im);
            const im: T = @floatCast(x1.re * x2.im + x1.im * x2.re);
            return .{ .re = re, .im = if (flip) -im else im };
        }

        fn calc(x0: usize, n: usize, ang: H) C(H) {
            var x = x0 << 3;
            if (x < 4 * n) { // first half
                if (x < 2 * n) { // first quadrant
                    if (x < n) return .{ .re = @cos(f(x) * ang), .im = @sin(f(x) * ang) };
                    return .{ .re = @sin(f(2 * n - x) * ang), .im = @cos(f(2 * n - x) * ang) };
                }
                // second quadrant
                x -= 2 * n;
                if (x < n) return .{ .re = -@sin(f(x) * ang), .im = @cos(f(x) * ang) };
                return .{ .re = -@cos(f(2 * n - x) * ang), .im = @sin(f(2 * n - x) * ang) };
            }
            x = 8 * n - x;
            if (x < 2 * n) { // third quadrant
                if (x < n) return .{ .re = @cos(f(x) * ang), .im = -@sin(f(x) * ang) };
                return .{ .re = @sin(f(2 * n - x) * ang), .im = -@cos(f(2 * n - x) * ang) };
            }
            // fourth quadrant
            x -= 2 * n;
            if (x < n) return .{ .re = -@sin(f(x) * ang), .im = -@cos(f(x) * ang) };
            return .{ .re = -@cos(f(2 * n - x) * ang), .im = -@sin(f(2 * n - x) * ang) };
        }

        inline fn f(x: usize) H {
            return @floatFromInt(x);
        }
    };
}
