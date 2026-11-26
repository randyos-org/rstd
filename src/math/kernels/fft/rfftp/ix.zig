//! `rfftp`'s accessors and small operations from PocketFFT
//! (`pocketfft_hdronly.h`, cpp branch, commit 4d865cba, BSD-3, LICENSE),
//! named as its lambdas and macros.

/// The forward passes': `CC` by `l1`, `CH` by the radix.
pub fn Forward(comptime T: type) type {
    return Access(T, true);
}

/// The backward passes': `CC` by the radix, `CH` by `l1`.
pub fn Backward(comptime T: type) type {
    return Access(T, false);
}

fn Access(comptime T: type, comptime fwd: bool) type {
    return struct {
        ido: usize,
        l1: usize,
        cdim: usize,
        cc: [*]const T,
        ch: [*]T,
        wa: [*]const T,

        const Self = @This();

        pub inline fn CC(s: Self, a: usize, b: usize, c: usize) T {
            return s.cc[a + s.ido * (b + (if (fwd) s.l1 else s.cdim) * c)];
        }
        pub inline fn CH(s: Self, a: usize, b: usize, c: usize) *T {
            return &s.ch[a + s.ido * (b + (if (fwd) s.cdim else s.l1) * c)];
        }
        pub inline fn WA(s: Self, x: usize, i: usize) T {
            return s.wa[i + x * (s.ido - 1)];
        }
    };
}

/// `PM`: `a = c + d`, `b = c - d`.
pub inline fn pm(comptime T: type, a: *T, b: *T, c: T, d: T) void {
    a.* = c + d;
    b.* = c - d;
}

/// `MULPM`: `(a + ib) = conj(c + id) (e + if)`.
pub inline fn mulpm(comptime T: type, c: T, d: T, e: T, f: T) struct { a: T, b: T } {
    return .{ .a = c * e + d * f, .b = c * f - d * e };
}

/// `POCKETFFT_REARRANGE`: `a2 = a + b`, `b2 = i (b - a)`.
pub inline fn rearrange(comptime T: type, rx: T, ix_: T, ry: T, iy: T) struct { rx: T, ix: T, ry: T, iy: T } {
    return .{ .rx = rx + ry, .ix = ix_ + iy, .ry = ix_ - iy, .iy = ry - rx };
}
