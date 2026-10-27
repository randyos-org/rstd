//! Shared pieces of the dot kernels: the reference's kernel ABI and the
//! drivers of `kernel/x86_64/{s,d,c,z}dot.c`. Nothing here depends on the BLAS
//! root, so the per-core kernel objects (`_kernel.zig`) can use it too.
//!
//! A driver is a block kernel for the contiguous part, the reference's scalar
//! tail, and its strided loops. A block leaves its sum in the memory the
//! driver hands it; here `Block.run(n1, x, y)` returns what it left there
//! (zero before). `fuse` says whether the compiler flags of the core being
//! built contract multiply-adds (see `fma.zig`); it is applied per statement
//! the way clang's contraction is: only within one expression, so a product
//! kept in its own statement stays un-fused, and a `float` product added to a
//! `double` is rounded first.
//!
//! The reference's `SMP` path, which splits a long vector over threads and
//! adds the partial sums, is the interface's concern, not the kernel's: these
//! are its single-threaded `dot_compute`.

const std = @import("std");
const util = @import("../../util.zig");
const types = @import("../../types.zig");
const fma = @import("../../fma.zig");

const BlasLong = types.BlasLong;
const CComplex = types.CComplex;
const Slot = util.Slot;

/// The kernels an architecture's `<arch>/kernel.zig` iterates over.
pub const slots = [_]Slot{ .sdot_k, .dsdot_k, .ddot_k, .cdotu_k, .cdotc_k, .zdotu_k, .zdotc_k };

/// What a kernel works on: `E` is its element type (a component for the
/// complex ones), `R` what it returns.
pub const Form = struct { E: type, R: type, complex: bool };

/// The form of `slot`, one of `slots`.
pub fn form(comptime slot: Slot) Form {
    return switch (slot) {
        .sdot_k => .{ .E = f32, .R = f32, .complex = false },
        .dsdot_k => .{ .E = f32, .R = f64, .complex = false },
        .ddot_k => .{ .E = f64, .R = f64, .complex = false },
        .cdotu_k, .cdotc_k => .{ .E = f32, .R = CComplex(f32), .complex = true },
        .zdotu_k, .zdotc_k => .{ .E = f64, .R = CComplex(f64), .complex = true },
        else => @compileError(@tagName(slot) ++ " is not a dot kernel"),
    };
}

// ---- kernel ABI ------------------------------------------------------------

/// Vtable slot type of `slot`; the elements are `form(slot).E`s, interleaved
/// pairs for the complex ones.
pub fn Kernel(comptime slot: Slot) type {
    const E = form(slot).E;
    return *const fn (BlasLong, [*]E, BlasLong, [*]E, BlasLong) callconv(.c) form(slot).R;
}

/// How many partial sums the complex kernels' blocks leave: `cdot` keeps two
/// sets of four, `zdot` one.
fn partials(comptime slot: Slot) comptime_int {
    return if (form(slot).E == f32) 8 else 4;
}

// ---- blocks ----------------------------------------------------------------

/// A block calling the reference's `f(n, x, y, d)`: for a real slot `d` is one
/// element, for a complex one the array of partial sums.
pub fn Asm(comptime slot: Slot, comptime f: anytype) type {
    const E = form(slot).E;
    if (form(slot).complex) return struct {
        pub fn run(n1: BlasLong, x: [*]E, y: [*]E) [partials(slot)]E {
            var d: [partials(slot)]E = @splat(0);
            f(n1, x, y, &d);
            return d;
        }
    };
    return struct {
        pub fn run(n1: BlasLong, x: [*]E, y: [*]E) E {
            var d: E = 0;
            f(n1, x, y, @ptrCast(&d));
            return d;
        }
    };
}

/// The block the reference's driver falls back on where a core has no
/// microkernel (`ddot_kernel_8`, `sdot_kernel_16`, `cdot_kernel_16`,
/// `zdot_kernel_8`): plain C, in its own summation order.
pub fn Fallback(comptime slot: Slot, comptime fuse: bool) type {
    const E = form(slot).E;
    if (form(slot).complex) {
        const N = partials(slot);
        // the slots of one complex element's four products (re*re, im*im,
        // re*im, im*re); `cdot` alternates two sets of them
        const places: [2][4]usize = if (N == 8) .{ .{ 0, 1, 4, 5 }, .{ 2, 3, 6, 7 } } else .{ .{ 0, 1, 2, 3 }, .{ 0, 1, 2, 3 } };
        return struct {
            pub fn run(n1: BlasLong, x: [*]E, y: [*]E) [N]E {
                var dot: [N]E = @splat(0);
                var j: usize = 0;
                var i: BlasLong = 0;
                while (i < n1) : (i += 4) {
                    inline for (0..4) |k| {
                        const s = places[k % 2];
                        const xr = x[j + 2 * k];
                        const xi = x[j + 2 * k + 1];
                        const yr = y[j + 2 * k];
                        const yi = y[j + 2 * k + 1];
                        dot[s[0]] = fma.madd(fuse, E, xr, yr, dot[s[0]]);
                        dot[s[1]] = fma.madd(fuse, E, xi, yi, dot[s[1]]);
                        dot[s[2]] = fma.madd(fuse, E, xr, yi, dot[s[2]]);
                        dot[s[3]] = fma.madd(fuse, E, xi, yr, dot[s[3]]);
                    }
                    j += 8;
                }
                return dot;
            }
        };
    }
    return struct {
        pub fn run(n1: BlasLong, x: [*]E, y: [*]E) E {
            var dot: E = 0;
            var i: usize = 0;
            const n: usize = @intCast(n1);
            while (i < n) : (i += 8) {
                // `y0*x0 + y1*x1 + ...`: clang fuses the first product, then
                // each next one into the running sum
                var r = fma.madd(fuse, E, y[i], x[i], y[i + 1] * x[i + 1]);
                inline for (2..8) |k| r = fma.madd(fuse, E, y[i + k], x[i + k], r);
                dot += r;
            }
            return dot;
        }
    };
}

// ---- drivers ---------------------------------------------------------------

/// `ddot.c`.
fn Ddot(comptime mask: comptime_int, comptime Block: type, comptime fuse: bool) type {
    return struct {
        pub fn call(n: BlasLong, x: [*]f64, incx: BlasLong, y: [*]f64, incy: BlasLong) f64 {
            if (n <= 0) return 0;

            if (incx == 1 and incy == 1) {
                var dot: f64 = 0;
                const n1 = n & -mask;
                if (n1 != 0) dot = Block.run(n1, x, y);
                var i = n1;
                while (i < n) : (i += 1) {
                    const idx: usize = @intCast(i);
                    dot = fma.madd(fuse, f64, y[idx], x[idx], dot);
                }
                return dot;
            }

            var temp1: f64 = 0;
            var temp2: f64 = 0;
            const n1 = n & -4;
            var i: BlasLong = 0;
            var ix: BlasLong = 0;
            var iy: BlasLong = 0;
            while (i < n1) : (i += 4) {
                const m1 = util.at(f64, y, iy).* * util.at(f64, x, ix).*;
                const m2 = util.at(f64, y, iy + incy).* * util.at(f64, x, ix + incx).*;
                const m3 = util.at(f64, y, iy + 2 * incy).* * util.at(f64, x, ix + 2 * incx).*;
                const m4 = util.at(f64, y, iy + 3 * incy).* * util.at(f64, x, ix + 3 * incx).*;
                ix += incx * 4;
                iy += incy * 4;
                temp1 += m1 + m3;
                temp2 += m2 + m4;
            }
            while (i < n) : (i += 1) {
                temp1 = fma.madd(fuse, f64, util.at(f64, y, iy).*, util.at(f64, x, ix).*, temp1);
                ix += incx;
                iy += incy;
            }
            return temp1 + temp2;
        }
    };
}

/// `sdot.c`: the sum is a `double`, the products `float`.
fn Sdot(comptime mask: comptime_int, comptime Block: type, comptime fuse: bool) type {
    return struct {
        pub fn call(n: BlasLong, x: [*]f32, incx: BlasLong, y: [*]f32, incy: BlasLong) f32 {
            if (n <= 0) return 0;
            var dot: f64 = 0;

            if (incx == 1 and incy == 1) {
                var mydot: f32 = 0;
                const n1 = n & -mask;
                if (n1 != 0) mydot = Block.run(n1, x, y);
                var i = n1;
                while (i < n) : (i += 1) {
                    const idx: usize = @intCast(i);
                    const m = y[idx] * x[idx];
                    dot += m;
                }
                dot += mydot;
                return @floatCast(dot);
            }

            const n1 = n & -2;
            var i: BlasLong = 0;
            var ix: BlasLong = 0;
            var iy: BlasLong = 0;
            while (i < n1) : (i += 2) {
                const s = fma.madd(fuse, f32, util.at(f32, y, iy).*, util.at(f32, x, ix).*, util.at(f32, y, iy + incy).* * util.at(f32, x, ix + incx).*);
                dot += s;
                ix += incx * 2;
                iy += incy * 2;
            }
            while (i < n) : (i += 1) {
                const m = util.at(f32, y, iy).* * util.at(f32, x, ix).*;
                dot += m;
                ix += incx;
                iy += incy;
            }
            return @floatCast(dot);
        }
    };
}

/// `sdot.c` built as `DSDOT`: `float` operands, `double` products and sum. Its
/// block gets chunks of 32 whatever the block size is.
fn Dsdot(comptime mask: comptime_int, comptime Block: type, comptime fuse: bool) type {
    return struct {
        pub fn call(n: BlasLong, x: [*]f32, incx: BlasLong, y: [*]f32, incy: BlasLong) f64 {
            if (n <= 0) return 0;
            var dot: f64 = 0;

            if (incx == 1 and incy == 1) {
                var mydot: f64 = 0;
                const n1 = n & -mask;
                var i: BlasLong = 0;
                while (i < n1) : (i += 32) {
                    const idx: usize = @intCast(i);
                    mydot += Block.run(32, x + idx, y + idx);
                }
                i = n1;
                while (i < n) : (i += 1) {
                    const idx: usize = @intCast(i);
                    dot = fma.madd(fuse, f64, y[idx], x[idx], dot);
                }
                dot += mydot;
                return dot;
            }

            const n1 = n & -2;
            var i: BlasLong = 0;
            var ix: BlasLong = 0;
            var iy: BlasLong = 0;
            while (i < n1) : (i += 2) {
                const s = fma.madd(fuse, f64, util.at(f32, y, iy).*, util.at(f32, x, ix).*, @as(f64, util.at(f32, y, iy + incy).*) * util.at(f32, x, ix + incx).*);
                dot += s;
                ix += incx * 2;
                iy += incy * 2;
            }
            while (i < n) : (i += 1) {
                dot = fma.madd(fuse, f64, util.at(f32, y, iy).*, util.at(f32, x, ix).*, dot);
                ix += incx;
                iy += incy;
            }
            return dot;
        }
    };
}

/// `cdot.c` (`F` is `f32`: two sets of four partial sums) and `zdot.c`
/// (`f64`: one); `conj` is `CONJ`, the `*dotc` form.
fn Cdot(comptime F: type, comptime conj: bool, comptime mask: comptime_int, comptime Block: type, comptime fuse: bool) type {
    const N = if (F == f32) 8 else 4;
    // re*re, im*im, re*im, im*re
    const s: [4]usize = if (F == f32) .{ 0, 1, 4, 5 } else .{ 0, 1, 2, 3 };
    return struct {
        inline fn acc(dot: *[N]F, xr: F, xi: F, yr: F, yi: F) void {
            dot[s[0]] = fma.madd(fuse, F, xr, yr, dot[s[0]]);
            dot[s[1]] = fma.madd(fuse, F, xi, yi, dot[s[1]]);
            dot[s[2]] = fma.madd(fuse, F, xr, yi, dot[s[2]]);
            dot[s[3]] = fma.madd(fuse, F, xi, yr, dot[s[3]]);
        }

        pub fn call(n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong) CComplex(F) {
            var dot: [N]F = @splat(0);
            if (n <= 0) return .{ .re = 0, .im = 0 };

            if (incx == 1 and incy == 1) {
                const n1 = n & -mask;
                if (n1 != 0) {
                    dot = Block.run(n1, x, y);
                    if (comptime N == 8) {
                        dot[0] += dot[2];
                        dot[1] += dot[3];
                        dot[4] += dot[6];
                        dot[5] += dot[7];
                    }
                }
                var i = n1;
                var j: usize = @intCast(i * 2);
                while (i < n) : (i += 1) {
                    acc(&dot, x[j], x[j + 1], y[j], y[j + 1]);
                    j += 2;
                }
            } else {
                var i: BlasLong = 0;
                var ix: BlasLong = 0;
                var iy: BlasLong = 0;
                while (i < n) : (i += 1) {
                    acc(&dot, util.at(F, x, ix).*, util.at(F, x, ix + 1).*, util.at(F, y, iy).*, util.at(F, y, iy + 1).*);
                    ix += incx * 2;
                    iy += incy * 2;
                }
            }

            return if (conj)
                .{ .re = dot[s[0]] + dot[s[1]], .im = dot[s[2]] - dot[s[3]] }
            else
                .{ .re = dot[s[0]] - dot[s[1]], .im = dot[s[2]] + dot[s[3]] };
        }
    };
}

/// The driver of `slot` as a struct with `call`, around `Block` (the block
/// the core's microkernel gives, or the reference's own fallback).
pub fn Driver(comptime slot: Slot, comptime mask: comptime_int, comptime Block: type, comptime fuse: bool) type {
    return switch (slot) {
        .sdot_k => Sdot(mask, Block, fuse),
        .dsdot_k => Dsdot(mask, Block, fuse),
        .ddot_k => Ddot(mask, Block, fuse),
        .cdotu_k => Cdot(f32, false, mask, Block, fuse),
        .cdotc_k => Cdot(f32, true, mask, Block, fuse),
        .zdotu_k => Cdot(f64, false, mask, Block, fuse),
        .zdotc_k => Cdot(f64, true, mask, Block, fuse),
        else => @compileError(@tagName(slot) ++ " is not a dot kernel"),
    };
}

// ---- the driver, as `uk_block.driven` takes it -----------------------------

/// The kernel `slot` around `Block` (null for the reference's fallback
/// block), blocks of `mask`.
pub fn entry(comptime slot: Slot, comptime mask: comptime_int, comptime Block: ?type) Kernel(slot) {
    const E = form(slot).E;
    const Impl = Driver(slot, mask, Block orelse Fallback(slot, fma.fuse), fma.fuse);
    return &struct {
        fn f(n: BlasLong, x: [*]E, incx: BlasLong, y: [*]E, incy: BlasLong) callconv(.c) form(slot).R {
            return Impl.call(n, x, incx, y, incy);
        }
    }.f;
}
