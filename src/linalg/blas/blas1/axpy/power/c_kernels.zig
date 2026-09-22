//! POWER axpy kernels that are C upstream: the drivers of
//! `kernel/power/{s,d,c,z}axpy.c` and their `_power10.c` forms around a block.
//!
//! The blocks are of two kinds. Where upstream is inline assembly (`d`, `z`
//! and the POWER10 ones) the block is the converted `uk/` file, called with
//! the scalar arguments by value as there. Where it is AltiVec intrinsics
//! (`saxpy.c`, `caxpy.c`) the block is the same loop on Zig vectors.

const std = @import("std");
const types = @import("../../../types.zig");
const fma = @import("../../../fma.zig");
const kernel = @import("../kernel.zig");
const uk = @import("../uk_block.zig");
const util = @import("../../../util.zig");

const BlasLong = types.BlasLong;
const Complex = std.math.Complex;

/// A real block that takes `alpha` by value: `f(n, x, y, alpha)`.
fn RealByValue(comptime T: type, comptime f: anytype) type {
    return struct {
        pub fn run(n1: BlasLong, alpha: T, x: [*]T, y: [*]T) void {
            f(n1, x, y, alpha);
        }
    };
}

/// A complex block: `f(n, x, y, alpha_r, alpha_i)` on the components.
fn CplxByValue(comptime F: type, comptime f: anytype) type {
    return struct {
        pub fn run(n1: BlasLong, alpha: Complex(F), x: [*]Complex(F), y: [*]Complex(F)) void {
            f(n1, @ptrCast(x), @ptrCast(y), alpha.re, alpha.im);
        }
    };
}

/// `saxpy_kernel_64` of `saxpy.c`: `y += alpha * x` in vectors of four, sixteen
/// at a time.
const SaxpyVec = struct {
    pub fn run(n1: BlasLong, alpha: f32, x: [*]f32, y: [*]f32) void {
        const V = @Vector(4, f32);
        const va: V = @splat(alpha);
        var i: usize = 0;
        const n: usize = @intCast(n1);
        while (i < n) : (i += 64) {
            inline for (0..16) |k| {
                const vx = util.load(V, x + i + 4 * k);
                const vy = util.load(V, y + i + 4 * k);
                util.store(V, y + i + 4 * k, if (fma.fuse) @mulAdd(V, vx, va, vy) else vx * va + vy);
            }
        }
    }
};

/// `caxpy_kernel_16` of `caxpy.c`: two multiply-adds per vector of two
/// complex values, the second on the real and imaginary parts swapped.
fn CaxpyVec(comptime conj: bool) type {
    return struct {
        pub fn run(n1: BlasLong, alpha: Complex(f32), x: [*]Complex(f32), y: [*]Complex(f32)) void {
            const V = @Vector(4, f32);
            const ar = alpha.re;
            const ai = alpha.im;
            const valpha_r: V = if (conj) .{ ar, -ar, ar, -ar } else .{ ar, ar, ar, ar };
            const valpha_i: V = if (conj) .{ ai, ai, ai, ai } else .{ -ai, ai, -ai, ai };
            const xf: [*]const f32 = @ptrCast(x);
            const yf: [*]f32 = @ptrCast(y);
            var i: usize = 0;
            const n: usize = @intCast(n1);
            while (i < 2 * n) : (i += 32) {
                inline for (0..8) |k| {
                    const vx = util.load(V, xf + i + 4 * k);
                    var vy = util.load(V, yf + i + 4 * k);
                    vy = if (fma.fuse) @mulAdd(V, vx, valpha_r, vy) else vx * valpha_r + vy;
                    const swapped = @shuffle(f32, vx, undefined, [4]i32{ 1, 0, 3, 2 });
                    vy = if (fma.fuse) @mulAdd(V, swapped, valpha_i, vy) else swapped * valpha_i + vy;
                    util.store(V, yf + i + 4 * k, vy);
                }
            }
        }
    };
}

/// The real drivers of `daxpy_power10.c` and `saxpy_power10.c`: as the plain
/// driver, but the contiguous case first steps to a 32-byte aligned `x`.
fn RealAligned(comptime T: type, comptime mask: comptime_int, comptime Block: type) type {
    const shift = if (T == f64) 3 else 2;
    const Strided = kernel.Real(T, mask, null, fma.fuse);
    return struct {
        pub fn call(n: BlasLong, alpha: T, x: [*]T, incx: BlasLong, y: [*]T, incy: BlasLong) void {
            if (n <= 0) return;
            if (incx != 1 or incy != 1) return Strided.call(n, alpha, x, incx, y, incy);

            var i: BlasLong = 0;
            if (n >= mask) {
                const misalign: BlasLong = @intCast(@intFromPtr(x) & 0x1F);
                const head = ((32 - misalign) >> shift) & (32 / @sizeOf(T) - 1);
                while (i < head) : (i += 1) {
                    const idx: usize = @intCast(i);
                    y[idx] = kernel.fmadd(fma.fuse, T, alpha, x[idx], y[idx]);
                }
            }
            const n1 = (n - i) & -mask;
            const at: usize = @intCast(i);
            if (n1 != 0) Block.run(n1, alpha, x + at, y + at);
            i += n1;
            while (i < n) : (i += 1) {
                const idx: usize = @intCast(i);
                y[idx] = kernel.fmadd(fma.fuse, T, alpha, x[idx], y[idx]);
            }
        }
    };
}

/// Where each role is, by core: which block, and whether it is one of the
/// aligned power10 drivers.
pub const Source = enum { power8, power10 };

/// The block module of `role` for `source`, by the file names of `uk/`.
fn ukModule(comptime role: []const u8, comptime source: Source) type {
    return switch (source) {
        .power8 => if (role[0] == 'd')
            @import("uk/daxpy_power8.zig")
        else if (role[0] == 'z')
            @import("uk/zaxpy_power8.zig")
        else
            struct {},
        .power10 => if (role[0] == 'd')
            @import("uk/daxpy_power10.zig")
        else if (role[0] == 's')
            @import("uk/saxpy_power10.zig")
        else if (role[0] == 'c')
            @import("uk/caxpy_power10.zig")
        else
            @import("uk/zaxpy_power10.zig"),
    };
}

/// The kernel `role` has on the cores of `source`, in the reference's ABI.
pub fn entry(comptime role: []const u8, comptime source: Source) uk.Kernel(role) {
    const conj = role.len == 6;
    const mod = ukModule(role, source);
    const fuse = fma.fuse;
    switch (role[0]) {
        'd' => {
            const Block = RealByValue(f64, mod.daxpy_kernel_8);
            if (source == .power10) return kernel.realEntry(f64, RealAligned(f64, 16, Block));
            return kernel.realEntry(f64, kernel.Real(f64, 16, Block, fuse));
        },
        's' => {
            if (source == .power10) {
                const Block = RealByValue(f32, mod.saxpy_kernel_64);
                return kernel.realEntry(f32, RealAligned(f32, 64, Block));
            }
            return kernel.realEntry(f32, kernel.Real(f32, 64, SaxpyVec, fuse));
        },
        'c' => {
            if (source == .power10) {
                const f = if (conj) mod.caxpyc_kernel_8 else mod.caxpy_kernel_8;
                return kernel.complexEntry(f32, kernel.Cplx(f32, conj, 16, CplxByValue(f32, f), fuse));
            }
            return kernel.complexEntry(f32, kernel.Cplx(f32, conj, 16, CaxpyVec(conj), fuse));
        },
        else => {
            const f = if (conj) mod.zaxpyc_kernel_4 else mod.zaxpy_kernel_4;
            return kernel.complexEntry(f64, kernel.Cplx(f64, conj, 16, CplxByValue(f64, f), fuse));
        },
    }
}
