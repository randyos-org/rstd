//! The drivers of `kernel/x86_64/{s,d,c,z}scal.c` around a core's `uk/`
//! blocks; every element is one product or a stored zero, so any block gives
//! the same bits. `flag` (upstream's `dummy2`): 1 multiplies even by a zero `alpha`.
//! The complex C tails, which branch on special values, take the blocks'
//! full product here (issue #231). No BLAS-root imports (`_kernel.zig`).

const std = @import("std");
const util = @import("../../util.zig");
const types = @import("../../types.zig");
const uk = @import("../../uk_block.zig");
const coretypes = @import("../../cpu/coretypes.zig");

const BlasLong = types.BlasLong;
const Slot = util.Slot;

/// The kernels an architecture's `<arch>/kernel.zig` iterates over.
pub const slots = [_]Slot{ .sscal_k, .dscal_k, .cscal_k, .zscal_k };

/// The float type of `slot`, one of `slots` (a component for the complex ones).
pub fn Elem(comptime slot: Slot) type {
    return switch (slot) {
        .sscal_k, .cscal_k => f32,
        .dscal_k, .zscal_k => f64,
        else => @compileError(@tagName(slot) ++ " is not a scal kernel"),
    };
}

/// Whether `slot` is a complex one.
pub fn complex(comptime slot: Slot) bool {
    return slot == .cscal_k or slot == .zscal_k;
}

/// Vtable slot type of `slot`.
pub fn Kernel(comptime slot: Slot) type {
    return util.SlotFn(slot);
}

/// The unit-stride block size: `sscal_kernel_16`, `dscal_kernel_8`, in elements.
pub fn width(comptime slot: Slot) comptime_int {
    return if (slot == .sscal_k or slot == .cscal_k) 16 else 8;
}

/// One complex element times `ar + ai*i`, unfused, as the blocks compute it.
pub inline fn product(comptime F: type, ar: F, ai: F, re: *F, im: *F) void {
    const xr = re.*;
    const xi = im.*;
    re.* = ar * xr - ai * xi;
    im.* = ar * xi + ai * xr;
}

/// The C blocks of `cscal.c` and `zscal.c` (`#if !defined(HAVE_KERNEL_8)`).
fn ComplexFallback(comptime F: type) type {
    return struct {
        pub fn run(n: BlasLong, alpha: [*]const F, x: [*]F) void {
            for (0..@intCast(n)) |i| product(F, alpha[0], alpha[1], &x[2 * i], &x[2 * i + 1]);
        }
        pub fn zero(n: BlasLong, _: [*]const F, x: [*]F) void {
            for (0..2 * @as(usize, @intCast(n))) |i| x[i] = 0;
        }
    };
}

/// The C blocks of the drivers (`#if !defined(HAVE_KERNEL_8)`), and
/// `sscal.c`'s strided one.
fn Fallback(comptime F: type) type {
    return struct {
        pub fn run(n: BlasLong, alpha: [*]const F, x: [*]F) void {
            for (0..@intCast(n)) |i| x[i] *= alpha[0];
        }
        pub fn zero(n: BlasLong, _: [*]const F, x: [*]F) void {
            for (0..@intCast(n)) |i| x[i] = 0;
        }
        pub fn inc(n: BlasLong, alpha: [*]const F, x: [*]F, inc_x: BlasLong) void {
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                const v = util.at(F, x, i * inc_x);
                v.* = alpha[0] * v.*;
            }
        }
    };
}

/// The blocks `slot` calls on the core whose entry is `e`: the microkernel's
/// functions, the architecture's `Vector(slot)` where it is intrinsics, the
/// C where the chain names none; the strided one from the driver's own
/// assembly (`e.own`) where it has it.
fn Blocks(comptime slot: Slot, comptime e: uk.Entry, comptime Vector: fn (comptime Slot) type) type {
    const F = Elem(slot);
    const p = @tagName(slot)[0..1];
    const name = std.fmt.comptimePrint("{s}scal_kernel_{d}", .{ p, width(slot) });
    const inc_name = p ++ "scal_kernel_inc_8";
    const B = switch (e.block) {
        .assembly => struct {
            pub const run = @field(e.module, name);
            pub const zero = @field(e.module, name ++ "_zero");
        },
        .vector => Vector(slot),
        .none => if (complex(slot)) ComplexFallback(F) else Fallback(F),
    };
    return struct {
        pub const run = B.run;
        pub const zero = B.zero;
        pub const inc = if (e.own != void and @hasDecl(e.own, inc_name)) @field(e.own, inc_name) else Fallback(F).inc;
    };
}

/// The driver (`CNAME`) of `slot` around `B`.
fn Driver(comptime slot: Slot, comptime B: type) type {
    const F = Elem(slot);
    const w = width(slot);
    return struct {
        /// `x[i..n)` times `da` from element `j` on, the strided tail.
        fn tail(n: BlasLong, da: F, x: [*]F, inc_x: BlasLong, start_i: BlasLong, start_j: BlasLong) void {
            var i = start_i;
            var j = start_j;
            while (j < n) : (j += 1) {
                util.at(F, x, i).* *= da;
                i += inc_x;
            }
        }

        fn strided(n: BlasLong, da: F, x: [*]F, inc_x: BlasLong) void {
            var alpha = da;
            const n1 = n & -8;
            var i: BlasLong = 0;
            var j: BlasLong = 0;
            if (n1 > 0) {
                B.inc(n1, @ptrCast(&alpha), x, inc_x);
                i = n1 * inc_x;
                j = n1;
            }
            tail(n, da, x, inc_x, i, j);
        }

        fn unit(n: BlasLong, da: F, x: [*]F, comptime zero: bool) void {
            var alpha = da;
            const n1 = n & -w;
            if (n1 > 0) (if (zero) B.zero else B.run)(n1, @ptrCast(&alpha), x);
            var i = n1;
            while (i < n) : (i += 1) {
                const idx: usize = @intCast(i);
                if (zero) x[idx] = 0 else x[idx] *= da;
            }
        }

        pub fn call(n: BlasLong, _: BlasLong, _: BlasLong, da: F, x: [*]F, inc_x: BlasLong, _: [*]F, _: BlasLong, _: [*]F, flag: BlasLong) callconv(.c) i32 {
            if (flag == 1 and da == 0) {
                if (inc_x != 1) strided(n, da, x, inc_x) else unit(n, da, x, false);
            } else if (inc_x != 1) {
                if (da == 0) {
                    const n1 = n & -2;
                    var i: BlasLong = 0;
                    var j: BlasLong = 0;
                    while (j < n1) : (j += 2) {
                        util.at(F, x, i).* = 0;
                        util.at(F, x, i + inc_x).* = 0;
                        i += 2 * inc_x;
                    }
                    while (j < n) : (j += 1) {
                        util.at(F, x, i).* = 0;
                        i += inc_x;
                    }
                } else strided(n, da, x, inc_x);
            } else if (da == 0) unit(n, da, x, true) else unit(n, da, x, false);
            return 0;
        }
    };
}

/// The driver of `cscal.c` and `zscal.c` around `B`: the block over whole
/// blocks of a unit stride, the product (or zeros) everywhere else.
fn ComplexDriver(comptime slot: Slot, comptime B: type) type {
    const F = Elem(slot);
    const w = width(slot);
    return struct {
        pub fn call(n: BlasLong, _: BlasLong, _: BlasLong, da_r: F, da_i: F, x: [*]F, inc_x: BlasLong, _: [*]F, _: BlasLong, _: [*]F, flag: BlasLong) callconv(.c) i32 {
            const zero = da_r == 0 and da_i == 0 and flag == 0;
            var i: BlasLong = 0;
            var j: BlasLong = 0;
            if (inc_x == 1) {
                const n1 = n & -w;
                if (n1 > 0) {
                    const alpha = [2]F{ da_r, da_i };
                    if (zero) B.zero(n1, &alpha, x) else B.run(n1, &alpha, x);
                    i = 2 * n1;
                    j = n1;
                }
            }
            while (j < n) : (j += 1) {
                const re = util.at(F, x, i);
                const im = util.at(F, x, i + 1);
                if (zero) {
                    re.* = 0;
                    im.* = 0;
                } else product(F, da_r, da_i, re, im);
                i += 2 * inc_x;
            }
            return 0;
        }
    };
}

/// The kernel `slot` of the core whose index entry is `e` (`uk_block.Kernels`).
pub fn entryOf(comptime slot: Slot, comptime core: coretypes.X86CoreType, comptime e: uk.Entry, comptime Vector: fn (comptime Slot) type) Kernel(slot) {
    _ = core;
    const B = Blocks(slot, e, Vector);
    return if (complex(slot)) &ComplexDriver(slot, B).call else &Driver(slot, B).call;
}
