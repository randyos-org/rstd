//! The drivers of `kernel/x86_64/sscal.c` and `dscal.c`, around the blocks a
//! core's `uk/` index names: the unit-stride block, its zeroing twin, and the
//! strided one (`dscal.c`'s own assembly, `sscal.c`'s C). Every element is one
//! product or a stored zero, so any block gives the same bits. Nothing here
//! depends on the BLAS root, so the per-core objects (`_kernel.zig`) can use it.
//!
//! The last argument (`dummy2`) is the flag of the portable kernels: 1, which
//! the interface passes, multiplies even by a zero `alpha` (NaN stays NaN);
//! otherwise a zero `alpha` stores zeros.

const std = @import("std");
const util = @import("../../util.zig");
const types = @import("../../types.zig");
const uk = @import("../uk_block.zig");

const BlasLong = types.BlasLong;
const Slot = util.Slot;

/// The kernels an architecture's `<arch>/kernel.zig` iterates over: the real
/// ones driven here, and the complex ones only converted assembly fills.
pub const slots = [_]Slot{ .sscal_k, .dscal_k, .cscal_k, .zscal_k };

/// The element type of `slot`, one of `slots`.
pub fn Elem(comptime slot: Slot) type {
    return switch (slot) {
        .sscal_k => f32,
        .dscal_k => f64,
        else => @compileError(@tagName(slot) ++ " is not a real scal kernel"),
    };
}

/// Vtable slot type of `slot`.
pub fn Kernel(comptime slot: Slot) type {
    return util.SlotFn(slot);
}

/// The unit-stride block size: `sscal_kernel_16`, `dscal_kernel_8`.
pub fn width(comptime slot: Slot) comptime_int {
    return if (slot == .sscal_k) 16 else 8;
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
    const p = if (F == f32) "s" else "d";
    const name = std.fmt.comptimePrint("{s}scal_kernel_{d}", .{ p, width(slot) });
    const inc_name = p ++ "scal_kernel_inc_8";
    const B = switch (e.block) {
        .assembly => struct {
            pub const run = @field(e.module, name);
            pub const zero = @field(e.module, name ++ "_zero");
        },
        .vector => Vector(slot),
        .none => Fallback(F),
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

        pub fn call(n: BlasLong, _: BlasLong, _: BlasLong, da: F, x: [*]F, inc_x: BlasLong, _: [*]F, _: BlasLong, _: [*]F, dummy2: BlasLong) callconv(.c) i32 {
            if (dummy2 == 1 and da == 0) {
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

/// The kernel `slot` of the core whose index entry is `e` (`uk_block.driven`).
pub fn entryOf(comptime slot: Slot, comptime e: uk.Entry, comptime Vector: fn (comptime Slot) type) Kernel(slot) {
    return &Driver(slot, Blocks(slot, e, Vector)).call;
}
