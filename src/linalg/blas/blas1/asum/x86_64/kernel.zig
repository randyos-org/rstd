//! x86_64 asum kernels: the `@Vector` ports of `kernel/x86_64/sasum.c` and
//! `dasum.c` with their Haswell and Skylake-X intrinsics blocks, and of
//! `casum.c` and `zasum.c` with Skylake-X's. The `.S` kernels of the other
//! cores peel to an aligned address, so their slots take `vector.zig`'s
//! (`vector.peels`; `asm/`'s `skipped`). Haswell's blocks peel too
//! (`sasum_microk_haswell-2.c`, `dasum_microk_haswell-2.c`); here they load
//! unaligned instead, as upstream's Skylake-X blocks do, so no sum depends
//! on an address. Tests live in `asum/__root__.zig`; the AVX-512 kernels run
//! only on such hardware (T43).

const std = @import("std");
const coretypes = @import("../../../cpu/_coretypes.zig");
const registry = @import("../../../asm_registry.zig");
const types = @import("../../../types.zig");
const util = @import("../../../util.zig");
const Vtable = @import("../../../Vtable.zig");

const BlasLong = types.BlasLong;
const Slot = util.Slot;

pub const Core = coretypes.X86CoreType;

/// The real kernels' first block size and the rows above which it runs
/// (`n & ~255`, `n >= 256`).
const big = 256;

/// `?asum.c`'s `asum_compute` around `sum`, the block for a unit stride; a
/// strided `x` sums in order. `w` floats an element.
fn Compute(comptime F: type, comptime w: comptime_int, comptime sum: fn (usize, [*]const F) F) type {
    return struct {
        pub fn k(n: BlasLong, x: [*]F, inc_x: BlasLong) callconv(.c) F {
            if (n <= 0 or inc_x <= 0) return 0;
            if (inc_x == 1) return sum(@intCast(n), x);
            var s: F = 0;
            var i: usize = 0;
            const step: usize = @intCast(w * inc_x);
            for (0..@intCast(n)) |_| {
                if (w == 1) s += @abs(x[i]) else s += @abs(x[i]) + @abs(x[i + 1]);
                i += step;
            }
            return s;
        }
    };
}

/// `n` vectors of `l` floats from `p` into `acc`, vector `q` into sum `q % 4`.
inline fn add(comptime F: type, comptime l: comptime_int, comptime n: comptime_int, acc: *[4]@Vector(l, F), p: [*]const F) void {
    inline for (0..n) |q| acc[q % 4] += @abs(util.load(@Vector(l, F), p + q * l));
}

/// The real blocks: four sums of `wide` floats over blocks of 256 (Haswell's
/// 256 bits, Skylake-X's 512), reduced by `hsum` or AVX-512's halving
/// `_mm512_reduce_add`; then 128-bit sums to a multiple of 8, then the rest.
fn RealBlock(comptime F: type, comptime wide: comptime_int, comptime halving: bool) type {
    return struct {
        pub fn sum(n: usize, x: [*]const F) F {
            var s: F = 0;
            const tail_sse = n & ~@as(usize, 7);
            const tail_big = n & ~@as(usize, big - 1);
            if (n >= big) {
                var acc: [4]@Vector(wide, F) = @splat(@splat(0));
                var i: usize = 0;
                while (i < tail_big) : (i += 4 * wide) add(F, wide, 4, &acc, x + i);
                const a = acc[0] + acc[1] + acc[2] + acc[3];
                s += if (halving) util.hsumFold(a) else util.hsum(a);
            }
            if (n >= 8) {
                // `sasum`'s two sums of four, `dasum`'s four of two
                const l = 16 / @sizeOf(F);
                const k = 8 / l;
                var acc: [4]@Vector(l, F) = @splat(@splat(0));
                var i = tail_big;
                while (i < tail_sse) : (i += 8) add(F, l, k, &acc, x + i);
                var a = acc[0];
                inline for (1..k) |q| a += acc[q];
                s += util.hsum(a);
            }
            for (tail_sse..n) |i| s += @abs(x[i]);
            return s;
        }
    };
}

/// `casum_microk_skylakex-2.c`, `zasum_microk_skylakex-2.c`: `n` elements,
/// `2n` floats. Under 4 vectors of 512 bits, 128-bit sums; else 512-bit ones
/// eight vectors a step, then four, two, one and a masked one.
fn ComplexBlock(comptime F: type) type {
    const big_l = 64 / @sizeOf(F);
    const l = 16 / @sizeOf(F);
    return struct {
        pub fn sum(n: usize, x0: [*]const F) F {
            var n2 = 2 * n;
            var x = x0;
            if (n2 < 4 * big_l) {
                var acc: [4]@Vector(l, F) = @splat(@splat(0));
                inline for (.{ 8, 4, 2, 1 }) |v| {
                    if (n2 >= v * l) {
                        add(F, l, v, &acc, x);
                        n2 -= v * l;
                        x += v * l;
                    }
                }
                // `casum`'s last element, apart; `zasum`'s took a vector
                const rest: F = if (n2 != 0) @abs(x[0]) + @abs(x[1]) else 0;
                return rest + util.hsum((acc[0] + acc[1]) + (acc[2] + acc[3]));
            }
            var acc: [4]@Vector(big_l, F) = @splat(@splat(0));
            while (n2 >= 8 * big_l) : (n2 -= 8 * big_l) {
                add(F, big_l, 8, &acc, x);
                x += 8 * big_l;
            }
            inline for (.{ 4, 2, 1 }) |v| {
                if (n2 >= v * big_l) {
                    add(F, big_l, v, &acc, x);
                    n2 -= v * big_l;
                    x += v * big_l;
                }
            }
            if (n2 != 0) {
                // `_mm512_maskz_loadu`: the lanes past the end are zero
                var t: [big_l]F = @splat(0);
                for (0..n2) |i| t[i] = x[i];
                acc[0] += @abs(@as(@Vector(big_l, F), t));
            }
            return util.hsumFold((acc[0] + acc[1]) + (acc[2] + acc[3]));
        }
    };
}

/// The kernels of the Haswell blocks, or the Skylake-X ones.
pub fn Ported(comptime avx512: bool) type {
    return struct {
        pub fn entry(comptime slot: Slot) util.SlotFn(slot) {
            const wide = if (avx512) 64 else 32;
            return switch (slot) {
                .sasum_k => &Compute(f32, 1, RealBlock(f32, wide / 4, avx512).sum).k,
                .dasum_k => &Compute(f64, 1, RealBlock(f64, wide / 8, avx512).sum).k,
                .casum_k => &Compute(f32, 2, ComplexBlock(f32).sum).k,
                .zasum_k => &Compute(f64, 2, ComplexBlock(f64).sum).k,
                else => @compileError(@tagName(slot) ++ " is not an asum slot"),
            };
        }
    };
}

/// The ports and where they serve (`KERNEL.<CORE>`'s `?ASUMKERNEL`).
pub const ports = [_]registry.Port(Core){
    .{ .module = Ported(false), .cores = &.{.haswell}, .slots = &.{ .sasum_k, .dasum_k }, .written = .vector },
    .{ .module = Ported(true), .cores = &.{ .skylakex, .cooperlake, .sapphirerapids }, .slots = &.{ .sasum_k, .dasum_k, .casum_k, .zasum_k }, .written = .vector },
};

const P = registry.Ports(Core, "x86_64", &ports);
const A = registry.Kernels(Core, @import("asm/__root__.zig"));

pub fn hasSlot(comptime core: Core, comptime slot: Slot) bool {
    return A.hasSlot(core, slot) or P.hasSlot(core, slot);
}

pub fn has(comptime core: Core) bool {
    return A.has(core) or P.has(core);
}

pub const skips = A.skips;

pub fn exportCore(comptime core: Core) void {
    A.exportCore(core);
    P.exportCore(core);
}

pub fn install(vt: *Vtable, comptime core: Core) void {
    A.install(vt, core);
    P.install(vt, core);
}
