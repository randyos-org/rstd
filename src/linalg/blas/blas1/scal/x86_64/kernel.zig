//! x86_64 scal kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`). A core's kernel is a `.S` kernel
//! (`asm/`) or the C driver `{s,d,c,z}scal.c` (`../driver.zig`) around the
//! blocks `uk/__root__.zig` names: the reference's inline-assembly
//! microkernels, or for the AVX-512 cores, whose blocks are intrinsics, the
//! same loops on `@Vector`s. Tests live in `scal/__root__.zig`.

const types = @import("../../../types.zig");
const coretypes = @import("../../../cpu/coretypes.zig");
const util = @import("../../../util.zig");
const kernel = @import("../driver.zig");
const uk = @import("../../../uk_block.zig");

const BlasLong = types.BlasLong;
const Slot = util.Slot;

pub const Core = coretypes.X86CoreType;

/// The reference compiles the blocks with 512-bit vectors when
/// `__AVX512CD__` is defined, else 256-bit ones.
const has_avx512cd = util.x86Has(.avx512cd);

/// `{s,d,c,z}scal_microk_skylakex-2.c`: `x *= alpha` and its zeroing twin
/// over a whole block, one vector at a time. A complex one multiplies each
/// element by `alpha` as `ar*x + (-ai, ai)*swap(x)`, which is the full product.
fn VectorBlock(comptime slot: Slot) type {
    const F = kernel.Elem(slot);
    const lanes = (if (has_avx512cd) 64 else 32) / @sizeOf(F);
    const V = @Vector(lanes, F);
    const cs = if (kernel.complex(slot)) 2 else 1;
    const swap: [lanes]i32 = blk: {
        var m: [lanes]i32 = undefined;
        for (&m, 0..) |*e, i| e.* = @intCast(i ^ 1);
        break :blk m;
    };
    const sign: V = blk: {
        var s: [lanes]F = undefined;
        for (&s, 0..) |*e, i| e.* = if (i % 2 == 0) -1 else 1;
        break :blk s;
    };
    return struct {
        pub fn run(n: BlasLong, alpha: [*]const F, x: [*]F) void {
            const len = cs * @as(usize, @intCast(n));
            if (comptime cs == 2) {
                const ar: V = @splat(alpha[0]);
                const ai = @as(V, @splat(alpha[1])) * sign;
                var i: usize = 0;
                while (i < len) : (i += lanes) {
                    const xv = util.load(V, x + i);
                    util.store(V, x + i, ar * xv + ai * @shuffle(F, xv, undefined, swap));
                }
                return;
            }
            const a: V = @splat(alpha[0]);
            var i: usize = 0;
            while (i < len) : (i += lanes) util.store(V, x + i, a * util.load(V, x + i));
        }
        pub fn zero(n: BlasLong, _: [*]const F, x: [*]F) void {
            const len = cs * @as(usize, @intCast(n));
            var i: usize = 0;
            while (i < len) : (i += lanes) util.store(V, x + i, @splat(0));
        }
    };
}

const K = uk.Kernels(Core, kernel, @import("asm/__root__.zig"), @import("uk/__root__.zig"), VectorBlock);

pub const hasSlot = K.hasSlot;
pub const skips = K.skips;
pub const has = K.has;
pub const exportCore = K.exportCore;
pub const install = K.install;
