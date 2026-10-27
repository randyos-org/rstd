//! x86_64 dot kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`), and how they get installed into
//! the vtable (`install`): cores without an entry keep the portable kernels
//! from `dot.installGeneric`. Per-core tests live in `dot/__root__.zig`.
//!
//! A core's kernel is one of two things upstream, chosen by its
//! `KERNEL.<CORE>`: a whole-function `.S` kernel (`asm/`), or the C driver
//! `{s,d,c,z}dot.c`, which is the shared driver in `../kernel.zig` around a
//! block from `uk/` (the reference's inline-assembly microkernel, as written)
//! or, for the AVX-512 cores whose upstream block is C intrinsics, the same
//! reductions on `@Vector`s. `uk/__root__.zig` says which block each core's
//! driver includes. `uk.Kernels` composes them.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../../../types.zig");
const coretypes = @import("../../../cpu/_coretypes.zig");
const fma = @import("../../../fma.zig");
const util = @import("../../../util.zig");
const kernel = @import("../kernel.zig");
const uk = @import("../../uk_block.zig");

const BlasLong = types.BlasLong;

pub const Core = coretypes.X86CoreType;

/// The reference compiles these blocks as intrinsics when `__AVX512CD__` is
/// defined, with a 512-bit pass ahead of the 256-bit one.
const has_avx512cd = std.Target.x86.featureSetHas(sysinfo.target.cpu.features, .avx512cd);

/// The lower (`high` false) or upper half of `v`.
inline fn half(comptime T: type, v: anytype, comptime high: bool) @Vector(@typeInfo(@TypeOf(v)).vector.len / 2, T) {
    const n = @typeInfo(@TypeOf(v)).vector.len / 2;
    const mask = comptime blk: {
        var m: [n]i32 = undefined;
        for (&m, 0..) |*e, i| e.* = @intCast(i + if (high) n else 0);
        break :blk m;
    };
    return @shuffle(T, v, undefined, mask);
}

/// `ddot_microk_skylakex-2.c` and `sdot_microk_skylakex-2.c`, for the roles
/// that use them (`sdot`, `dsdot`, `ddot`): four accumulators of 256 bits, fed
/// ahead of that by four of 512 when `wide` (folded into the 256-bit ones), the
/// accumulators added in order, then the halves, then `hadd`. `fuse` is whether
/// `acc += a*b` contracts.
pub fn VectorBlockWith(comptime role: []const u8, comptime wide: bool, comptime fuse: bool) type {
    const T = uk.Float(role);
    const l256 = 32 / @sizeOf(T);
    const l512 = 64 / @sizeOf(T);
    const V = @Vector(l256, T);
    const W = @Vector(l512, T);
    return struct {
        pub fn run(n1: BlasLong, x: [*]T, y: [*]T) T {
            const n: usize = @intCast(n1);
            var acc: [4]V = @splat(@splat(0));
            var i: usize = 0;
            if (comptime wide) {
                var w: [4]W = @splat(@splat(0));
                const step = 4 * l512;
                const nw = n & ~@as(usize, step - 1);
                while (i < nw) : (i += step) {
                    inline for (0..4) |g| {
                        w[g] = fma.madd(fuse, W, util.load(W, x + i + g * l512), util.load(W, y + i + g * l512), w[g]);
                    }
                }
                inline for (0..4) |g| acc[g] = half(T, w[g], false) + half(T, w[g], true);
            }
            while (i < n) : (i += 4 * l256) {
                inline for (0..4) |g| {
                    acc[g] = fma.madd(fuse, V, util.load(V, x + i + g * l256), util.load(V, y + i + g * l256), acc[g]);
                }
            }
            const a = ((acc[0] + acc[1]) + acc[2]) + acc[3];
            const h = half(T, a, false) + half(T, a, true);
            // `_mm_hadd_pd` once, `_mm_hadd_ps` twice
            if (comptime T == f64) return h[0] + h[1];
            return (h[0] + h[1]) + (h[2] + h[3]);
        }
    };
}

/// The block for the core being built.
fn VectorBlock(comptime role: []const u8) type {
    return VectorBlockWith(role, has_avx512cd, fma.fuse);
}

const K = uk.Kernels(Core, kernel, @import("asm/__root__.zig"), @import("uk/__root__.zig"), VectorBlock);

pub const hasRole = K.hasRole;
pub const has = K.has;
pub const symbol = K.symbol;
pub const exportCore = K.exportCore;
pub const install = K.install;
