//! x86_64 axpy kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`), and how they get installed into
//! the vtable (`install`): cores without an entry keep the portable kernels
//! from `axpy.installGeneric`. Per-core tests live in `axpy/__root__.zig`,
//! not here -- they need the full axpy interface (`checkReal`, ...), which
//! imports outside `src/linalg/blas/` and so cannot appear in this file (see
//! `kernels/__root__.zig`).
//!
//! A core's kernel is one of two things upstream, chosen by its
//! `KERNEL.<CORE>`:
//!
//!   * a whole-function `.S` kernel, in `asm/`, one Zig file per source and
//!     precision (`asm/__root__.zig` says which core uses which), exported as the
//!     kernel directly;
//!   * the C driver `{s,d,c,z}axpy.c` around a block kernel, which is the
//!     shared driver in `../kernel.zig` around a block from `uk/` (the
//!     reference's inline-assembly microkernel, as written), or, for the
//!     AVX-512 cores whose upstream block is C intrinsics, the same loops on
//!     `@Vector`s. `uk/__root__.zig` says which block each core's driver
//!     includes and at what block size. `uk.Kernels` composes them.
//!
//! The build gives the object the core's CPU features; the assembler needs
//! them to accept the mnemonics, while whether the C-derived multiply-adds are
//! contracted is the separate policy in `fma.zig`.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../../../types.zig");
const coretypes = @import("../../../cpu/_coretypes.zig");
const fma = @import("../../../fma.zig");
const kernel = @import("../kernel.zig");
const uk = @import("../../uk_block.zig");

const BlasLong = types.BlasLong;

pub const Core = coretypes.X86CoreType;

const has_avx512f = std.Target.x86.featureSetHas(sysinfo.target.cpu.features, .avx512f);

/// `y += alpha * x` over a whole block, on vectors as wide as the core's
/// widest: what the reference's intrinsics block does (real roles only).
fn VectorBlock(comptime role: []const u8) type {
    const T = uk.Float(role);
    const width = (if (has_avx512f) 64 else 32) / @sizeOf(T);
    const V = @Vector(width, T);
    return struct {
        pub fn run(n1: BlasLong, alpha: T, x: [*]T, y: [*]T) void {
            const av: V = @splat(alpha);
            const len: usize = @intCast(n1);
            var i: usize = 0;
            while (i + width <= len) : (i += width) {
                const xv: V = x[i..][0..width].*;
                const yv: V = y[i..][0..width].*;
                y[i..][0..width].* = fma.madd(fma.fuse, V, av, xv, yv);
            }
        }
    };
}

const K = uk.Kernels(Core, kernel, @import("asm/__root__.zig"), @import("uk/__root__.zig"), VectorBlock);

pub const hasRole = K.hasRole;
pub const has = K.has;
pub const symbol = K.symbol;
pub const exportCore = K.exportCore;
pub const install = K.install;
