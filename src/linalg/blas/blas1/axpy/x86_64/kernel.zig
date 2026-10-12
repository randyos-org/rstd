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
//!     includes and at what block size.
//!
//! The build gives the object the core's CPU features; the assembler needs
//! them to accept the mnemonics, while whether the C-derived multiply-adds are
//! contracted is the separate policy in `fma.zig`.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../../../types.zig");
const coretypes = @import("../../../cpu/_coretypes.zig");
const fma = @import("../../../fma.zig");
const util = @import("../../../util.zig");
const registry = @import("../../../asm_registry.zig");
const kernel = @import("../kernel.zig");
const uk = @import("../uk_block.zig");
const ukidx = @import("uk/__root__.zig");
const Vtable = @import("../../../Vtable.zig");

const BlasLong = types.BlasLong;

pub const Core = coretypes.X86CoreType;
const roles = kernel.roles;

const asm_kernels = registry.Registry(Core, @import("asm/__root__.zig"), "", roles);

/// The C-driven entry for `role` on `core`, if it has one.
fn cEntry(comptime role: []const u8, comptime core: Core) ?ukidx.Entry {
    return ukidx.find(role, @tagName(core));
}

const has_avx512f = std.Target.x86.featureSetHas(sysinfo.cpu.features, .avx512f);

/// `y += alpha * x` over a whole block, on vectors as wide as the core's
/// widest: what the reference's intrinsics block does.
fn VectorBlock(comptime T: type) type {
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
                y[i..][0..width].* = kernel.fmadd(fma.fuse, V, av, xv, yv);
            }
        }
    };
}

/// Whether `role` has a kernel of its own on `core`, either kind.
pub fn hasRole(comptime core: Core, comptime role: []const u8) bool {
    return asm_kernels.hasRole(core, role) or cEntry(role, core) != null;
}

/// Whether `core` has any axpy kernel of its own.
pub fn has(comptime core: Core) bool {
    inline for (roles) |role| {
        if (comptime hasRole(core, role)) return true;
    }
    return false;
}

/// Exported symbol of `role` for `core`.
pub fn symbol(comptime role: []const u8, comptime core: Core) []const u8 {
    return util.kernelSymbol(role, @tagName(core));
}

/// Exports every axpy kernel `core` has. Call from a container-level `comptime`
/// block of the object's root file.
pub fn exportCore(comptime core: Core) void {
    inline for (roles) |role| {
        if (comptime asm_kernels.hasRole(core, role)) {
            asm_kernels.exportOne(role, core);
        } else if (comptime cEntry(role, core)) |e| {
            @export(uk.driven(role, e, VectorBlock), .{ .name = symbol(role, core) });
        }
    }
}

fn ext(comptime K: type, comptime name: []const u8) K {
    return @extern(K, .{ .name = name });
}

fn installCore(vt: *Vtable, comptime core: Core) void {
    if (comptime hasRole(core, "daxpy")) vt.daxpy_k = ext(kernel.RealKernel(f64), symbol("daxpy", core));
    if (comptime hasRole(core, "saxpy")) vt.saxpy_k = ext(kernel.RealKernel(f32), symbol("saxpy", core));
    if (comptime hasRole(core, "caxpy")) {
        vt.caxpy_k = ext(kernel.ComplexKernel(f32), symbol("caxpy", core));
        vt.caxpyc_k = ext(kernel.ComplexKernel(f32), symbol("caxpyc", core));
    }
    if (comptime hasRole(core, "zaxpy")) {
        vt.zaxpy_k = ext(kernel.ComplexKernel(f64), symbol("zaxpy", core));
        vt.zaxpyc_k = ext(kernel.ComplexKernel(f64), symbol("zaxpyc", core));
    }
}

/// Overrides the generic kernels installed by `axpy.install` for `core`.
pub fn install(vt: *Vtable, core: Core) void {
    switch (core) {
        inline .prescott,
        .core2,
        .penryn,
        .dunnington,
        .opteron,
        .opteron_sse3,
        .barcelona,
        .bobcat,
        .nano,
        .atom,
        .nehalem,
        .sandybridge,
        .bulldozer,
        .piledriver,
        .steamroller,
        .excavator,
        .haswell,
        .zen,
        .skylakex,
        .cooperlake,
        .sapphirerapids,
        => |c| installCore(vt, c),
        else => {},
    }
}
