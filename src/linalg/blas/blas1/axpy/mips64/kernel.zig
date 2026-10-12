//! MIPS axpy kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`), for both word sizes, and how they
//! get installed into the vtable (`install`): cores without an entry keep
//! the portable kernels from `axpy.installGeneric`. Per-core tests live in
//! `axpy/__root__.zig`, not here -- they need the full axpy interface
//! (`checkReal`, ...), which imports outside `src/linalg/blas/` and so
//! cannot appear in this file (see `kernels/__root__.zig`).
//!
//! What a core uses is its `KERNEL.<CORE>` choice, made in the reference by
//! `NO_MSA` where a core may or may not have the SIMD architecture: the SIMD
//! cores use `mips/{s,d,c,z}axpy_msa.c` (`msa.zig`), the Loongson 3A R3 the
//! whole-function `.S` kernels in `asm/` (`asm/__root__.zig`). The cores with
//! neither keep the portable kernels.

const sysinfo = @import("builtin");
const coretypes = @import("../../../cpu/_coretypes.zig");
const registry = @import("../../../asm_registry.zig");
const Vtable = @import("../../../Vtable.zig");
const kregistry = @import("../../../kernels/mips.zig");
const kernel = @import("../kernel.zig");
const roles = kernel.roles;
const msa = @import("../mips/msa.zig");

const R = registry.Registry(coretypes.MipsCoreType, @import("asm/__root__.zig"), "", roles);

pub const Core = R.Core;

/// Whether `core` uses the SIMD kernels.
fn usesMsa(comptime core: Core) bool {
    return switch (core) {
        .loongson3r4, .p5600, .i6400, .i6500, .p6600 => true,
        else => false,
    };
}

/// Whether `role` has a kernel of its own on `core`.
pub fn hasRole(comptime core: Core, comptime role: []const u8) bool {
    return usesMsa(core) or R.hasRole(core, role);
}

/// Whether `core` has any axpy kernel of its own.
pub fn has(comptime core: Core) bool {
    return usesMsa(core) or R.hasAny(core);
}

/// Exported symbol of `role` for `core`.
pub const symbol = R.symbol;

/// Exports every axpy kernel `core` has. Call from a container-level `comptime`
/// block of the object's root file.
pub fn exportCore(comptime core: Core) void {
    if (comptime usesMsa(core)) {
        @export(&msa.Real(f32).call, .{ .name = symbol("saxpy", core) });
        @export(&msa.Real(f64).call, .{ .name = symbol("daxpy", core) });
        @export(&msa.Cplx(f32, false).call, .{ .name = symbol("caxpy", core) });
        @export(&msa.Cplx(f32, true).call, .{ .name = symbol("caxpyc", core) });
        @export(&msa.Cplx(f64, false).call, .{ .name = symbol("zaxpy", core) });
        @export(&msa.Cplx(f64, true).call, .{ .name = symbol("zaxpyc", core) });
        return;
    }
    R.exportCore(core);
}

fn ext(comptime K: type, comptime name: []const u8) K {
    return @extern(K, .{ .name = name });
}

fn installCore(vt: *Vtable, comptime core: Core) void {
    // no object is linked for a core that does not run on this CPU
    if (comptime !kregistry.runsOn(sysinfo.cpu, core)) return;
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
        inline else => |c| installCore(vt, c),
    }
}
