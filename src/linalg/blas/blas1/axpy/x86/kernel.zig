//! 32-bit x86 axpy kernels, one set per core, built into the per-core kernel objects
//! (`src/linalg/blas/_kernel.zig`), and how they get installed into the
//! vtable (`install`): cores without an entry keep the portable kernels from
//! `axpy.installGeneric`. Per-core tests live in `axpy/__root__.zig`, not
//! here -- they need the full axpy interface (`checkReal`, ...), which
//! imports outside `src/linalg/blas/` and so cannot appear in this file (see
//! `kernels/__root__.zig`).
//!
//! The kernels that are assembly upstream are in `asm/`, one Zig file per
//! source and precision with its `#if`s and macros resolved (`asm/__root__.zig`
//! says which core uses which). A core whose reference kernel is C keeps the
//! portable kernels until that is ported.

const coretypes = @import("../../../cpu/_coretypes.zig");
const registry = @import("../../../asm_registry.zig");
const Vtable = @import("../../../Vtable.zig");
const kernel = @import("../kernel.zig");
const roles = kernel.roles;

const R = registry.Registry(coretypes.X86CoreType, @import("asm/__root__.zig"), "", roles);

pub const Core = R.Core;
pub const has = R.hasAny;
pub const hasRole = R.hasRole;
pub const symbol = R.symbol;
pub const exportCore = R.exportCore;

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
        inline else => |c| installCore(vt, c),
    }
}
