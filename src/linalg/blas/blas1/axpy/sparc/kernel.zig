//! SPARC axpy kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`), and how they get installed into
//! the vtable (`install`): cores without an entry keep the portable kernels
//! from `axpy.installGeneric`. Per-core tests live in `axpy/__root__.zig`,
//! not here -- they need the full axpy interface (`checkReal`, ...), which
//! imports outside `src/linalg/blas/` and so cannot appear in this file (see
//! `kernels/__root__.zig`).
//!
//! All of the reference's SPARC axpy kernels are whole-function `.S`, in
//! `asm/` (`asm/__root__.zig` says which core uses which), built once per pointer
//! width since their register choices differ.

const coretypes = @import("../../../cpu/_coretypes.zig");
const registry = @import("../../../asm_registry.zig");
const Vtable = @import("../../../Vtable.zig");
const kernel = @import("../kernel.zig");
const roles = kernel.roles;

const R = registry.Kernels(coretypes.SparcCoreType, @import("asm/__root__.zig"));

pub const Core = coretypes.SparcCoreType;
pub const has = R.has;
pub const hasRole = R.hasRole;
pub const symbol = R.symbol;
pub const exportCore = R.exportCore;

/// Points the vtable at the kernels `core` has.
pub fn install(vt: *Vtable, comptime core: Core) void {
    registry.installRoles(vt, core, roles, hasRole);
}
