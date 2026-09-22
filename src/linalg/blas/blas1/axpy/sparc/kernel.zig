//! SPARC axpy kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`).
//!
//! All of the reference's SPARC axpy kernels are whole-function `.S`, in
//! `asm/` (`asm/__root__.zig` says which core uses which), built once per pointer
//! width since their register choices differ.

const coretypes = @import("../../../cpu/_coretypes.zig");
const registry = @import("../../../asm_registry.zig");
const roles = @import("../kernel.zig").roles;

const width = if (@sizeOf(usize) == 8) "64" else "32";

const R = registry.Registry(coretypes.SparcCoreType, @import("asm/__root__.zig"), width, roles);

pub const Core = R.Core;
pub const has = R.hasAny;
pub const hasRole = R.hasRole;
pub const symbol = R.symbol;
pub const exportCore = R.exportCore;
