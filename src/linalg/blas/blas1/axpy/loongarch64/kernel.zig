//! LoongArch64 axpy kernels, one set per core, built into the per-core kernel objects
//! (`src/linalg/blas/_kernel.zig`).
//!
//! The kernels that are assembly upstream are in `asm/`, one Zig file per
//! source and precision with its `#if`s and macros resolved (`asm/__root__.zig`
//! says which core uses which). A core whose reference kernel is C keeps the
//! portable kernels until that is ported.

const coretypes = @import("../../../cpu/_coretypes.zig");
const registry = @import("../../../asm_registry.zig");
const roles = @import("../kernel.zig").roles;

const R = registry.Registry(coretypes.LoongarchCoreType, @import("asm/__root__.zig"), "", roles);

pub const Core = R.Core;
pub const has = R.hasAny;
pub const hasRole = R.hasRole;
pub const symbol = R.symbol;
pub const exportCore = R.exportCore;
