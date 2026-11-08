//! BLAS kernel registry for architectures without a registry of their own: no kernels are ported yet, so no per-core
//! object is built and every core keeps the portable kernels. See
//! `x86.zig` for what a registry provides and how to fill one in.

const coretypes = @import("../cpu/coretypes.zig");
const registry = @import("../registry.zig");

pub const Core = coretypes.GenericCoreType;

/// No kernels; the portable ones stay.
pub const Registry = registry.Arch(Core, &.{}, .{}, false);
