//! BLAS kernel registry for s390x: which routines have kernels for which
//! cores, how the per-core objects (`_kernel.zig`) export them, and how
//! `cpu.Arch` points the vtable at them. The build (`build/blaskern.zig`)
//! asks it which cores need an object at all.
//!
//! Add a routine to `modules` to wire it in.

const coretypes = @import("../cpu/coretypes.zig");
const registry = @import("../registry.zig");

pub const Core = coretypes.S390xCoreType;

/// Every routine with kernels here: `axpy`'s own module, which adds the
/// kernels upstream writes in C, and gemm's.
const modules = .{
    @import("../blas1/axpy/s390x/kernel.zig"),
    @import("../blas3/kernels.zig").Registry(Core),
};

/// The registry (`registry.Arch`).
pub const Registry = registry.Arch(Core, &.{.s390x}, modules, false);
