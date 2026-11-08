//! BLAS kernel registry for wasm: which routines have kernels for which
//! cores, how the per-core objects (`_kernel.zig`) export them, and how
//! `cpu.Arch` points the vtable at them. The build (`build/blaskern.zig`)
//! asks it which cores need an object at all.
//!
//! Add a routine to `modules` to wire it in.

const coretypes = @import("../cpu/coretypes.zig");
const registry = @import("../registry.zig");

pub const Core = coretypes.WasmCoreType;

/// Every routine with kernels here: `axpy`'s own module, which adds the
/// kernels upstream writes in C, and gemm's.
const modules = .{
    @import("../blas1/axpy/wasm/kernel.zig"),
    @import("../blas3/kernels.zig").Registry(Core),
};

/// The registry (`registry.Arch`).
pub const Registry = registry.Arch(Core, &.{ .wasm32, .wasm64 }, modules, false);
