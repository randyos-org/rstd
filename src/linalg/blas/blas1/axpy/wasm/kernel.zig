//! WebAssembly axpy kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`).
//!
//! Only `saxpy` and `daxpy` have a kernel of their own (`simd.zig`, from
//! `kernel/wasm/axpy.c`, for the core with SIMD128); the complex kernels
//! upstream are the plain C loops, which the portable kernels are.

const coretypes = @import("../../../cpu/_coretypes.zig");
const util = @import("../../../util.zig");
const simd = @import("simd.zig");

pub const Core = coretypes.WasmCoreType;

/// Whether `role` has a kernel of its own on `core`.
pub fn hasRole(comptime core: Core, comptime role: []const u8) bool {
    return core == .wasm128_generic and (role[0] == 's' or role[0] == 'd') and role.len == 5;
}

/// Whether `core` has any axpy kernel of its own.
pub fn has(comptime core: Core) bool {
    return core == .wasm128_generic;
}

/// Exported symbol of `role` for `core`.
pub fn symbol(comptime role: []const u8, comptime core: Core) []const u8 {
    return util.kernelSymbol(role, @tagName(core));
}

/// Exports every axpy kernel `core` has. Call from a container-level `comptime`
/// block of the object's root file.
pub fn exportCore(comptime core: Core) void {
    if (comptime core != .wasm128_generic) return;
    @export(&simd.Real(f32).call, .{ .name = symbol("saxpy", core) });
    @export(&simd.Real(f64).call, .{ .name = symbol("daxpy", core) });
}
