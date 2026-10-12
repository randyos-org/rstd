//! WebAssembly axpy kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`), and how they get installed into
//! the vtable (`install`): cores without an entry keep the portable kernels
//! from `axpy.installGeneric`. Per-core tests live in `axpy/__root__.zig`,
//! not here -- they need the full axpy interface (`checkReal`, ...), which
//! imports outside `src/linalg/blas/` and so cannot appear in this file (see
//! `kernels/__root__.zig`).
//!
//! Only `saxpy` and `daxpy` have a kernel of their own (`simd.zig`, from
//! `kernel/wasm/axpy.c`, for the core with SIMD128); the complex kernels
//! upstream are the plain C loops, which the portable kernels are.

const coretypes = @import("../../../cpu/_coretypes.zig");
const registry = @import("../../../asm_registry.zig");
const Vtable = @import("../../../Vtable.zig");
const util = @import("../../../util.zig");
const kernel = @import("../kernel.zig");
const roles = kernel.roles;
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
    @export(&simd.Real(f32).call, .{ .name = symbol("saxpy", core), .visibility = .hidden });
    @export(&simd.Real(f64).call, .{ .name = symbol("daxpy", core), .visibility = .hidden });
}

/// Points the vtable at the kernels `core` has.
pub fn install(vt: *Vtable, comptime core: Core) void {
    registry.installRoles(vt, core, roles, hasRole);
}
