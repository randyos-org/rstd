//! RISC-V axpy kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`).
//!
//! The reference's kernels here are C over the vector extension's
//! intrinsics, `axpy_rvv.c` and `zaxpy_rvv.c` for the cores whose
//! `KERNEL.<CORE>` names them; they are in `rvv.zig` as the assembly the
//! intrinsics stand for. The cores that use the plain C loops
//! (`riscv64_generic`, `u74`) keep the portable kernels.

const coretypes = @import("../../../cpu/_coretypes.zig");
const util = @import("../../../util.zig");
const rvv = @import("rvv.zig");

pub const Core = coretypes.RiscvCoreType;

/// Whether `core` uses the `_rvv` kernels.
fn usesRvv(comptime core: Core) bool {
    return switch (core) {
        .riscv64_zvl128b, .x280 => true,
        else => false,
    };
}

/// Whether `role` has a kernel of its own on `core`.
pub fn hasRole(comptime core: Core, comptime role: []const u8) bool {
    _ = role;
    return usesRvv(core);
}

/// Whether `core` has any axpy kernel of its own.
pub fn has(comptime core: Core) bool {
    return usesRvv(core);
}

/// Exported symbol of `role` for `core`.
pub fn symbol(comptime role: []const u8, comptime core: Core) []const u8 {
    return util.kernelSymbol(role, @tagName(core));
}

/// Exports every axpy kernel `core` has. Call from a container-level `comptime`
/// block of the object's root file.
pub fn exportCore(comptime core: Core) void {
    if (comptime !usesRvv(core)) return;
    @export(&rvv.Real(f32).call, .{ .name = symbol("saxpy", core) });
    @export(&rvv.Real(f64).call, .{ .name = symbol("daxpy", core) });
    @export(&rvv.Cplx(f32, false).call, .{ .name = symbol("caxpy", core) });
    @export(&rvv.Cplx(f32, true).call, .{ .name = symbol("caxpyc", core) });
    @export(&rvv.Cplx(f64, false).call, .{ .name = symbol("zaxpy", core) });
    @export(&rvv.Cplx(f64, true).call, .{ .name = symbol("zaxpyc", core) });
}
