//! RISC-V axpy kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`), and how they get installed into
//! the vtable (`install`): cores without an entry keep the portable kernels
//! from `axpy.installGeneric`. Per-core tests live in `axpy/__root__.zig`,
//! not here -- they need the full axpy interface (`checkReal`, ...), which
//! imports outside `src/linalg/blas/` and so cannot appear in this file (see
//! `kernels/__root__.zig`).
//!
//! The reference's kernels here are C over the vector extension's
//! intrinsics, `axpy_rvv.c` and `zaxpy_rvv.c` for the cores whose
//! `KERNEL.<CORE>` names them; they are in `rvv.zig` as the assembly the
//! intrinsics stand for. The cores that use the plain C loops
//! (`riscv64_generic`, `u74`) keep the portable kernels.

const coretypes = @import("../../../cpu/_coretypes.zig");
const Vtable = @import("../../../Vtable.zig");
const util = @import("../../../util.zig");
const kernel = @import("../kernel.zig");
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
