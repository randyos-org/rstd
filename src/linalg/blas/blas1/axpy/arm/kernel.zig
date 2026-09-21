//! 32-bit ARM axpy kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`).
//!
//! The reference's one assembly axpy kernel for this architecture,
//! `axpy_vfp.S` (all four precisions), is in `asm/` as a Zig file with its
//! `#if`s and macros resolved for a precision and float ABI. It is emitted as
//! a naked function whose body is that assembly and reads its arguments
//! straight from the ABI registers, as the reference's does. Cores whose
//! reference kernel is C (ARMv5) keep the portable kernels.

const std = @import("std");
const sysinfo = @import("builtin");
const coretypes = @import("../../../cpu/_coretypes.zig");
const asmexport = @import("../../../asm/exports.zig");
const axpy_vfp = @import("asm/axpy_vfp.zig");

pub const Core = coretypes.ArmCoreType;

const roles = [_][]const u8{ "saxpy", "daxpy", "caxpy", "caxpyc", "zaxpy", "zaxpyc" };

/// `__ARM_PCS_VFP`: floating-point arguments are passed in VFP registers.
const hard_float = switch (sysinfo.abi) {
    .eabihf, .gnueabihf, .musleabihf => true,
    else => false,
};

/// Whether `core` has the VFP kernel: its `KERNEL.<CORE>` choice, with
/// `KERNEL.ARMV7` including `KERNEL.ARMV6`.
fn hasVfp(comptime core: Core) bool {
    return core != .armv5;
}

/// The kernel's assembly for `role`.
fn code(comptime role: []const u8, comptime core: Core) []const u8 {
    comptime {
        const double = role[0] == 'd' or role[0] == 'z';
        const complex = role[0] == 'c' or role[0] == 'z';
        const conj = role[role.len - 1] == 'c';
        return axpy_vfp.text(double, complex, conj, hard_float, asmexport.labelSuffix(role, @tagName(core)));
    }
}

/// Whether `core` has any assembly axpy kernel.
pub fn has(comptime core: Core) bool {
    return hasVfp(core);
}

/// Whether `role` has a kernel on `core`.
pub fn hasRole(comptime core: Core, comptime role: []const u8) bool {
    _ = role;
    return hasVfp(core);
}

/// Exported symbol of `role` for `core`.
pub fn symbol(comptime role: []const u8, comptime core: Core) []const u8 {
    return asmexport.symbol(role, @tagName(core));
}

/// Exports every axpy kernel `core` has. Call from a container-level `comptime`
/// block of the object's root file.
pub fn exportCore(comptime core: Core) void {
    if (comptime !hasVfp(core)) return;
    inline for (roles) |role| {
        @export(&asmexport.Naked(code(role, core), 4).f, .{ .name = symbol(role, core) });
    }
}
