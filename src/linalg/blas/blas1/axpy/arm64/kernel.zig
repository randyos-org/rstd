//! AArch64 axpy kernels, one set per core, built into the per-core kernel objects
//! (`src/linalg/blas/_kernel.zig`), and how they get installed into the
//! vtable (`install`): cores without an entry keep the portable kernels from
//! `axpy.installGeneric`. Per-core tests live in `axpy/__root__.zig`, not
//! here -- they need the full axpy interface (`checkReal`, ...), which
//! imports outside `src/linalg/blas/` and so cannot appear in this file (see
//! `kernels/__root__.zig`).
//!
//! The kernels that are assembly upstream are in `asm/`, one Zig file per
//! source and precision with its `#if`s and macros resolved (`asm/__root__.zig`
//! says which core uses which). `A64FX` and `NEOVERSEV1` use `axpy_sve.c`
//! (`sve.zig`); `THUNDERX` uses `daxpy_thunderx.c`, its block (`thunderx.zig`)
//! under the shared driver. Other cores keep the portable kernels.

const std = @import("std");
const sysinfo = @import("builtin");
const coretypes = @import("../../../cpu/_coretypes.zig");
const registry = @import("../../../asm_registry.zig");
const Vtable = @import("../../../Vtable.zig");
const kernel = @import("../kernel.zig");
const roles = kernel.roles;
const uk = @import("../../uk_block.zig");
const sve = @import("sve.zig");
const thunderx = @import("thunderx.zig");

const R = registry.Kernels(coretypes.Arm64CoreType, @import("asm/__root__.zig"));

pub const Core = coretypes.Arm64CoreType;
pub const symbol = R.symbol;

const thunderx_entry = .{
    .mask = 32,
    .block = uk.Block.assembly,
    .kernel = "daxpy_kernel_8",
    .module = thunderx,
};

fn usesSve(comptime core: Core, comptime role: []const u8) bool {
    return (core == .a64fx or core == .neoversev1) and role.len == 5 and (role[0] == 's' or role[0] == 'd');
}

fn usesThunderx(comptime core: Core, comptime role: []const u8) bool {
    return core == .thunderx and std.mem.eql(u8, role, "daxpy");
}

/// Whether `role` has a kernel of its own on `core`.
pub fn hasRole(comptime core: Core, comptime role: []const u8) bool {
    return R.hasRole(core, role) or usesSve(core, role) or usesThunderx(core, role);
}

/// Whether `core` has any axpy kernel of its own.
pub fn has(comptime core: Core) bool {
    inline for (roles) |role| {
        if (comptime hasRole(core, role)) return true;
    }
    return false;
}

/// Exports every axpy kernel `core` has. Call from a container-level `comptime`
/// block of the object's root file.
pub fn exportCore(comptime core: Core) void {
    inline for (roles) |role| {
        if (comptime R.hasRole(core, role)) {
            R.exportOne(role, core);
        } else if (comptime usesSve(core, role)) {
            const T = if (role[0] == 'd') f64 else f32;
            @export(&sve.Real(T).call, .{ .name = symbol(role, core), .visibility = .hidden });
        } else if (comptime usesThunderx(core, role)) {
            @export(uk.driven(kernel, role, thunderx_entry, uk.NoVector), .{ .name = symbol(role, core), .visibility = .hidden });
        }
    }
}

/// Points the vtable at the kernels `core` has.
pub fn install(vt: *Vtable, comptime core: Core) void {
    registry.installRoles(vt, core, roles, hasRole);
}
