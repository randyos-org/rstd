//! POWER axpy kernels, one set per core, built into the per-core kernel objects
//! (`src/linalg/blas/_kernel.zig`), and how they get installed into the
//! vtable (`install`): cores without an entry keep the portable kernels from
//! `axpy.installGeneric`. Per-core tests live in `axpy/__root__.zig`, not
//! here -- they need the full axpy interface (`checkReal`, ...), which
//! imports outside `src/linalg/blas/` and so cannot appear in this file (see
//! `kernels/__root__.zig`).
//!
//! The kernels that are assembly upstream are in `asm/`, one Zig file per
//! source and precision with its `#if`s and macros resolved (`asm/__root__.zig`
//! says which core uses which). `POWER8` and up use the C kernels of
//! `c_kernels.zig`.

const std = @import("std");
const sysinfo = @import("builtin");
const coretypes = @import("../../../cpu/_coretypes.zig");
const registry = @import("../../../asm_registry.zig");
const Vtable = @import("../../../Vtable.zig");
const kernel = @import("../kernel.zig");
const roles = kernel.roles;
const c_kernels = @import("c_kernels.zig");

const R = registry.Kernels(coretypes.PowerCoreType, @import("asm/__root__.zig"));

pub const Core = coretypes.PowerCoreType;
pub const symbol = R.symbol;

/// Which C kernels `core` uses (`KERNEL.POWER8`, `.POWER9`, `.POWER10`); the
/// assembly `caxpy` alternatives of the first two are for older compilers.
fn cSource(comptime core: Core) ?c_kernels.Source {
    return switch (core) {
        .power8, .power9 => .power8,
        .power10 => .power10,
        else => null,
    };
}

/// Exports every axpy kernel `core` has.
pub fn exportCore(comptime core: Core) void {
    inline for (roles) |role| {
        if (comptime cSource(core)) |source| {
            @export(c_kernels.entry(role, source), .{ .name = symbol(role, core), .visibility = .hidden });
        } else if (comptime R.hasRole(core, role)) {
            R.exportOne(role, core);
        }
    }
}

/// Whether `core` has any axpy kernel. (The 32-bit cores have no 64-bit
/// code; `kernels/power.zig`'s `runsOn` keeps them off 64-bit targets.)
pub fn has(comptime core: Core) bool {
    return cSource(core) != null or R.has(core);
}

/// Whether `role` has a kernel on `core` for this target.
pub fn hasRole(comptime core: Core, comptime role: []const u8) bool {
    return cSource(core) != null or R.hasRole(core, role);
}

/// Points the vtable at the kernels `core` has.
pub fn install(vt: *Vtable, comptime core: Core) void {
    registry.installRoles(vt, core, roles, hasRole);
}
