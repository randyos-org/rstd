//! s390x axpy kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`), and how they get installed into
//! the vtable (`install`): cores without an entry keep the portable kernels
//! from `axpy.installGeneric`. Per-core tests live in `axpy/__root__.zig`,
//! not here -- they need the full axpy interface (`checkReal`, ...), which
//! imports outside `src/linalg/blas/` and so cannot appear in this file (see
//! `kernels/__root__.zig`).
//!
//! The reference's `zarch/{s,d,c,z}axpy.c` are the shared C driver around a
//! block that is inline assembly on the vector facility: the block is in
//! `uk/`, as written, and composed with the driver in `../kernel.zig`.
//! `uk/__root__.zig` says which core uses which. The cores whose
//! `KERNEL.<CORE>` names the plain C loops keep the portable kernels.

const coretypes = @import("../../../cpu/_coretypes.zig");
const util = @import("../../../util.zig");
const Vtable = @import("../../../Vtable.zig");
const kernel = @import("../kernel.zig");
const roles = kernel.roles;
const uk = @import("../uk_block.zig");
const ukidx = @import("uk/__root__.zig");

pub const Core = coretypes.S390xCoreType;

/// The entry for `role` on `core`, if it has one.
fn entry(comptime role: []const u8, comptime core: Core) ?ukidx.Entry {
    return ukidx.find(role, @tagName(core));
}

/// Whether `role` has a kernel of its own on `core`.
pub fn hasRole(comptime core: Core, comptime role: []const u8) bool {
    return entry(role, core) != null;
}

/// Whether `core` has any axpy kernel of its own.
pub fn has(comptime core: Core) bool {
    inline for (roles) |role| {
        if (comptime hasRole(core, role)) return true;
    }
    return false;
}

/// Exported symbol of `role` for `core`.
pub fn symbol(comptime role: []const u8, comptime core: Core) []const u8 {
    return util.kernelSymbol(role, @tagName(core));
}

/// Exports every axpy kernel `core` has. Call from a container-level `comptime`
/// block of the object's root file.
pub fn exportCore(comptime core: Core) void {
    inline for (roles) |role| {
        if (comptime entry(role, core)) |e| {
            @export(uk.driven(role, e, uk.NoVector), .{ .name = symbol(role, core) });
        }
    }
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
