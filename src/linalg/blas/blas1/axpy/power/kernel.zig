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
const kregistry = @import("../../../kernels/power.zig");
const kernel = @import("../kernel.zig");
const roles = kernel.roles;
const c_kernels = @import("c_kernels.zig");

/// The sources are built once per pointer width, as their register choices
/// differ.
const width = if (@sizeOf(usize) == 8) "64" else "32";

const R = registry.Registry(coretypes.PowerCoreType, @import("asm/__root__.zig"), width, roles);

pub const Core = R.Core;
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

/// Exports every axpy kernel `core` has, for `arch`.
pub fn exportCore(comptime arch: std.Target.Cpu.Arch, comptime core: Core) void {
    _ = arch;
    inline for (roles) |role| {
        if (comptime cSource(core)) |source| {
            @export(c_kernels.entry(role, source), .{ .name = symbol(role, core) });
        } else if (comptime R.hasRole(core, role)) {
            R.exportOne(role, core);
        }
    }
}

/// Whether `core` has any axpy kernel on `arch`.
pub fn has(comptime arch: std.Target.Cpu.Arch, comptime core: Core) bool {
    return exists(arch, core) and (cSource(core) != null or R.hasAny(core));
}

/// The 32-bit cores have no 64-bit code, so no kernels there.
fn exists(comptime arch: std.Target.Cpu.Arch, comptime core: Core) bool {
    return switch (core) {
        .ppc440, .ppc440fp2, .ppcg4 => !arch.isPowerPC64(),
        else => true,
    };
}

/// Whether `role` has a kernel on `core` for `arch`. Both byte orders and the
/// 32-bit ABIs share the files.
pub fn hasRole(comptime arch: std.Target.Cpu.Arch, comptime core: Core, comptime role: []const u8) bool {
    return exists(arch, core) and (cSource(core) != null or R.hasRole(core, role));
}

fn ext(comptime K: type, comptime name: []const u8) K {
    return @extern(K, .{ .name = name });
}

fn installCore(vt: *Vtable, comptime core: Core) void {
    // no object is linked for a core that does not run on this CPU
    if (comptime !kregistry.runsOn(sysinfo.cpu, core)) return;
    if (comptime hasRole(sysinfo.cpu.arch, core, "daxpy")) vt.daxpy_k = ext(kernel.RealKernel(f64), symbol("daxpy", core));
    if (comptime hasRole(sysinfo.cpu.arch, core, "saxpy")) vt.saxpy_k = ext(kernel.RealKernel(f32), symbol("saxpy", core));
    if (comptime hasRole(sysinfo.cpu.arch, core, "caxpy")) {
        vt.caxpy_k = ext(kernel.ComplexKernel(f32), symbol("caxpy", core));
        vt.caxpyc_k = ext(kernel.ComplexKernel(f32), symbol("caxpyc", core));
    }
    if (comptime hasRole(sysinfo.cpu.arch, core, "zaxpy")) {
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
