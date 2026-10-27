//! MIPS axpy kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`), for both word sizes, and how they
//! get installed into the vtable (`install`): cores without an entry keep
//! the portable kernels from `axpy.installGeneric`. Per-core tests live in
//! `axpy/__root__.zig`, not here -- they need the full axpy interface
//! (`checkReal`, ...), which imports outside `src/linalg/blas/` and so
//! cannot appear in this file (see `kernels/__root__.zig`).
//!
//! What a core uses is its `KERNEL.<CORE>` choice, made in the reference by
//! `NO_MSA` where a core may or may not have the SIMD architecture: the SIMD
//! cores use `mips/{s,d,c,z}axpy_msa.c` (`msa.zig`), the Loongson 3A R3 the
//! whole-function `.S` kernels in `asm/` (`asm/__root__.zig`). The cores with
//! neither keep the portable kernels. On MIPS64 the architecture's `KERNEL`
//! then resets the complex ones to `arm/zaxpy.c` (it is read after the core's
//! file), so only `saxpy` and `daxpy` use the SIMD kernels there.

const sysinfo = @import("builtin");
const coretypes = @import("../../../cpu/_coretypes.zig");
const registry = @import("../../../asm_registry.zig");
const Vtable = @import("../../../Vtable.zig");
const kernel = @import("../kernel.zig");
const roles = kernel.roles;
const msa = @import("../mips/msa.zig");

const R = registry.Kernels(coretypes.MipsCoreType, @import("asm/__root__.zig"));

pub const Core = coretypes.MipsCoreType;

/// Whether `core` uses the SIMD kernel for `role`.
fn usesMsa(comptime core: Core, comptime role: []const u8) bool {
    const complex = comptime role[0] == 'c' or role[0] == 'z';
    if (complex and sysinfo.target.cpu.arch.isMips64()) return false;
    return switch (core) {
        .loongson3r4, .p5600, .i6400, .i6500, .p6600 => true,
        else => false,
    };
}

/// Whether `core` uses the SIMD kernel for any role.
fn usesAnyMsa(comptime core: Core) bool {
    inline for (roles) |r| if (usesMsa(core, r)) return true;
    return false;
}

/// Whether `role` has a kernel of its own on `core`.
pub fn hasRole(comptime core: Core, comptime role: []const u8) bool {
    return usesMsa(core, role) or R.hasRole(core, role);
}

/// Whether `core` has any axpy kernel of its own.
pub fn has(comptime core: Core) bool {
    return usesAnyMsa(core) or R.has(core);
}

/// Exported symbol of `role` for `core`.
pub const symbol = R.symbol;

/// Exports every axpy kernel `core` has. Call from a container-level `comptime`
/// block of the object's root file.
pub fn exportCore(comptime core: Core) void {
    if (comptime usesAnyMsa(core)) {
        inline for (.{
            .{ "saxpy", &msa.Real(f32).call },
            .{ "daxpy", &msa.Real(f64).call },
            .{ "caxpy", &msa.Cplx(f32, false).call },
            .{ "caxpyc", &msa.Cplx(f32, true).call },
            .{ "zaxpy", &msa.Cplx(f64, false).call },
            .{ "zaxpyc", &msa.Cplx(f64, true).call },
        }) |e| {
            if (comptime usesMsa(core, e[0])) @export(e[1], .{ .name = symbol(e[0], core), .visibility = .hidden });
        }
    }
    R.exportCore(core);
}

/// Points the vtable at the kernels `core` has.
pub fn install(vt: *Vtable, comptime core: Core) void {
    registry.installRoles(vt, core, roles, hasRole);
}
