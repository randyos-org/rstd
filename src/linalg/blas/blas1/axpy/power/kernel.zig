//! POWER axpy kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`).
//!
//! The reference's assembly axpy kernels for POWER: `axpy.S` and `zaxpy.S`
//! (the default for the older POWER cores and the 970), `axpy_ppc440.S` and
//! `zaxpy_ppc440.S` (PPC440 and G4), all in `asm/` as Zig files with
//! their `#if`s and macros resolved for a precision, core and ABI. They are
//! emitted as naked functions whose body is that assembly and read their
//! arguments straight from the ABI registers, as the reference's do.
//!
//! The PPC440FP2's kernels (`*_hummer.S`) use the Double Hummer paired
//! floating-point instructions, which LLVM cannot assemble; they are
//! generated but not exported. Otherwise the POWER8, POWER9 and POWER10
//! kernels are C, except `caxpy` on little-endian POWER8 and POWER9, which is
//! GCC output kept as assembly.

const std = @import("std");
const coretypes = @import("../../../cpu/_coretypes.zig");
const asmexport = @import("../../../asm/exports.zig");
const axpy_s = @import("asm/axpy.zig");
const zaxpy_s = @import("asm/zaxpy.zig");
const axpy_ppc440 = @import("asm/axpy_ppc440.zig");
const zaxpy_ppc440 = @import("asm/zaxpy_ppc440.zig");
const caxpy_power8 = @import("asm/caxpy_power8.zig");
const caxpy_power9 = @import("asm/caxpy_power9.zig");

pub const Core = coretypes.PowerCoreType;

const roles = [_][]const u8{ "saxpy", "daxpy", "caxpy", "caxpyc", "zaxpy", "zaxpyc" };

/// The 32-bit cores only exist as 32-bit targets.
fn existsOn(comptime arch: std.Target.Cpu.Arch, comptime core: Core) bool {
    return switch (core) {
        .ppc440, .ppc440fp2, .ppcg4 => !arch.isPowerPC64(),
        else => true,
    };
}

/// The source `role` uses on `core`: its `KERNEL.<CORE>` choice, else the
/// default of `Makefile.L1`.
fn source(comptime arch: std.Target.Cpu.Arch, comptime role: []const u8, comptime core: Core) ?type {
    if (!existsOn(arch, core)) return null;
    const real = std.mem.eql(u8, role, "saxpy") or std.mem.eql(u8, role, "daxpy");
    return switch (core) {
        .power3, .power4, .power5, .power6, .ppc970, .cell => if (real) axpy_s else zaxpy_s,
        // the complex ones are `../arm/zaxpy.c` on big-endian
        .ppc440, .ppcg4 => if (real) axpy_ppc440 else if (arch.endian() == .little) zaxpy_ppc440 else null,
        // `caxpy_power8.S` and `caxpy_power9.S` are what the reference uses
        // where the compiler is not GCC 9 or later, which includes clang; on
        // big-endian it is `caxpy.c` (C) and the double complex is C too
        .power8, .power9 => if (std.mem.eql(u8, role, "caxpy") or std.mem.eql(u8, role, "caxpyc"))
            (if (arch.endian() != .little) null else if (core == .power8) caxpy_power8 else caxpy_power9)
        else
            null,
        .power10, .ppc440fp2, .unknown => null,
    };
}

/// The kernel's assembly for `role` on `core`.
fn code(comptime arch: std.Target.Cpu.Arch, comptime role: []const u8, comptime core: Core) []const u8 {
    comptime {
        const double = role[0] == 'd' or role[0] == 'z';
        const conj = role[role.len - 1] == 'c';
        const sfx = asmexport.labelSuffix(role, @tagName(core));
        return source(arch, role, core).?.text(double, conj, arch, core, symbol(role, core), sfx);
    }
}

/// Whether `core` has any assembly axpy kernel on `arch`.
pub fn has(comptime arch: std.Target.Cpu.Arch, comptime core: Core) bool {
    inline for (roles) |role| {
        if (comptime source(arch, role, core) != null) return true;
    }
    return false;
}

/// Whether `role` has a kernel on `core`.
pub fn hasRole(comptime arch: std.Target.Cpu.Arch, comptime core: Core, comptime role: []const u8) bool {
    return comptime source(arch, role, core) != null;
}

/// Exported symbol of `role` for `core`.
pub fn symbol(comptime role: []const u8, comptime core: Core) []const u8 {
    return asmexport.symbol(role, @tagName(core));
}

/// Exports every axpy kernel `core` has. Call from a container-level `comptime`
/// block of the object's root file.
pub fn exportCore(comptime arch: std.Target.Cpu.Arch, comptime core: Core) void {
    inline for (roles) |role| {
        if (comptime source(arch, role, core) != null) {
            // `PROLOGUE`: `.align 6`, which on POWER is 2^6 bytes
            @export(&asmexport.Naked(code(arch, role, core), 64).f, .{ .name = symbol(role, core) });
        }
    }
}
