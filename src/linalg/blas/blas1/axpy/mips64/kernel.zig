//! MIPS64 axpy kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`).
//!
//! The reference's MIPS64 axpy kernels that are assembly are for the Loongson
//! 3A (`KERNEL.LOONGSON3R3`): `axpy_loongson3a.S` for `saxpy`,
//! `daxpy_loongson3a_simd.S` for `daxpy` and the generic `zaxpy.S` for the
//! complex ones. They are in `asm/` as Zig files with their `#if`s and
//! macros resolved for a precision, and are emitted as naked functions whose
//! body is that assembly, reading their arguments straight from the ABI
//! registers as the reference's do.
//!
//! These have not been run: the sources prefetch by loading into `$0` and
//! `daxpy` encodes Loongson's paired loads and stores as `.word`s, which
//! neither qemu's MIPS CPUs nor other hardware accept, so this is here for
//! whoever has a Loongson 3A. The Loongson 3A R4 uses MSA C kernels where
//! the compiler has MSA, and the other cores use C, so they keep the portable
//! kernels.

const std = @import("std");
const coretypes = @import("../../../cpu/_coretypes.zig");
const asmexport = @import("../../../asm/exports.zig");
const axpy_loongson3a = @import("asm/axpy_loongson3a.zig");
const daxpy_loongson3a_simd = @import("asm/daxpy_loongson3a_simd.zig");
const zaxpy_s = @import("asm/zaxpy.zig");

pub const Core = coretypes.MipsCoreType;

const roles = [_][]const u8{ "saxpy", "daxpy", "caxpy", "caxpyc", "zaxpy", "zaxpyc" };

/// The source `role` uses on `core`.
fn source(comptime role: []const u8, comptime core: Core) ?type {
    if (core != .loongson3r3) return null;
    if (std.mem.eql(u8, role, "saxpy")) return axpy_loongson3a;
    if (std.mem.eql(u8, role, "daxpy")) return daxpy_loongson3a_simd;
    return zaxpy_s;
}

/// The kernel's assembly for `role` on `core`.
fn code(comptime role: []const u8, comptime core: Core) []const u8 {
    comptime {
        const double = role[0] == 'd' or role[0] == 'z';
        const conj = role[role.len - 1] == 'c';
        return source(role, core).?.text(double, conj, asmexport.labelSuffix(role, @tagName(core)));
    }
}

/// Whether `core` has any assembly axpy kernel.
pub fn has(comptime core: Core) bool {
    return comptime core == .loongson3r3;
}

/// Whether `role` has a kernel on `core`.
pub fn hasRole(comptime core: Core, comptime role: []const u8) bool {
    _ = role;
    return comptime core == .loongson3r3;
}

/// Exported symbol of `role` for `core`.
pub fn symbol(comptime role: []const u8, comptime core: Core) []const u8 {
    return asmexport.symbol(role, @tagName(core));
}

/// Exports every axpy kernel `core` has. Call from a container-level `comptime`
/// block of the object's root file.
pub fn exportCore(comptime core: Core) void {
    if (comptime core != .loongson3r3) return;
    inline for (roles) |role| {
        // `PROLOGUE`: `.align 5`, which on MIPS is 2^5 bytes
        @export(&asmexport.Naked(code(role, core), 32).f, .{ .name = symbol(role, core) });
    }
}
