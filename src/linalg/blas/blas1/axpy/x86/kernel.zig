//! 32-bit x86 axpy kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`).
//!
//! The reference's 32-bit x86 axpy kernels are the whole-function `.S`
//! `axpy.S` and `zaxpy.S`, which work on the x87 stack; no core has anything
//! else (`KERNEL.<CORE>` in `kernel/x86/` names none, so the default of
//! `Makefile.L1` applies). They are in `asm/` as Zig files with their
//! `#if`s and macros resolved for a precision, and are emitted as naked
//! functions whose body is that assembly, reading their arguments off the
//! stack as the reference's do. The x87 `qaxpy`/`xaxpy` are not ported: they
//! only exist upstream in `EXPRECISION` builds.

const std = @import("std");
const coretypes = @import("../../../cpu/_coretypes.zig");
const asmexport = @import("../../../asm/exports.zig");
const axpy_x87 = @import("asm/axpy.zig");
const zaxpy_x87 = @import("asm/zaxpy.zig");

pub const Core = coretypes.X86CoreType;

const roles = [_][]const u8{ "saxpy", "daxpy", "caxpy", "caxpyc", "zaxpy", "zaxpyc" };

/// `HAVE_3DNOW` of `getarch.c`: the cores that prefetch with `prefetch` and
/// `prefetchw`.
fn have3dnow(comptime core: Core) bool {
    return switch (core) {
        .athlon, .opteron, .opteron_sse3 => true,
        else => false,
    };
}

/// The kernel's assembly for `role` on `core`.
fn code(comptime role: []const u8, comptime core: Core) []const u8 {
    comptime {
        const sfx = asmexport.labelSuffix(role, @tagName(core));
        const real = std.mem.eql(u8, role, "saxpy") or std.mem.eql(u8, role, "daxpy");
        const double = role[0] == 'd' or role[0] == 'z';
        const conj = role[role.len - 1] == 'c';
        return if (real)
            axpy_x87.text(double, conj, have3dnow(core), sfx)
        else
            zaxpy_x87.text(double, conj, have3dnow(core), sfx);
    }
}

/// Whether `core` has an assembly axpy kernel: all but `unknown`, which
/// keeps the portable ones.
pub fn has(comptime core: Core) bool {
    return core != .unknown;
}

/// Whether `role` has a kernel on `core`.
pub fn hasRole(comptime core: Core, comptime role: []const u8) bool {
    _ = role;
    return has(core);
}

/// Exported symbol of `role` for `core`.
pub fn symbol(comptime role: []const u8, comptime core: Core) []const u8 {
    return asmexport.symbol(role, @tagName(core));
}

/// Exports every axpy kernel `core` has. Call from a container-level `comptime`
/// block of the object's root file.
pub fn exportCore(comptime core: Core) void {
    if (comptime !has(core)) return;
    inline for (roles) |role| {
        // `PROLOGUE`: `.align 16`
        @export(&asmexport.Naked(code(role, core), 16).f, .{ .name = symbol(role, core) });
    }
}
