//! LoongArch axpy kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`).
//!
//! The reference's LoongArch axpy kernels are all assembly, one source per
//! vector extension (LSX and LASX) and per data kind, in `asm/` as
//! Zig files with their `#if`s and macros resolved for a precision. They are
//! emitted as naked functions whose body is that assembly and read their
//! arguments straight from the ABI registers, as the reference's do. Cores
//! with neither extension keep the portable kernels.
//!
//! Which source a core uses is its `KERNEL.<CORE>` choice in
//! `kernel/loongarch64/`: LA264 uses the LSX ones, LA464 the LASX ones. The
//! other cores with the extensions (`corespec.loongarch`) follow.

const std = @import("std");
const coretypes = @import("../../../cpu/_coretypes.zig");
const asmexport = @import("../../../asm/exports.zig");
const axpy_lsx = @import("asm/axpy_lsx.zig");
const axpy_lasx = @import("asm/axpy_lasx.zig");
const caxpy_lsx = @import("asm/caxpy_lsx.zig");
const caxpy_lasx = @import("asm/caxpy_lasx.zig");

pub const Core = coretypes.LoongarchCoreType;

const roles = [_][]const u8{ "saxpy", "daxpy", "caxpy", "caxpyc", "zaxpy", "zaxpyc" };

const Ext = enum { none, lsx, lasx };

fn ext(comptime core: Core) Ext {
    return switch (core) {
        .la264, .loongson2k1000 => .lsx,
        .la464, .loongson3r5 => .lasx,
        .la64_generic, .loongsongeneric, .unknown => .none,
    };
}

/// The source `role` uses on `core`.
fn source(comptime role: []const u8, comptime core: Core) ?type {
    const real = std.mem.eql(u8, role, "saxpy") or std.mem.eql(u8, role, "daxpy");
    return switch (ext(core)) {
        .none => null,
        .lsx => if (real) axpy_lsx else caxpy_lsx,
        .lasx => if (real) axpy_lasx else caxpy_lasx,
    };
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
    return comptime ext(core) != .none;
}

/// Whether `role` has a kernel on `core`.
pub fn hasRole(comptime core: Core, comptime role: []const u8) bool {
    _ = role;
    return comptime ext(core) != .none;
}

/// Exported symbol of `role` for `core`.
pub fn symbol(comptime role: []const u8, comptime core: Core) []const u8 {
    return asmexport.symbol(role, @tagName(core));
}

/// Exports every axpy kernel `core` has. Call from a container-level `comptime`
/// block of the object's root file.
pub fn exportCore(comptime core: Core) void {
    if (comptime ext(core) == .none) return;
    inline for (roles) |role| {
        // `PROLOGUE`: `.align 5`, which on LoongArch is 2^5 bytes
        @export(&asmexport.Naked(code(role, core), 32).f, .{ .name = symbol(role, core) });
    }
}
