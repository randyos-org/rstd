//! AArch64 axpy kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`).
//!
//! All of the reference's AArch64 axpy kernels that are assembly are here, as
//! the whole-function `.S` kernels in `asm/`, one Zig file per source
//! with its `#if`s and macros resolved for a precision. They are emitted as
//! naked functions whose body is that assembly and read their arguments
//! straight from the ABI registers, as the reference's do. The cores whose
//! reference kernel is C (`axpy_sve.c` for Neoverse V1 and A64FX,
//! `daxpy_thunderx.c` for ThunderX) are not ported and keep the portable
//! kernels.
//!
//! Which source a core uses is its `KERNEL.<CORE>` choice in
//! `kernel/arm64/`.

const std = @import("std");
const coretypes = @import("../../../cpu/_coretypes.zig");
const asmexport = @import("../../../asm/exports.zig");
const axpy_s = @import("asm/axpy.zig");
const zaxpy_s = @import("asm/zaxpy.zig");
const daxpy_thunderx2t99 = @import("asm/daxpy_thunderx2t99.zig");

pub const Core = coretypes.Arm64CoreType;

const roles = [_][]const u8{ "saxpy", "daxpy", "caxpy", "caxpyc", "zaxpy", "zaxpyc" };

/// The source `role` uses on `core`, if it is assembly there.
fn source(comptime role: []const u8, comptime core: Core) ?type {
    const real = std.mem.eql(u8, role, "saxpy") or std.mem.eql(u8, role, "daxpy");
    if (!real) return if (core == .unknown) null else zaxpy_s;
    if (std.mem.eql(u8, role, "saxpy")) return switch (core) {
        .neoversev1, .a64fx, .unknown => null,
        else => axpy_s,
    };
    return switch (core) {
        .armv8, .cortexa53, .cortexa55, .cortexa57, .cortexa72, .cortexa73, .falkor, .tsv110, .emag8180 => axpy_s,
        .thunderx2t99, .thunderx3t110, .neoversen1, .neoversen2, .neoversev2, .armv8sve, .armv9sme, .vortex, .vortexm4 => daxpy_thunderx2t99,
        .thunderx, .neoversev1, .a64fx, .unknown => null,
    };
}

/// The kernel's assembly for `role` on `core`.
fn code(comptime role: []const u8, comptime core: Core) []const u8 {
    comptime {
        const sfx = asmexport.labelSuffix(role, @tagName(core));
        const conj = role[role.len - 1] == 'c';
        const double = role[0] == 'd' or role[0] == 'z';
        const file = source(role, core).?;
        return if (file == daxpy_thunderx2t99) file.text(true, false, sfx) else file.text(double, conj, sfx);
    }
}

/// Whether `core` has any assembly axpy kernel.
pub fn has(comptime core: Core) bool {
    inline for (roles) |role| {
        if (comptime source(role, core) != null) return true;
    }
    return false;
}

/// Exported symbol of `role` for `core`.
pub fn symbol(comptime role: []const u8, comptime core: Core) []const u8 {
    return asmexport.symbol(role, @tagName(core));
}

/// Whether `role` has a kernel on `core`.
pub fn hasRole(comptime core: Core, comptime role: []const u8) bool {
    return comptime source(role, core) != null;
}

/// Exports every axpy kernel `core` has. Call from a container-level `comptime`
/// block of the object's root file.
pub fn exportCore(comptime core: Core) void {
    inline for (roles) |role| {
        if (comptime source(role, core) != null) {
            // `PROLOGUE`: `.p2align 2`
            @export(&asmexport.Naked(code(role, core), 4).f, .{ .name = symbol(role, core) });
        }
    }
}
