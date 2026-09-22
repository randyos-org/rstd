//! POWER side of axpy: picks the per-core kernels for a core type.
//!
//! The kernels themselves (`kernel.zig`) are compiled per core into
//! objects by the build (`build/blaskern.zig`) and linked in; this file only
//! names their symbols and points the vtable at them, following the
//! reference's `KERNEL.<CORE>` choice. Cores without an entry keep the
//! portable kernels from `axpy.installGeneric`.

const std = @import("std");
const sysinfo = @import("builtin");
const Vtable = @import("../../../Vtable.zig");
const power = @import("../../../cpu/power.zig");
const axpy = @import("../__root__.zig");
const kx = @import("kernel.zig");
const registry = @import("../../../kernels/power.zig");

fn ext(comptime K: type, comptime name: []const u8) K {
    return @extern(K, .{ .name = name });
}

fn installCore(vt: *Vtable, comptime core: power.CoreType) void {
    // no object is linked for a core that does not run on this CPU
    if (comptime !registry.runsOn(sysinfo.cpu, core)) return;
    if (comptime kx.hasRole(sysinfo.cpu.arch, core, "daxpy")) vt.daxpy_k = ext(axpy.RealKernel(f64), kx.symbol("daxpy", core));
    if (comptime kx.hasRole(sysinfo.cpu.arch, core, "saxpy")) vt.saxpy_k = ext(axpy.RealKernel(f32), kx.symbol("saxpy", core));
    if (comptime kx.hasRole(sysinfo.cpu.arch, core, "caxpy")) {
        vt.caxpy_k = ext(axpy.ComplexKernel(f32), kx.symbol("caxpy", core));
        vt.caxpyc_k = ext(axpy.ComplexKernel(f32), kx.symbol("caxpyc", core));
    }
    if (comptime kx.hasRole(sysinfo.cpu.arch, core, "zaxpy")) {
        vt.zaxpy_k = ext(axpy.ComplexKernel(f64), kx.symbol("zaxpy", core));
        vt.zaxpyc_k = ext(axpy.ComplexKernel(f64), kx.symbol("zaxpyc", core));
    }
}

/// Overrides the generic kernels installed by `axpy.install` for `core`.
pub fn install(vt: *Vtable, core: power.CoreType) void {
    switch (core) {
        inline else => |c| installCore(vt, c),
    }
}

// ---- tests -----------------------------------------------------------------

test "per-core kernels match a scalar loop" {
    inline for (@typeInfo(power.CoreType).@"enum".field_names) |name| {
        const core = @field(power.CoreType, name);
        var vt: Vtable = .{};
        axpy.installGeneric(&vt);
        install(&vt, core);
        try axpy.checkReal(f32, vt.saxpy_k.?);
        try axpy.checkReal(f64, vt.daxpy_k.?);
        try axpy.checkComplex(f32, false, vt.caxpy_k.?);
        try axpy.checkComplex(f32, true, vt.caxpyc_k.?);
        try axpy.checkComplex(f64, false, vt.zaxpy_k.?);
        try axpy.checkComplex(f64, true, vt.zaxpyc_k.?);
    }
}
