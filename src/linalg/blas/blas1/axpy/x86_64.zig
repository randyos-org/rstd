//! x86_64 side of axpy: picks the per-core kernels for a core type.
//!
//! The kernels themselves (`kernel_x86_64.zig`) are compiled per core into
//! objects with that core's CPU features, by the build (`build/blaskern.zig`),
//! and linked in; this file only names their symbols and
//! points the vtable at them, following the reference's `KERNEL.<CORE>`
//! choice. Cores without an entry keep the portable kernels from
//! `axpy.installGeneric`.

const std = @import("std");
const sysinfo = @import("builtin");
const zigstd = @import("../../../../zigstd.zig");
const Vtable = @import("../../Vtable.zig");
const x86 = @import("../../cpu/x86.zig");
const axpy = @import("__root__.zig");
const kx = @import("kernel_x86_64.zig");

fn ext(comptime K: type, comptime name: []const u8) K {
    return @extern(K, .{ .name = name });
}

fn installCore(vt: *Vtable, comptime core: x86.CoreType) void {
    if (comptime kx.has(core, "daxpy")) vt.daxpy_k = ext(axpy.RealKernel(f64), kx.symbol("daxpy", core));
    if (comptime kx.has(core, "saxpy")) vt.saxpy_k = ext(axpy.RealKernel(f32), kx.symbol("saxpy", core));
    if (comptime kx.has(core, "caxpy")) {
        vt.caxpy_k = ext(axpy.ComplexKernel(f32), kx.symbol("caxpy", core));
        vt.caxpyc_k = ext(axpy.ComplexKernel(f32), kx.symbol("caxpyc", core));
    }
    if (comptime kx.has(core, "zaxpy")) {
        vt.zaxpy_k = ext(axpy.ComplexKernel(f64), kx.symbol("zaxpy", core));
        vt.zaxpyc_k = ext(axpy.ComplexKernel(f64), kx.symbol("zaxpyc", core));
    }
}

/// Overrides the generic kernels installed by `axpy.install` for `core`.
pub fn install(vt: *Vtable, core: x86.CoreType) void {
    switch (core) {
        inline .prescott,
        .core2,
        .penryn,
        .dunnington,
        .opteron,
        .opteron_sse3,
        .barcelona,
        .bobcat,
        .nano,
        .atom,
        .nehalem,
        .sandybridge,
        .bulldozer,
        .piledriver,
        .steamroller,
        .excavator,
        .haswell,
        .zen,
        .skylakex,
        .cooperlake,
        .sapphirerapids,
        => |c| installCore(vt, c),
        else => {},
    }
}

// ---- tests -----------------------------------------------------------------
// Each core's kernels run only where the host CPU can execute them.

test "per-core kernels match a scalar loop where the host supports them" {
    if (comptime sysinfo.cpu.arch != .x86_64) return error.SkipZigTest;
    inline for (@typeInfo(x86.CoreType).@"enum".field_names) |name| {
        const core = @field(x86.CoreType, name);
        var vt: Vtable = .{};
        axpy.installGeneric(&vt);
        install(&vt, core);
        if (x86.supportsCore(core)) {
            try axpy.checkReal(f32, vt.saxpy_k.?);
            try axpy.checkReal(f64, vt.daxpy_k.?);
            try axpy.checkComplex(f32, false, vt.caxpy_k.?);
            try axpy.checkComplex(f32, true, vt.caxpyc_k.?);
            try axpy.checkComplex(f64, false, vt.zaxpy_k.?);
            try axpy.checkComplex(f64, true, vt.zaxpyc_k.?);
        }
    }
}
