//! PowerPC CPU detection, mirroring `driver/others/dynamic_power.c`.
//!
//! The reference uses GCC's `__builtin_cpu_is`, which glibc backs with
//! `AT_PLATFORM` from the auxiliary vector, and falls back to reading the
//! PVR register directly on AIX. We read `AT_PLATFORM` ourselves, since
//! Zig has no equivalent sysinfo.

const std = @import("std");
const coretypes = @import("_coretypes.zig");
const sysinfo = @import("builtin");
const auxv = @import("_auxv.zig");

pub const CoreType = coretypes.PowerCoreType;
pub const kernel = @import("../kernels/power.zig");

/// `AT_PLATFORM` points at a NUL-terminated string such as "power9".
pub fn platform() ?[]const u8 {
    const ptr = auxv.get(auxv.AT_PLATFORM);
    if (ptr == 0) return null;
    return std.mem.sliceTo(@as([*:0]const u8, @ptrFromInt(ptr)), 0);
}

pub fn detect() CoreType {
    const name = platform() orelse return .unknown;

    // Newest first: a "power11" chip running in compatibility mode still
    // reports the mode it is in, so an exact match is what we want.
    if (std.mem.eql(u8, name, "power11")) return .power10;
    if (std.mem.eql(u8, name, "power10")) return .power10;
    if (std.mem.eql(u8, name, "power9")) return .power9;
    if (std.mem.eql(u8, name, "power8")) return .power8;
    if (std.mem.eql(u8, name, "power7")) return .power6;
    if (std.mem.eql(u8, name, "power6") or std.mem.eql(u8, name, "power6x")) return .power6;
    return .unknown;
}

const Vtable = @import("../Vtable.zig");
const axpy = @import("../blas1/axpy/__root__.zig");
pub fn vtableFromCoreType(core_type: CoreType) Vtable {
    var vt: Vtable = .{};
    axpy.install(&vt, core_type);
    return vt;
}

test "detect runs on PowerPC targets" {
    if (comptime !sysinfo.cpu.arch.isPowerPC()) return error.SkipZigTest;
    _ = detect();
}
