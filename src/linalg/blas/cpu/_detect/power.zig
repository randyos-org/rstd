//! PowerPC CPU detection, mirroring `driver/others/dynamic_power.c`.
//!
//! The reference uses GCC's `__builtin_cpu_is`, which glibc backs with
//! `AT_PLATFORM` from the auxiliary vector, and falls back to reading the
//! PVR register directly on AIX. We read `AT_PLATFORM` ourselves, since
//! Zig has no equivalent sysinfo.

const std = @import("std");
const coretypes = @import("../_coretypes.zig");
const sysinfo = @import("builtin");
const auxv = @import("../_auxv.zig");

pub const CoreType = coretypes.PowerCoreType;

/// The `AT_PLATFORM` names `detect` distinguishes; tags match the strings.
pub const Platform = enum { power11, power10, power9, power8, power7, power6, power6x };

/// `AT_PLATFORM` (a NUL-terminated string such as "power9"), if it is one we
/// know. A chip in compatibility mode reports the mode it is in.
pub fn platform() ?Platform {
    const ptr = auxv.get(auxv.AT_PLATFORM);
    if (ptr == 0) return null;
    return std.meta.stringToEnum(Platform, std.mem.sliceTo(@as([*:0]const u8, @ptrFromInt(ptr)), 0));
}

pub fn detect() CoreType {
    return switch (platform() orelse return .unknown) {
        .power11, .power10 => .power10,
        .power9 => .power9,
        .power8 => .power8,
        .power7, .power6, .power6x => .power6,
    };
}

test "detect runs on PowerPC targets" {
    if (comptime !sysinfo.cpu.arch.isPowerPC()) return error.SkipZigTest;
    _ = detect();
}
