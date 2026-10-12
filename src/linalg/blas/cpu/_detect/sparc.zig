//! SPARC CPU selection. The reference has two static targets and detects
//! nothing at run time, so the core type follows the target: V9 (any 64-bit
//! target, or the 32-bit one built for V9) or the generic one.

const std = @import("std");
const coretypes = @import("../_coretypes.zig");
const sysinfo = @import("builtin");

pub const CoreType = coretypes.SparcCoreType;

pub fn detect() CoreType {
    if (comptime sysinfo.cpu.arch.isSPARC() and (sysinfo.cpu.arch == .sparc64 or sysinfo.cpu.has(.sparc, .v9))) return .sparc;
    return .sparcv7;
}

test "detect runs" {
    _ = detect();
}
