//! SPARC CPU selection. The reference has two static targets and detects
//! nothing at run time, so the core type follows the target: V9 (any 64-bit
//! target, or the 32-bit one built for V9) or the generic one.

const std = @import("std");
const coretypes = @import("../_coretypes.zig");
const select = @import("../_select.zig");
const sysinfo = @import("builtin");

pub const CoreType = coretypes.SparcCoreType;
pub const Decision = select.Decision(CoreType);

/// Nothing is probed; the build's target is the answer.
pub const Caps = struct {
    v9: bool = sysinfo.cpu.arch.isSPARC() and (sysinfo.cpu.arch == .sparc64 or sysinfo.cpu.has(.sparc, .v9)),
};

pub fn probe() Caps {
    return .{};
}

pub fn runnable(core: CoreType, c: Caps) bool {
    return core == .sparcv7 or c.v9;
}

pub fn decide(c: Caps) Decision {
    return .{ .core = if (c.v9) .sparc else .sparcv7 };
}

/// Whether the kernels of `core` run on this CPU.
pub fn supportsCore(core: CoreType) bool {
    return runnable(core, probe());
}

pub fn detect() CoreType {
    return decide(probe()).core;
}

pub fn detectNoted() Decision {
    return decide(probe());
}

test "V9 targets get the V9 core" {
    try std.testing.expectEqual(CoreType.sparc, decide(.{ .v9 = true }).core);
    try std.testing.expectEqual(CoreType.sparcv7, decide(.{ .v9 = false }).core);
    try std.testing.expect(!runnable(.sparc, .{ .v9 = false }));
}

test "detect runs" {
    try std.testing.expect(supportsCore(detect()));
}
