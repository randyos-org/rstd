//! SPARC CPU selection: the reference's two static targets, by the build's
//! target (V9, or the generic one); nothing is detected at run time.

const std = @import("std");
const coretypes = @import("../coretypes.zig");
const select = @import("../select.zig");
const sysinfo = @import("builtin");

pub const CoreType = coretypes.SparcCoreType;
pub const Decision = select.Decision(CoreType);

/// Nothing is probed; the build's target is the answer.
pub const Caps = struct {
    v9: bool = sysinfo.target.cpu.arch.isSparc() and (sysinfo.target.cpu.arch == .sparc64 or sysinfo.target.cpu.has(.sparc, .v9)),
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

test "V9 targets get the V9 core" {
    try std.testing.expectEqual(CoreType.sparc, decide(.{ .v9 = true }).core);
    try std.testing.expectEqual(CoreType.sparcv7, decide(.{ .v9 = false }).core);
    try std.testing.expect(!runnable(.sparc, .{ .v9 = false }));
}

test "detect runs" {
    try std.testing.expect(select.canRun(@This(), select.detect(@This(), null)));
}
