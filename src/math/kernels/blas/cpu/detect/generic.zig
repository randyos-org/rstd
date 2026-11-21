//! Fallback for architectures with no detector of their own, as the
//! reference's portable `<ARCH>_GENERIC` kernel set.

const std = @import("std");
const coretypes = @import("../coretypes.zig");
const select = @import("../select.zig");

pub const CoreType = coretypes.GenericCoreType;
pub const Decision = select.Decision(CoreType);

pub const Caps = struct {};

pub fn probe() Caps {
    return .{};
}

pub fn runnable(core: CoreType, c: Caps) bool {
    _ = core;
    _ = c;
    return true;
}

pub fn decide(c: Caps) Decision {
    _ = c;
    return .{ .core = .generic };
}

test "generic detect is always available" {
    try std.testing.expectEqual(CoreType.generic, select.detect(@This(), null));
}
