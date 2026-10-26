//! Fallback for architectures with no detector of their own.
//!
//! The reference has the same notion: a `<ARCH>_GENERIC` kernel set built
//! from portable C, used whenever nothing more specific applies.

const coretypes = @import("../_coretypes.zig");
const select = @import("../_select.zig");

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

pub fn detect() CoreType {
    return .generic;
}

pub fn detectNoted() Decision {
    return decide(probe());
}

const std = @import("std");

test "generic detect is always available" {
    try std.testing.expectEqual(CoreType.generic, detect());
}
