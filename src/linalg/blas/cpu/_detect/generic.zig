//! Fallback for architectures with no detector of their own.
//!
//! The reference has the same notion: a `<ARCH>_GENERIC` kernel set built
//! from portable C, used whenever nothing more specific applies.

const std = @import("std");
const coretypes = @import("../_coretypes.zig");
const zigstd = @import("../../../../zigstd.zig");

pub const CoreType = coretypes.GenericCoreType;

pub fn detect() CoreType {
    return .generic;
}

test "generic detect is always available" {
    try zigstd.testExpectEqual(CoreType.generic, detect());
}
