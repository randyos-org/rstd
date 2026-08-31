//! std lib aliases for convenience

const std = @import("std");

pub const Allocator = std.mem.Allocator;
pub const MainInit = std.process.Init;
pub const StdIo = std.Io;

pub const equal = std.mem.eql;
pub const sprintf = std.mem.print;

pub const testExpectTrue = std.testing.expect;
pub const testExpectEqual = std.testing.expectEqual;
pub const testExpectEqualSlices = std.testing.expectEqualSlices;
pub const testExpectError = std.testing.expectError;
pub const testAllocator = std.testing.allocator;
