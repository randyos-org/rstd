//! std lib aliases for convenience

const std = @import("std");
const zigstd = @import("zigstd.zig");

pub const Allocator = std.mem.Allocator;
pub const MainInit = std.process.Init;
pub const StdIo = std.Io;
pub const SrcInfo = std.lang.SourceLocation;

pub const equal = std.mem.eql;
pub const sprintf = std.mem.print;
pub const log2Int = std.math.log2_int;
pub const assert = std.debug.assert;

pub const testExpectTrue = std.testing.expect;
pub const testExpectEqual = std.testing.expectEqual;
pub const testExpectEqualSlices = std.testing.expectEqualSlices;
pub const testExpectEqualStrings = std.testing.expectEqualStrings;
pub const testExpectError = std.testing.expectError;
pub const testExpectApproxEqAbs = std.testing.expectApproxEqAbs;
pub const testExpectApproxEqRel = std.testing.expectApproxEqRel;
pub const testRefAllDecls = std.testing.refAllDecls;
pub const testAllocator = std.testing.allocator;
pub const TestFailingAllocator = std.testing.FailingAllocator;
