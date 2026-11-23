//! What the Accelerate tests share.

const std = @import("std");

pub fn fill(comptime F: type, rng: std.Random, out: []F) void {
    for (out) |*v| v.* = rng.float(F) * 2 - 1;
}
