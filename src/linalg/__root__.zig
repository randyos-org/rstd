pub const matlab = @import("matlab/__root__.zig");
pub const vector = @import("vector.zig");

test {
    const std = @import("std");
    std.testing.refAllDecls(@This());
}
