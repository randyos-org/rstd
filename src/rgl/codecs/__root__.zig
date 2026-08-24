const std = @import("std");

pub const ico = @import("ico.zig");
pub const png = @import("png.zig");

test {
    std.testing.refAllDecls(@This());
}
