const zigstd = @import("../../zigstd.zig");

pub const ico = @import("ico.zig");
pub const png = @import("png.zig");

test {
    zigstd.testRefAllDecls(@This());
}
