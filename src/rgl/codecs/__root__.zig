pub const ico = @import("ico.zig");
pub const png = @import("png.zig");

const zigstd = @import("../../zigstd.zig");

test {
    zigstd.testRefAllDecls(@This());
}
