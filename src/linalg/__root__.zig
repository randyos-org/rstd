pub const matlab = @import("matlab/__root__.zig");
pub const vector = @import("vector.zig");
pub const blas = @import("blas/__root__.zig");

test {
    const zigstd = @import("../zigstd.zig");
    zigstd.testRefAllDecls(@This());
}
