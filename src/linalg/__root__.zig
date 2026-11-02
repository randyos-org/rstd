pub const matlab = @import("matlab/__root__.zig");
pub const vector = @import("vector.zig");
pub const matrix = @import("matrix.zig");
pub const blas = @import("blas/__root__.zig");
pub const lapack = @import("lapack/__root__.zig");

test {
    const zigstd = @import("../zigstd.zig");
    zigstd.testRefAllDecls(@This());
}
