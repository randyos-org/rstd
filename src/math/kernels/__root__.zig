//! `rstd.math.kernels`: the kernel layer, every optimized routine with every
//! knob explicit (`../design.md`, issue #606).

pub const blas = @import("blas/__root__.zig");
pub const lapack = @import("lapack/__root__.zig");
pub const sparse = @import("sparse/__root__.zig");

test {
    _ = blas;
    _ = lapack;
    _ = sparse;
}
