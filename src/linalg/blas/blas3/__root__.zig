//! BLAS Level 3 (matrix-matrix) routines.
//!
//! One file per routine, holding the vtable wrappers for that routine.

const zigstd = @import("../../../zigstd.zig");

pub const gemm = @import("gemm/__root__.zig");
pub const gemm3m = @import("gemm3m/__root__.zig");
pub const hemm = @import("hemm/__root__.zig");
pub const hemm3m = @import("hemm3m/__root__.zig");
pub const symm = @import("symm/__root__.zig");
pub const symm3m = @import("symm3m/__root__.zig");
pub const syr2k = @import("syr2k/__root__.zig");
pub const syrk = @import("syrk/__root__.zig");
pub const trmm = @import("trmm/__root__.zig");
pub const trsm = @import("trsm/__root__.zig");

test {
    zigstd.testRefAllDecls(@This());
}
