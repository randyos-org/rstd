//! BLAS Level 3 (matrix-matrix) routines.
//!
//! One directory per routine, holding its vtable wrappers and kernels;
//! `interface.zig` is the checks in front of `level3.zig`, the blocked driver
//! (OpenBLAS's `interface/` and `driver/level3/`), and `check.zig` the tests'
//! reference.

const zigstd = @import("../../../zigstd.zig");

pub const gemm = @import("gemm/__root__.zig");
pub const gemm3m = @import("gemm3m/__root__.zig");
pub const hemm = @import("hemm/__root__.zig");
pub const hemm3m = @import("hemm3m/__root__.zig");
pub const interface = @import("interface.zig");
pub const level3 = @import("level3.zig");
pub const symm = @import("symm/__root__.zig");
pub const symm3m = @import("symm3m/__root__.zig");
pub const syr2k = @import("syr2k/__root__.zig");
pub const syrk = @import("syrk/__root__.zig");
pub const trmm = @import("trmm/__root__.zig");
pub const trsm = @import("trsm/__root__.zig");

const Vtable = @import("../Vtable.zig");

/// Points the vtable at every level-3 routine's portable kernels, where it
/// has them.
pub fn installGeneric(vt: *Vtable) void {
    gemm.installGeneric(vt);
}

test {
    zigstd.testRefAllDecls(@This());
}
