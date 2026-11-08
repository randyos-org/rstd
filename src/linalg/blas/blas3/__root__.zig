//! BLAS Level 3 (matrix-matrix) routines.
//!
//! One directory per routine, holding its kernels and tests; `interface.zig`
//! checks the arguments in front of the blocked drivers (`driver.zig`,
//! `driver3m.zig`, `triangular.zig`; OpenBLAS's `interface/` and
//! `driver/level3/`), and `_test/blas3/check.zig` is the tests' reference.

pub const gemm = @import("gemm/__root__.zig");
pub const gemm3m = @import("gemm3m/__root__.zig");
pub const hemm3m = @import("hemm3m/__root__.zig");
pub const symm = @import("symm/__root__.zig");
pub const symm3m = @import("symm3m/__root__.zig");
pub const syr2k = @import("syr2k/__root__.zig");
pub const syrk = @import("syrk/__root__.zig");
pub const trmm = @import("trmm/__root__.zig");
pub const trsm = @import("trsm/__root__.zig");
/// `sbgemm`, `bgemm`, `shgemm` and `sbgemmt`.
pub const half = @import("half.zig");
/// `?gemm_batch` and `?gemm_batch_strided`, the half gemms' too.
pub const batch = @import("batch.zig");
pub const interface = @import("interface.zig");
pub const driver = @import("driver.zig");
/// The 3M drivers (`gemm3m`, `symm3m`, `hemm3m`).
pub const driver3m = @import("driver3m.zig");
/// The level-3 kernel registry (`gemm/kernel.zig`'s tile and the sets over it).
pub const kernels = @import("kernels.zig");
/// The `trmm` and `trsm` drivers.
pub const triangular = @import("triangular.zig");
/// The threaded level-3 driver (`threading.md`).
pub const threaded = @import("threaded.zig");

const Vtable = @import("../Vtable.zig");
const zigstd = @import("../../../zigstd.zig");

/// Points the vtable at every level-3 routine's portable kernels, where it
/// has them.
pub fn installGeneric(vt: *Vtable) void {
    gemm.installGeneric(vt);
    half.installGeneric(vt);
}

test {
    zigstd.testRefAllDecls(@This());
}
