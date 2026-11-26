//! Level 3's reference-shaped core in one place (pointers, sizes, leading
//! dimensions, the reference's named errors), for the LAPACK routines, the
//! ported tests and the other internal callers; `api.zig` is the public layer.

const iface = @import("interface.zig");
const half = @import("half.zig");
const batch = @import("batch.zig");

pub const Error = iface.Error;
pub const ArgError = iface.ArgError;
pub const HalfError = half.HalfError;
pub const BatchError = batch.Error;
pub const HalfPrecision = half.HalfPrecision;
pub const GemmGroup = batch.GemmGroup;

pub const gemm = iface.gemm;
pub const gemm3m = iface.gemm3m;
pub const gemmt = iface.gemmt;
pub const symm = iface.symm;
pub const hemm = iface.hemm;
pub const symm3m = iface.symm3m;
pub const hemm3m = iface.hemm3m;
pub const syrk = iface.syrk;
pub const herk = iface.herk;
pub const syr2k = iface.syr2k;
pub const her2k = iface.her2k;
pub const trmm = iface.trmm;
pub const trsm = iface.trsm;
pub const gemmHalf = half.gemmHalf;
pub const gemmtHalf = half.gemmtHalf;
pub const gemmBatch = batch.gemmBatch;
pub const gemmBatchStrided = batch.gemmBatchStrided;
pub const gemmHalfBatch = batch.gemmHalfBatch;
pub const gemmHalfBatchStrided = batch.gemmHalfBatchStrided;
