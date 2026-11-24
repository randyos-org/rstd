const threading = @import("threading.zig");
const test_io = @import("_test/io.zig");
const select = @import("select.zig");

pub const types = @import("types.zig");
pub const cpu = @import("cpu/__root__.zig");
pub const Vtable = @import("Vtable.zig");
/// Strided vector views and the shared error pieces (`strided.zig`).
pub const strided = @import("strided.zig");
pub const Strided = strided.Strided;
pub const ConstStrided = strided.ConstStrided;
pub const Config = @import("Config.zig");
/// The `@Vector` kernels some cores install; `tools/blasbench` races them.
pub const vector = @import("vector.zig");
const blas1 = @import("blas1/__root__.zig");
const blas2 = @import("blas2/__root__.zig");
const blas3 = @import("blas3/__root__.zig");
const ext = @import("ext/__root__.zig");
const l1 = blas1.api;
const l2 = blas2.interface;
const l2api = blas2.api;
const tri = blas2.triangular;
const band = blas2.banded;
const rank = blas2.rank;
const l3 = blas3.interface;
const l3api = blas3.api;

// ---- the BLAS routines ---------------------------------------------------------
//
// Level 1 on strided views (`blas1/api.zig`); levels 2 and 3 in Fortran
// argument order, column-major, sizes and strides `isize`. The element type
// explicit. Each takes `io` after its allocator, or first without one
// (`threading.md`); bad arguments return named errors in `xerbla`'s order.

pub const Trans = types.Trans;
pub const Uplo = types.Uplo;
pub const Diag = types.Diag;
pub const Side = types.Side;
pub const Extreme = l1.Extreme;
pub const HalfPrecision = l2.HalfPrecision;

/// Level 1, two vectors: counts that disagree, or `error.Canceled`.
pub const Level1Error = l1.PairError;
/// Level 2: counts that disagree, a slice too short for its matrix or a
/// leading dimension too small; `error.OutOfMemory` or `error.Canceled`.
pub const Level2Error = l2api.Error;
/// `gemvHalf`, which allocates nothing: `Level2Error` less `OutOfMemory`.
pub const GemvHalfError = l2api.HalfError;
/// Level 3: a slice too short for its matrix or a leading dimension too
/// small; `error.OutOfMemory` or `error.Canceled`.
pub const Level3Error = l3api.Error;
/// `gemmtHalf`, which allocates nothing: `Level3Error` less `OutOfMemory`.
pub const GemmtHalfError = l3api.HalfError;
/// The batched gemms: `Level3Error` and `error.LengthMismatch`.
pub const BatchError = l3api.BatchError;
/// One group of `gemmBatch` or `gemmHalfBatch`, by input and output type.
pub const GemmGroup = l3api.GemmGroup;
/// `laswp`: a row or pivot outside the matrix, a slice too short or a leading
/// dimension too small.
pub const LaswpError = ext.api.LaswpError;
/// `geadd`: an unsupported letter, a slice too short or a leading dimension
/// too small.
pub const GeaddError = ext.api.GeaddError;
/// `omatcopy`: `GeaddError`, `error.Overlap` where `A` and `B` share memory,
/// and `error.Canceled` from `io`, which large real copies thread on.
pub const OmatcopyError = ext.api.OmatcopyError;
/// `imatcopy`: `GeaddError`, and `error.OutOfMemory`.
pub const ImatcopyError = ext.api.ImatcopyError;

pub const copy = l1.copy;
pub const swap = l1.swap;
pub const scal = l1.scal;
pub const scalReal = l1.scalReal;
pub const dot = l1.dot;
pub const dotc = l1.dotc;
pub const dsdot = l1.dsdot;
pub const sdsdot = l1.sdsdot;
pub const nrm2 = l1.nrm2;
pub const asum = l1.asum;
pub const sum = l1.sum;
pub const sumWith = l1.sumWith;
pub const asumWith = l1.asumWith;
pub const dotWith = l1.dotWith;
pub const dotcWith = l1.dotcWith;
pub const Summation = l1.Summation;
pub const index = l1.index;
pub const extreme = l1.extreme;
pub const rot = l1.rot;
pub const rotm = l1.rotm;
pub const axpy = l1.axpy;
pub const axpyc = l1.axpyc;
pub const axpby = l1.axpby;
pub const sbdot = l1.sbdot;
pub const bscal = l1.bscal;
pub const sbtobf16 = l1.sbtobf16;
pub const bf16to = l1.bf16to;
pub const rotg = l1.rotg;
pub const zrotg = l1.zrotg;
pub const rotmg = l1.rotmg;
pub const cabs1 = l1.cabs1;

pub const gemv = l2api.gemv;
pub const gemvHalf = l2api.gemvHalf;
pub const ger = l2api.ger;
pub const gerc = l2api.gerc;
pub const symv = l2api.symv;
pub const hemv = l2api.hemv;
pub const trmv = l2api.trmv;
pub const trsv = l2api.trsv;
pub const tbmv = l2api.tbmv;
pub const tbsv = l2api.tbsv;
pub const tpmv = l2api.tpmv;
pub const tpsv = l2api.tpsv;
pub const gbmv = l2api.gbmv;
pub const sbmv = l2api.sbmv;
pub const hbmv = l2api.hbmv;
pub const spmv = l2api.spmv;
pub const hpmv = l2api.hpmv;
pub const syr = l2api.syr;
pub const spr = l2api.spr;
pub const her = l2api.her;
pub const hpr = l2api.hpr;
pub const syr2 = l2api.syr2;
pub const spr2 = l2api.spr2;
pub const her2 = l2api.her2;
pub const hpr2 = l2api.hpr2;

pub const gemm = l3api.gemm;
pub const symm = l3api.symm;
pub const hemm = l3api.hemm;
pub const syrk = l3api.syrk;
pub const herk = l3api.herk;
pub const syr2k = l3api.syr2k;
pub const her2k = l3api.her2k;
pub const trmm = l3api.trmm;
pub const trsm = l3api.trsm;
pub const gemmt = l3api.gemmt;
pub const gemm3m = l3api.gemm3m;
pub const symm3m = l3api.symm3m;
pub const hemm3m = l3api.hemm3m;

pub const geadd = ext.api.geadd;
pub const omatcopy = ext.api.omatcopy;
pub const imatcopy = ext.api.imatcopy;
pub const laswp = ext.api.laswp;

pub const gemmHalf = l3api.gemmHalf;
pub const gemmtHalf = l3api.gemmtHalf;
pub const gemmBatch = l3api.gemmBatch;
pub const gemmBatchStrided = l3api.gemmBatchStrided;
pub const gemmHalfBatch = l3api.gemmHalfBatch;
pub const gemmHalfBatchStrided = l3api.gemmHalfBatchStrided;

// ---- selection and the vtable (`select.zig`) --------------------------------------

pub var config: Config = .{};

pub const Slot = select.Slot;
pub const vtable = select.vtable;
pub const table = select.table;
pub const tableOf = select.tableOf;
pub const kernel = select.kernel;
pub const init = select.init;
pub const shutdown = select.shutdown;
pub const useAccelerate = select.useAccelerate;
pub const coreType = select.coreType;
pub const coreName = select.coreName;
pub const TestSaved = select.TestSaved;
pub const testUse = select.testUse;
pub const testTable = select.testTable;
pub const testRestore = select.testRestore;

const zigstd = @import("../../../zigstd.zig");

test {
    _ = @import("abi.zig");
    _ = @import("FpState.zig");
    _ = @import("apple_accel.zig");
    // hosted only: the rigs' `math.kernels.blas._test.` filter would take it
    if (comptime test_io.hosted) _ = @import("_test/apple_accel/__root__.zig");
    _ = @import("threading.zig");
    _ = @import("threading/_test/plan.zig");
    _ = @import("threading/_test/region.zig");
    _ = @import("select.zig");
    _ = @import("_test/select.zig");
    _ = @import("strided.zig");
    _ = @import("_test/__root__.zig");
    zigstd.testRefAllDecls(@This());
}
