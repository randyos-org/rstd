//! `rstd.math.kernels.sparse` (`../../design.md`, #610): compressed sparse
//! columns first, as MATLAB and SuiteSparse, with 0-based `u32` indices;
//! compressed rows its row-major twin over the same arrays. The kernels a
//! sparse SQP solver needs: `A x`, `Aᵀ y`, sparse `axpy` and `dot`, `scal`,
//! extraction, Gustavson `spgemm` and triangular solves (T150, #84). The
//! iterative solvers: conjugate gradients and restarted GMRES.

const csc = @import("csc.zig");
const csr = @import("csr.zig");
const triplets = @import("triplets.zig");
const vector = @import("vector.zig");
const iterative = @import("iterative.zig");
const cg_solver = @import("cg.zig");
const gmres_solver = @import("gmres.zig");
const zigstd = @import("../../../zigstd.zig");

pub const CompressedColumns = csc.CompressedColumns;
pub const CompressedRows = csr.CompressedRows;
pub const Triplets = triplets.Triplets;
pub const SparseVector = vector.SparseVector;

pub const Error = csc.Error;
pub const FromDenseError = csc.FromDenseError;
pub const ToDenseError = csc.ToDenseError;
pub const SubError = csc.SubError;
pub const SelectError = csc.SelectError;
pub const MulError = csc.MulError;
pub const SpgemmError = csc.SpgemmError;
pub const NumericError = csc.NumericError;
pub const TriSolveError = csc.TriSolveError;
pub const VectorError = vector.VectorError;

pub const LinearOperator = iterative.LinearOperator;
pub const SolveReport = iterative.SolveReport;
pub const SolveError = iterative.SolveError;
pub const cg = cg_solver.cg;
pub const gmres = gmres_solver.gmres;

test {
    zigstd.testRefAllDecls(@This());
    _ = @import("mul.zig");
    _ = @import("spgemm.zig");
    _ = @import("triangular.zig");
    _ = @import("_test/build.zig");
    _ = @import("_test/iterative.zig");
    _ = @import("_test/mul.zig");
    _ = @import("_test/spgemm.zig");
    _ = @import("_test/triangular.zig");
}
