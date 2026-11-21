//! Sparse matrices (src/linalg/design.md, "Sparse matrices"): CSC first, as
//! MATLAB and SuiteSparse, with 0-based `u32` indices; CSR its row-major
//! twin over the same arrays. The kernels a sparse SQP solver needs: `A x`,
//! `Aᵀ y`, sparse `axpy` and `dot`, `scal`, extraction, Gustavson `spgemm`
//! and triangular solves (T150, #84). The iterative solvers are stubs (N2).

const csc = @import("csc.zig");
const csr = @import("csr.zig");
const triplets = @import("triplets.zig");
const vector = @import("vector.zig");
const iterative = @import("iterative.zig");

pub const Csc = csc.Csc;
pub const Csr = csr.Csr;
pub const Triplets = triplets.Triplets;
pub const Vector = vector.Vector;

pub const Error = csc.Error;
pub const SubError = csc.SubError;
pub const SelectError = csc.SelectError;
pub const MulError = csc.MulError;
pub const SpgemmError = csc.SpgemmError;
pub const NumericError = csc.NumericError;
pub const TriSolveError = csc.TriSolveError;
pub const VectorError = vector.VectorError;

pub const Operator = iterative.Operator;
pub const Options = iterative.Options;
pub const Stats = iterative.Stats;
pub const SolveError = iterative.SolveError;
pub const cg = iterative.cg;
pub const gmres = iterative.gmres;

test {
    const zigstd = @import("../../zigstd.zig");
    zigstd.testRefAllDecls(@This());
    _ = @import("mul.zig");
    _ = @import("spgemm.zig");
    _ = @import("triangular.zig");
    _ = @import("_test/build.zig");
    _ = @import("_test/mul.zig");
    _ = @import("_test/spgemm.zig");
    _ = @import("_test/triangular.zig");
}
