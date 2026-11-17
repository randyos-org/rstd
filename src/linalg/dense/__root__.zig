//! `rstd.linalg`'s run-time-shape functions (design.md, "Operations over
//! views"), re-exported flat at `rstd.linalg`. Each operand is a `View` or
//! `ConstView`, a `Mat`, a pointer to a `StaticMat` or a slice (a column),
//! in any mix; a vector is a view of one column or one row. Shapes at or
//! below the inline limit take the fixed-size forms' operation order, so a
//! `Mat` and a `StaticMat` holding the same numbers give the same bits (not
//! `f64` on x87, which rounds excess precision per code path; issue #264); past
//! it, the BLAS and LAPACK. Outputs must not overlap inputs.

const reductions = @import("reductions.zig");
pub const VectorError = reductions.VectorError;
pub const axpy = reductions.axpy;
pub const dot = reductions.dot;
pub const max = reductions.max;
pub const min = reductions.min;
pub const norm = reductions.norm;
pub const scale = reductions.scale;
pub const sum = reductions.sum;

const products = @import("products.zig");
pub const MulError = products.MulError;
pub const TransposeError = products.TransposeError;
pub const ctranspose = products.ctranspose;
pub const ctransposeInto = products.ctransposeInto;
pub const mul = products.mul;
pub const mulInto = products.mulInto;
pub const transpose = products.transpose;
pub const transposeInto = products.transposeInto;

const inverse = @import("inverse.zig");
pub const DetError = inverse.DetError;
pub const InvError = inverse.InvError;
pub const SolveError = inverse.SolveError;
pub const det = inverse.det;
pub const inv = inverse.inv;
pub const invInto = inverse.invInto;
pub const solve = inverse.solve;
pub const solveInto = inverse.solveInto;

const factor = @import("factor.zig");
pub const CholError = factor.CholError;
pub const Lu = factor.Lu;
pub const Qr = factor.Qr;
pub const QrMode = factor.QrMode;
pub const chol = factor.chol;
pub const lu = factor.lu;
pub const qr = factor.qr;

const singular = @import("singular.zig");
pub const MatrixNorm = singular.MatrixNorm;
pub const MatrixNormError = singular.MatrixNormError;
pub const Svd = singular.Svd;
pub const SvdError = singular.SvdError;
pub const SvdMode = singular.SvdMode;
pub const cond = singular.cond;
pub const matrixNorm = singular.matrixNorm;
pub const pinv = singular.pinv;
pub const rank = singular.rank;
pub const svd = singular.svd;

const eigen = @import("eigen.zig");
pub const Eig = eigen.Eig;
pub const EigError = eigen.EigError;
pub const eig = eigen.eig;

test {
    _ = @import("reductions.zig");
    _ = @import("products.zig");
    _ = @import("inverse.zig");
    _ = @import("factor.zig");
    _ = @import("singular.zig");
    _ = @import("eigen.zig");
    _ = @import("util.zig");
    _ = @import("_test/products.zig");
    _ = @import("_test/inverse.zig");
    _ = @import("_test/factor.zig");
    _ = @import("_test/paths.zig");
}
