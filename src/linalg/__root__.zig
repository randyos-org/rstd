pub const matlab = @import("matlab/__root__.zig");
pub const vector = @import("vector.zig");
pub const matrix = @import("matrix.zig");
pub const blas = @import("blas/__root__.zig");
pub const lapack = @import("lapack/__root__.zig");

// the public layer (design.md, "Planned: the public `rstd.linalg` layer"):
// its types, and `dense.zig`'s run-time-shape functions, flat
const dense = @import("dense.zig");
pub const Vec = vector.Vec;
pub const StaticMat = matrix.StaticMat;
pub const Mat = matrix.Mat;
pub const View = matrix.View;
pub const ConstView = matrix.ConstView;
pub const VectorError = dense.VectorError;
pub const MulError = dense.MulError;
pub const TransposeError = dense.TransposeError;
pub const DetError = dense.DetError;
pub const InvError = dense.InvError;
pub const SolveError = dense.SolveError;
pub const CholError = dense.CholError;
pub const MatrixNormError = dense.MatrixNormError;
pub const Lu = dense.Lu;
pub const Qr = dense.Qr;
pub const QrMode = dense.QrMode;
pub const MatrixNorm = dense.MatrixNorm;
pub const dot = dense.dot;
pub const norm = dense.norm;
pub const sum = dense.sum;
pub const max = dense.max;
pub const min = dense.min;
pub const axpy = dense.axpy;
pub const scale = dense.scale;
pub const mul = dense.mul;
pub const mulInto = dense.mulInto;
pub const transpose = dense.transpose;
pub const ctranspose = dense.ctranspose;
pub const transposeInto = dense.transposeInto;
pub const ctransposeInto = dense.ctransposeInto;
pub const det = dense.det;
pub const inv = dense.inv;
pub const invInto = dense.invInto;
pub const solve = dense.solve;
pub const solveInto = dense.solveInto;
pub const lu = dense.lu;
pub const chol = dense.chol;
pub const qr = dense.qr;
pub const matrixNorm = dense.matrixNorm;

test {
    const zigstd = @import("../zigstd.zig");
    zigstd.testRefAllDecls(@This());
}
