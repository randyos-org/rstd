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

const subspace = @import("subspace.zig");
pub const orth = subspace.orth;
pub const nullSpace = subspace.nullSpace;

const condition = @import("condition.zig");
pub const RcondError = condition.RcondError;
pub const rcond = condition.rcond;

const vector_norm = @import("vector_norm.zig");
pub const VectorNorm = vector_norm.VectorNorm;
pub const VectorNormError = vector_norm.VectorNormError;
pub const vectorNorm = vector_norm.vectorNorm;

const elementwise = @import("elementwise.zig");
pub const ElementwiseError = elementwise.ElementwiseError;
pub const ElementwiseIntoError = elementwise.ElementwiseIntoError;
pub const add = elementwise.add;
pub const sub = elementwise.sub;
pub const times = elementwise.times;
pub const rdivide = elementwise.rdivide;
pub const addInto = elementwise.addInto;
pub const subInto = elementwise.subInto;
pub const timesInto = elementwise.timesInto;
pub const rdivideInto = elementwise.rdivideInto;
pub const prod = elementwise.prod;
pub const map = elementwise.map;
pub const mapInto = elementwise.mapInto;
pub const realPart = elementwise.realPart;
pub const imagPart = elementwise.imagPart;
pub const conj = elementwise.conj;

const schur_form = @import("schur_form.zig");
pub const HessError = schur_form.HessError;
pub const Hess = schur_form.Hess;
pub const hess = schur_form.hess;
pub const SchurError = schur_form.SchurError;
pub const Schur = schur_form.Schur;
pub const schur = schur_form.schur;
pub const schurComplex = schur_form.schurComplex;

const expm_ = @import("expm.zig");
pub const ExpmError = expm_.ExpmError;
pub const expm = expm_.expm;

test {
    _ = @import("expm.zig");
    _ = @import("schur_form.zig");
    _ = @import("elementwise.zig");
    _ = @import("vector_norm.zig");
    _ = @import("condition.zig");
    _ = @import("subspace.zig");
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
