//! Vtable wrappers for the `copy` kernels.
//!
//! Argument names are copied from the OpenBLAS kernel sources (see the
//! `CNAME(...)` definitions under `kernel/<arch>/*.c`). Several are named
//! `dummy*` in the reference: kernel-ABI padding slots this shape does not use.

const blas = @import("../../__root__.zig");
const types = @import("../../types.zig");

const BlasLong = types.BlasLong;
const BlasInt = types.BlasInt;
const HF16 = types.HF16;
const BF16 = types.BF16;
const XDouble = types.XDouble;
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;
const ComplexXD = types.ComplexXD;

pub fn sbcopy_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong) i32 {
    return blas.kernel(.sbcopy_k)(n, x, inc_x, y, inc_y);
}
pub fn scopy_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong) i32 {
    return blas.kernel(.scopy_k)(n, x, inc_x, y, inc_y);
}
pub fn dcopy_k(n: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong) i32 {
    return blas.kernel(.dcopy_k)(n, x, inc_x, y, inc_y);
}
pub fn qcopy_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong) i32 {
    return blas.kernel(.qcopy_k)(n, x, inc_x, y, inc_y);
}
pub fn ccopy_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong) i32 {
    return blas.kernel(.ccopy_k)(n, x, inc_x, y, inc_y);
}
pub fn zcopy_k(n: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong) i32 {
    return blas.kernel(.zcopy_k)(n, x, inc_x, y, inc_y);
}
pub fn xcopy_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong) i32 {
    return blas.kernel(.xcopy_k)(n, x, inc_x, y, inc_y);
}
