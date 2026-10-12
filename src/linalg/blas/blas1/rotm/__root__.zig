//! Vtable wrappers for the `rotm` kernels.
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

pub fn sbrotm_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, param: [*]f32) i32 {
    return blas.kernel(.sbrotm_k)(n, x, inc_x, y, inc_y, param);
}
pub fn srotm_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, param: [*]f32) i32 {
    return blas.kernel(.srotm_k)(n, x, inc_x, y, inc_y, param);
}
pub fn drotm_k(n: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, param: [*]f64) i32 {
    return blas.kernel(.drotm_k)(n, x, inc_x, y, inc_y, param);
}
pub fn qrotm_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, param: [*]XDouble) i32 {
    return blas.kernel(.qrotm_k)(n, x, inc_x, y, inc_y, param);
}
