//! Vtable wrappers for the `dot` kernels.
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

pub fn sbdot_k(n: BlasLong, x: [*]BF16, inc_x: BlasLong, y: [*]BF16, inc_y: BlasLong) f32 {
    return blas.kernel(.sbdot_k)(n, x, inc_x, y, inc_y);
}
pub fn dsbdot_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong) f64 {
    return blas.kernel(.dsbdot_k)(n, x, inc_x, y, inc_y);
}
pub fn sdot_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong) f32 {
    return blas.kernel(.sdot_k)(n, x, inc_x, y, inc_y);
}
pub fn ddot_k(n: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong) f64 {
    return blas.kernel(.ddot_k)(n, x, inc_x, y, inc_y);
}
pub fn dsdot_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong) f64 {
    return blas.kernel(.dsdot_k)(n, x, inc_x, y, inc_y);
}
pub fn qdot_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong) XDouble {
    return blas.kernel(.qdot_k)(n, x, inc_x, y, inc_y);
}
pub fn cdotu_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong) Complex32 {
    return blas.kernel(.cdotu_k)(n, x, inc_x, y, inc_y).toComplex();
}
pub fn cdotc_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong) Complex32 {
    return blas.kernel(.cdotc_k)(n, x, inc_x, y, inc_y).toComplex();
}
pub fn zdotu_k(n: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong) Complex64 {
    return blas.kernel(.zdotu_k)(n, x, inc_x, y, inc_y).toComplex();
}
pub fn zdotc_k(n: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong) Complex64 {
    return blas.kernel(.zdotc_k)(n, x, inc_x, y, inc_y).toComplex();
}
pub fn xdotu_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong) ComplexXD {
    return blas.kernel(.xdotu_k)(n, x, inc_x, y, inc_y).toComplex();
}
pub fn xdotc_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong) ComplexXD {
    return blas.kernel(.xdotc_k)(n, x, inc_x, y, inc_y).toComplex();
}
