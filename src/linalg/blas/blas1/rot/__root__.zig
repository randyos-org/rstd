//! Vtable wrappers for the `rot` kernels.
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

pub fn sbrot_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, c: f32, s: f32) i32 {
    return blas.kernel(.sbrot_k)(n, x, inc_x, y, inc_y, c, s);
}
pub fn srot_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, c: f32, s: f32) i32 {
    return blas.kernel(.srot_k)(n, x, inc_x, y, inc_y, c, s);
}
pub fn drot_k(n: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, c: f64, s: f64) i32 {
    return blas.kernel(.drot_k)(n, x, inc_x, y, inc_y, c, s);
}
pub fn qrot_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, c: XDouble, s: XDouble) i32 {
    return blas.kernel(.qrot_k)(n, x, inc_x, y, inc_y, c, s);
}
pub fn csrot_k(n: BlasLong, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, c: f32, s: f32) i32 {
    return blas.kernel(.csrot_k)(n, x, inc_x, y, inc_y, c, s);
}
pub fn zdrot_k(n: BlasLong, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, c: f64, s: f64) i32 {
    return blas.kernel(.zdrot_k)(n, x, inc_x, y, inc_y, c, s);
}
pub fn xqrot_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, c: XDouble, s: XDouble) i32 {
    return blas.kernel(.xqrot_k)(n, x, inc_x, y, inc_y, c, s);
}
