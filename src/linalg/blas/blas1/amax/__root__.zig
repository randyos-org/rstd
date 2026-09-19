//! Vtable wrappers for the `amax` kernels.
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

pub fn sbamax_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas._vtable.sbamax_k.?(n, x, inc_x);
}
pub fn samax_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas._vtable.samax_k.?(n, x, inc_x);
}
pub fn damax_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) f64 {
    return blas._vtable.damax_k.?(n, x, inc_x);
}
pub fn qamax_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) XDouble {
    return blas._vtable.qamax_k.?(n, x, inc_x);
}
pub fn camax_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas._vtable.camax_k.?(n, x, inc_x);
}
pub fn zamax_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) f64 {
    return blas._vtable.zamax_k.?(n, x, inc_x);
}
pub fn xamax_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) XDouble {
    return blas._vtable.xamax_k.?(n, x, inc_x);
}
