//! Vtable wrappers for the `sum` kernels.
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

pub fn sbsum_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas._vtable.sbsum_k.?(n, x, inc_x);
}
pub fn ssum_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas._vtable.ssum_k.?(n, x, inc_x);
}
pub fn dsum_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) f64 {
    return blas._vtable.dsum_k.?(n, x, inc_x);
}
pub fn qsum_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) XDouble {
    return blas._vtable.qsum_k.?(n, x, inc_x);
}
pub fn csum_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas._vtable.csum_k.?(n, x, inc_x);
}
pub fn zsum_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) f64 {
    return blas._vtable.zsum_k.?(n, x, inc_x);
}
pub fn xsum_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) XDouble {
    return blas._vtable.xsum_k.?(n, x, inc_x);
}
