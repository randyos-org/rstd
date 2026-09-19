//! Vtable wrappers for the `max` kernels.
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

pub fn sbmax_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas._vtable.sbmax_k.?(n, x, inc_x);
}
pub fn smax_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas._vtable.smax_k.?(n, x, inc_x);
}
pub fn dmax_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) f64 {
    return blas._vtable.dmax_k.?(n, x, inc_x);
}
pub fn qmax_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) XDouble {
    return blas._vtable.qmax_k.?(n, x, inc_x);
}
