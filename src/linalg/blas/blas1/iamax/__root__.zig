//! Vtable wrappers for the `iamax` kernels.
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

pub fn isbamax_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) BlasLong {
    return blas._vtable.isbamax_k.?(n, x, inc_x);
}
pub fn isamax_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) BlasLong {
    return blas._vtable.isamax_k.?(n, x, inc_x);
}
pub fn idamax_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) BlasLong {
    return blas._vtable.idamax_k.?(n, x, inc_x);
}
pub fn iqamax_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) BlasLong {
    return blas._vtable.iqamax_k.?(n, x, inc_x);
}
pub fn icamax_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) BlasLong {
    return blas._vtable.icamax_k.?(n, x, inc_x);
}
pub fn izamax_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) BlasLong {
    return blas._vtable.izamax_k.?(n, x, inc_x);
}
pub fn ixamax_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) BlasLong {
    return blas._vtable.ixamax_k.?(n, x, inc_x);
}
