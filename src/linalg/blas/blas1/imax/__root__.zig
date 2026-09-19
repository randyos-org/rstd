//! Vtable wrappers for the `imax` kernels.
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

pub fn isbmax_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) BlasLong {
    return blas._vtable.isbmax_k.?(n, x, inc_x);
}
pub fn ismax_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) BlasLong {
    return blas._vtable.ismax_k.?(n, x, inc_x);
}
pub fn idmax_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) BlasLong {
    return blas._vtable.idmax_k.?(n, x, inc_x);
}
pub fn iqmax_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) BlasLong {
    return blas._vtable.iqmax_k.?(n, x, inc_x);
}
