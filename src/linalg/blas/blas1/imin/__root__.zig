//! Vtable wrappers for the `imin` kernels.
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

pub fn isbmin_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) BlasLong {
    return blas.kernel(.isbmin_k)(n, x, inc_x);
}
pub fn ismin_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) BlasLong {
    return blas.kernel(.ismin_k)(n, x, inc_x);
}
pub fn idmin_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) BlasLong {
    return blas.kernel(.idmin_k)(n, x, inc_x);
}
pub fn iqmin_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) BlasLong {
    return blas.kernel(.iqmin_k)(n, x, inc_x);
}
