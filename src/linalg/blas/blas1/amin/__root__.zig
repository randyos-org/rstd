//! Vtable wrappers for the `amin` kernels.
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

pub fn sbamin_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas.kernel(.sbamin_k)(n, x, inc_x);
}
pub fn samin_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas.kernel(.samin_k)(n, x, inc_x);
}
pub fn damin_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) f64 {
    return blas.kernel(.damin_k)(n, x, inc_x);
}
pub fn qamin_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) XDouble {
    return blas.kernel(.qamin_k)(n, x, inc_x);
}
pub fn camin_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) f32 {
    return blas.kernel(.camin_k)(n, x, inc_x);
}
pub fn zamin_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) f64 {
    return blas.kernel(.zamin_k)(n, x, inc_x);
}
pub fn xamin_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) XDouble {
    return blas.kernel(.xamin_k)(n, x, inc_x);
}
