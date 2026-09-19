//! Vtable wrappers for the `iamin` kernels.
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

pub fn isbamin_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) BlasLong {
    return blas._vtable.isbamin_k.?(n, x, inc_x);
}
pub fn isamin_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) BlasLong {
    return blas._vtable.isamin_k.?(n, x, inc_x);
}
pub fn idamin_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) BlasLong {
    return blas._vtable.idamin_k.?(n, x, inc_x);
}
pub fn iqamin_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) BlasLong {
    return blas._vtable.iqamin_k.?(n, x, inc_x);
}
pub fn icamin_k(n: BlasLong, x: [*]f32, inc_x: BlasLong) BlasLong {
    return blas._vtable.icamin_k.?(n, x, inc_x);
}
pub fn izamin_k(n: BlasLong, x: [*]f64, inc_x: BlasLong) BlasLong {
    return blas._vtable.izamin_k.?(n, x, inc_x);
}
pub fn ixamin_k(n: BlasLong, x: [*]XDouble, inc_x: BlasLong) BlasLong {
    return blas._vtable.ixamin_k.?(n, x, inc_x);
}
