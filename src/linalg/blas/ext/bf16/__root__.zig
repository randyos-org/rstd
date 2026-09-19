//! Vtable wrappers for the `bf16` kernels.
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

pub fn sbstobf16_k(n: BlasLong, in: [*]f32, inc_in: BlasLong, out: [*]BF16, inc_out: BlasLong) void {
    blas._vtable.sbstobf16_k.?(n, in, inc_in, out, inc_out);
}
pub fn sbdtobf16_k(n: BlasLong, in: [*]f64, inc_in: BlasLong, out: [*]BF16, inc_out: BlasLong) void {
    blas._vtable.sbdtobf16_k.?(n, in, inc_in, out, inc_out);
}
pub fn sbf16tos_k(n: BlasLong, in: [*]BF16, inc_in: BlasLong, out: [*]f32, inc_out: BlasLong) void {
    blas._vtable.sbf16tos_k.?(n, in, inc_in, out, inc_out);
}
pub fn dbf16tod_k(n: BlasLong, in: [*]BF16, inc_in: BlasLong, out: [*]f64, inc_out: BlasLong) void {
    blas._vtable.dbf16tod_k.?(n, in, inc_in, out, inc_out);
}
