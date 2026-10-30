//! Vtable wrappers for the `symm3m` kernels.
//!
//! Argument names are copied from the OpenBLAS kernel sources (see the
//! `CNAME(...)` definitions under `kernel/<arch>/*.c`). Several are named
//! `dummy*` in the reference: kernel-ABI padding slots this shape does not use.

const blas = @import("../../__root__.zig");
const types = @import("../../types.zig");

const BlasLong = types.BlasLong;
const HF16 = types.HF16;
const BF16 = types.BF16;
const XDouble = types.XDouble;
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;
const ComplexXD = types.ComplexXD;

pub fn csymm3m_iucopyb(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.csymm3m_iucopyb)(m, n, a, lda, posX, posY, b);
}
pub fn csymm3m_ilcopyb(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.csymm3m_ilcopyb)(m, n, a, lda, posX, posY, b);
}
pub fn csymm3m_iucopyr(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.csymm3m_iucopyr)(m, n, a, lda, posX, posY, b);
}
pub fn csymm3m_ilcopyr(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.csymm3m_ilcopyr)(m, n, a, lda, posX, posY, b);
}
pub fn csymm3m_iucopyi(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.csymm3m_iucopyi)(m, n, a, lda, posX, posY, b);
}
pub fn csymm3m_ilcopyi(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f32) i32 {
    return blas.kernel(.csymm3m_ilcopyi)(m, n, a, lda, posX, posY, b);
}
pub fn csymm3m_oucopyb(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: f32, alpha_i: f32, b: [*]f32) i32 {
    return blas.kernel(.csymm3m_oucopyb)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn csymm3m_olcopyb(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: f32, alpha_i: f32, b: [*]f32) i32 {
    return blas.kernel(.csymm3m_olcopyb)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn csymm3m_oucopyr(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: f32, alpha_i: f32, b: [*]f32) i32 {
    return blas.kernel(.csymm3m_oucopyr)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn csymm3m_olcopyr(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: f32, alpha_i: f32, b: [*]f32) i32 {
    return blas.kernel(.csymm3m_olcopyr)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn csymm3m_oucopyi(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: f32, alpha_i: f32, b: [*]f32) i32 {
    return blas.kernel(.csymm3m_oucopyi)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn csymm3m_olcopyi(m: BlasLong, n: BlasLong, a: [*]f32, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: f32, alpha_i: f32, b: [*]f32) i32 {
    return blas.kernel(.csymm3m_olcopyi)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn zsymm3m_iucopyb(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zsymm3m_iucopyb)(m, n, a, lda, posX, posY, b);
}
pub fn zsymm3m_ilcopyb(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zsymm3m_ilcopyb)(m, n, a, lda, posX, posY, b);
}
pub fn zsymm3m_iucopyr(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zsymm3m_iucopyr)(m, n, a, lda, posX, posY, b);
}
pub fn zsymm3m_ilcopyr(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zsymm3m_ilcopyr)(m, n, a, lda, posX, posY, b);
}
pub fn zsymm3m_iucopyi(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zsymm3m_iucopyi)(m, n, a, lda, posX, posY, b);
}
pub fn zsymm3m_ilcopyi(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]f64) i32 {
    return blas.kernel(.zsymm3m_ilcopyi)(m, n, a, lda, posX, posY, b);
}
pub fn zsymm3m_oucopyb(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: f64, alpha_i: f64, b: [*]f64) i32 {
    return blas.kernel(.zsymm3m_oucopyb)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn zsymm3m_olcopyb(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: f64, alpha_i: f64, b: [*]f64) i32 {
    return blas.kernel(.zsymm3m_olcopyb)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn zsymm3m_oucopyr(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: f64, alpha_i: f64, b: [*]f64) i32 {
    return blas.kernel(.zsymm3m_oucopyr)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn zsymm3m_olcopyr(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: f64, alpha_i: f64, b: [*]f64) i32 {
    return blas.kernel(.zsymm3m_olcopyr)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn zsymm3m_oucopyi(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: f64, alpha_i: f64, b: [*]f64) i32 {
    return blas.kernel(.zsymm3m_oucopyi)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn zsymm3m_olcopyi(m: BlasLong, n: BlasLong, a: [*]f64, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: f64, alpha_i: f64, b: [*]f64) i32 {
    return blas.kernel(.zsymm3m_olcopyi)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn xsymm3m_iucopyb(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xsymm3m_iucopyb)(m, n, a, lda, posX, posY, b);
}
pub fn xsymm3m_ilcopyb(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xsymm3m_ilcopyb)(m, n, a, lda, posX, posY, b);
}
pub fn xsymm3m_iucopyr(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xsymm3m_iucopyr)(m, n, a, lda, posX, posY, b);
}
pub fn xsymm3m_ilcopyr(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xsymm3m_ilcopyr)(m, n, a, lda, posX, posY, b);
}
pub fn xsymm3m_iucopyi(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xsymm3m_iucopyi)(m, n, a, lda, posX, posY, b);
}
pub fn xsymm3m_ilcopyi(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]XDouble) i32 {
    return blas.kernel(.xsymm3m_ilcopyi)(m, n, a, lda, posX, posY, b);
}
pub fn xsymm3m_oucopyb(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: XDouble, alpha_i: XDouble, b: [*]XDouble) i32 {
    return blas.kernel(.xsymm3m_oucopyb)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn xsymm3m_olcopyb(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: XDouble, alpha_i: XDouble, b: [*]XDouble) i32 {
    return blas.kernel(.xsymm3m_olcopyb)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn xsymm3m_oucopyr(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: XDouble, alpha_i: XDouble, b: [*]XDouble) i32 {
    return blas.kernel(.xsymm3m_oucopyr)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn xsymm3m_olcopyr(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: XDouble, alpha_i: XDouble, b: [*]XDouble) i32 {
    return blas.kernel(.xsymm3m_olcopyr)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn xsymm3m_oucopyi(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: XDouble, alpha_i: XDouble, b: [*]XDouble) i32 {
    return blas.kernel(.xsymm3m_oucopyi)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}
pub fn xsymm3m_olcopyi(m: BlasLong, n: BlasLong, a: [*]XDouble, lda: BlasLong, posX: BlasLong, posY: BlasLong, alpha_r: XDouble, alpha_i: XDouble, b: [*]XDouble) i32 {
    return blas.kernel(.xsymm3m_olcopyi)(m, n, a, lda, posX, posY, alpha_r, alpha_i, b);
}

// ---- tests -----------------------------------------------------------------

const check = @import("../../_test/blas3/check.zig");
const blas1check = @import("../../_test/blas1/check.zig");
const Vtable = @import("../../Vtable.zig");
const portable = @import("../../portable.zig");

fn cases(vt: *const Vtable) !void {
    try check.routineAll(vt, .symm3m);
}

test "portable symm3m kernels match the wide reference through the driver" {
    var vt: Vtable = .{};
    portable.install(&vt);
    try cases(&vt);
}

test "every linked core's symm3m kernels match the wide reference, where this machine runs them" {
    try blas1check.eachCore(cases);
}
