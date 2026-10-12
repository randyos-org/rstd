//! Vtable wrappers for the `axpby` kernels.
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

pub fn saxpby_k(n: BlasLong, alpha: f32, x: [*]f32, inc_x: BlasLong, beta: f32, y: [*]f32, inc_y: BlasLong) i32 {
    return blas._vtable.saxpby_k.?(n, alpha, x, inc_x, beta, y, inc_y);
}
pub fn daxpby_k(n: BlasLong, alpha: f64, x: [*]f64, inc_x: BlasLong, beta: f64, y: [*]f64, inc_y: BlasLong) i32 {
    return blas._vtable.daxpby_k.?(n, alpha, x, inc_x, beta, y, inc_y);
}
pub fn caxpby_k(n: BlasLong, alpha_r: f32, alpha_i: f32, x: [*]f32, inc_x: BlasLong, beta_r: f32, beta_i: f32, y: [*]f32, inc_y: BlasLong) i32 {
    return blas._vtable.caxpby_k.?(n, alpha_r, alpha_i, x, inc_x, beta_r, beta_i, y, inc_y);
}
pub fn zaxpby_k(n: BlasLong, alpha_r: f64, alpha_i: f64, x: [*]f64, inc_x: BlasLong, beta_r: f64, beta_i: f64, y: [*]f64, inc_y: BlasLong) i32 {
    return blas._vtable.zaxpby_k.?(n, alpha_r, alpha_i, x, inc_x, beta_r, beta_i, y, inc_y);
}

const cpu = @import("../../cpu/__root__.zig");
const sysinfo = @import("builtin");
const x86_64 = @import("x86_64/kernel.zig");
const x86 = @import("x86/kernel.zig");
const arm64 = @import("arm64/kernel.zig");
const arm = @import("arm/kernel.zig");
const loongarch64 = @import("loongarch64/kernel.zig");
const mips64 = @import("mips64/kernel.zig");
const power = @import("power/kernel.zig");
const sparc = @import("sparc/kernel.zig");

/// Points the vtable's `axpby` slots at the assembly kernels the host core
/// has, where any are converted; other cores/architectures keep whatever
/// `installGeneric`-equivalent (not yet written for this routine) would set.
pub fn install(vt: *blas.Vtable, core: cpu.CoreType) void {
    switch (sysinfo.cpu.arch) {
        .x86_64 => x86_64.install(vt, core),
        .x86 => x86.install(vt, core),
        .aarch64 => arm64.install(vt, core),
        .arm => arm.install(vt, core),
        .loongarch64 => loongarch64.install(vt, core),
        .mips, .mipsel, .mips64, .mips64el => mips64.install(vt, core),
        .powerpc, .powerpcle, .powerpc64, .powerpc64le => power.install(vt, core),
        .sparc, .sparc64 => sparc.install(vt, core),
        else => {},
    }
}
