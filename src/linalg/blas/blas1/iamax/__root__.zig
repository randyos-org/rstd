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

/// Points the vtable's `iamax` slots at the assembly kernels the host core
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
