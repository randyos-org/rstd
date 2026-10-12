//! Vtable wrappers for the `swap` kernels.
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

pub fn sbswap_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, dummy3: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, dummy: [*]f32, dummy2: BlasLong) i32 {
    return blas._vtable.sbswap_k.?(n, dummy0, dummy1, dummy3, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn sswap_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, dummy3: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, dummy: [*]f32, dummy2: BlasLong) i32 {
    return blas._vtable.sswap_k.?(n, dummy0, dummy1, dummy3, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn dswap_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, dummy3: f64, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, dummy: [*]f64, dummy2: BlasLong) i32 {
    return blas._vtable.dswap_k.?(n, dummy0, dummy1, dummy3, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn qswap_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, dummy3: XDouble, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, dummy: [*]XDouble, dummy2: BlasLong) i32 {
    return blas._vtable.qswap_k.?(n, dummy0, dummy1, dummy3, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn cswap_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, dummy3: f32, dummy4: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, dummy: [*]f32, dummy2: BlasLong) i32 {
    return blas._vtable.cswap_k.?(n, dummy0, dummy1, dummy3, dummy4, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn zswap_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, dummy3: f64, dummy4: f64, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, dummy: [*]f64, dummy2: BlasLong) i32 {
    return blas._vtable.zswap_k.?(n, dummy0, dummy1, dummy3, dummy4, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn xswap_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, dummy3: XDouble, dummy4: XDouble, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, dummy: [*]XDouble, dummy2: BlasLong) i32 {
    return blas._vtable.xswap_k.?(n, dummy0, dummy1, dummy3, dummy4, x, inc_x, y, inc_y, dummy, dummy2);
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

/// Points the vtable's `swap` slots at the assembly kernels the host core
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
