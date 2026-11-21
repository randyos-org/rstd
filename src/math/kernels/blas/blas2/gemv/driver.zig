//! The drivers of `kernel/x86_64/{s,d,c,z}gemv_n_4.c` and `gemv_t_4.c` around
//! the blocks `uk/index_n.zig` and `uk/index_t.zig` name (microkernel, the
//! driver's assembly, else its C, as `HAVE_KERNEL_*` picks). The C contracts
//! by the core's flags (`fma.madd`); callers pass `m`, `n` of at least 1. No
//! BLAS-root imports (`_kernel.zig`).

const std = @import("std");
const sysinfo = @import("builtin");
const util = @import("../../util.zig");
const types = @import("../../types.zig");
const fma = @import("../../fma.zig");
const uk = @import("../../uk_block.zig");
const coretypes = @import("../../cpu/coretypes.zig");
const real_driver = @import("driver_real.zig");
const complex_driver = @import("driver_complex.zig");

const BlasLong = types.BlasLong;
const Slot = util.Slot;
const Core = coretypes.X86CoreType;

/// The kernels an architecture's `<arch>/kernel.zig` iterates over.
pub const slots = [_]Slot{
    .sgemv_n, .sgemv_t, .dgemv_n, .dgemv_t,
    .cgemv_n, .cgemv_t, .cgemv_r, .cgemv_c,
    .cgemv_o, .cgemv_u, .cgemv_s, .cgemv_d,
    .zgemv_n, .zgemv_t, .zgemv_r, .zgemv_c,
    .zgemv_o, .zgemv_u, .zgemv_s, .zgemv_d,
};

/// The element type of `slot`, one of the real `slots`.
pub fn Elem(comptime slot: Slot) type {
    return switch (slot) {
        .sgemv_n, .sgemv_t => f32,
        .dgemv_n, .dgemv_t => f64,
        else => @compileError(@tagName(slot) ++ " is not a real gemv kernel"),
    };
}

/// Vtable slot type of `slot`.
pub fn Kernel(comptime slot: Slot) type {
    return util.SlotFn(slot);
}

/// `NBMAX`, the rows a block takes: `sgemv` halves it on Steamroller and
/// Excavator.
pub fn nbmax(comptime slot: Slot, comptime core: Core) BlasLong {
    if (Elem(slot) == f64) return 2048;
    return switch (core) {
        .steamroller, .excavator => 2048,
        else => 4096,
    };
}

/// `p` moved by `k` elements, either way.
pub inline fn off(comptime F: type, p: [*]F, k: BlasLong) [*]F {
    return @ptrCast(util.at(F, p, k));
}

/// The block types, as the C declares them.
pub fn N4x8(comptime F: type) type {
    return *const fn (BlasLong, [*]const [*]F, [*]F, [*]F, BlasLong, [*]const F) void;
}
pub fn N4(comptime F: type) type {
    return *const fn (BlasLong, [*]const [*]F, [*]F, [*]F, [*]const F) void;
}
pub fn N1(comptime F: type) type {
    return *const fn (BlasLong, [*]F, [*]F, [*]F, [*]const F) void;
}
pub fn T4(comptime F: type) type {
    return *const fn (BlasLong, [*]const [*]F, [*]F, [*]F) void;
}
pub fn T2(comptime F: type) type {
    return *const fn (BlasLong, [*]F, [*]F, [*]F, [*]F) void;
}
pub fn T1(comptime F: type) type {
    return *const fn (BlasLong, [*]F, [*]F, [*]F) void;
}

/// `<p>gemv_kernel_<shape>` from the microkernel `M`, else from the driver's
/// own assembly `O`, else `fallback`, the driver's C.
pub fn choose(comptime Fn: type, comptime M: type, comptime O: type, comptime name: []const u8, comptime fallback: ?Fn) Fn {
    if (@hasDecl(M, name)) return &@field(M, name);
    if (O != void and @hasDecl(O, name)) return &@field(O, name);
    return fallback orelse @compileError("no block " ++ name);
}

pub const Fallback = real_driver.Fallback;
const GemvN = real_driver.GemvN;
const GemvT = real_driver.GemvT;
const ZgemvN = complex_driver.ZgemvN;
const ZgemvT = complex_driver.ZgemvT;

/// The kernel `slot` of `core`, whose index entry is `e` (`uk_block.Kernels`).
pub fn entryOf(comptime slot: Slot, comptime core: Core, comptime e: uk.Entry, comptime Vector: fn (comptime Slot) type) Kernel(slot) {
    const M = switch (e.block) {
        .assembly => e.module,
        .vector => Vector(slot),
        .none => struct {},
    };
    return switch (slot) {
        .sgemv_n, .dgemv_n => &GemvN(slot, core, M, e.own, fma.fuse, void).call,
        .sgemv_t, .dgemv_t => &GemvT(slot, core, M, e.own, fma.fuse, void).call,
        .cgemv_n, .cgemv_r, .cgemv_o, .cgemv_s, .zgemv_n, .zgemv_r, .zgemv_o, .zgemv_s => &ZgemvN(slot, M, fma.fuse).call,
        .cgemv_t, .cgemv_c, .cgemv_u, .cgemv_d, .zgemv_t, .zgemv_c, .zgemv_u, .zgemv_d => &ZgemvT(slot, M, fma.fuse).call,
        else => @compileError(@tagName(slot) ++ " is not driven here"),
    };
}

/// Skylake-X's `sgemv_n`: `sgemv_n_4.c` around Haswell's block `e` (its
/// index entry) and `S`, the intrinsics block (`sgemv_n_microk_skylakex-8.c`).
pub fn skylakeN(comptime e: uk.Entry, comptime S: type) Kernel(.sgemv_n) {
    return &GemvN(.sgemv_n, .skylakex, e.module, e.own, fma.fuse, S).call;
}

/// Skylake-X's `sgemv_t`: `sgemv_t_4.c` around Haswell's block `e` and `S`
/// (`sgemv_t_microk_skylakex.c`).
pub fn skylakeT(comptime e: uk.Entry, comptime S: type) Kernel(.sgemv_t) {
    return &GemvT(.sgemv_t, .skylakex, e.module, e.own, fma.fuse, S).call;
}
