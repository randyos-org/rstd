//! x86_64 `symv` kernels: the converted `.S` kernels (`asm/`), and the
//! `@Vector` ports of `kernel/x86_64/ssymv_L.c`, `ssymv_U.c`, `dsymv_L.c` and
//! `dsymv_U.c` (`../driver.zig`) with the arithmetic of each core's
//! `*_microk_<core>-2.c` block. `dsymv_L_microk_skylakex-2.c` is intrinsics,
//! the others inline assembly; Skylake-X's other three take Haswell's.
//! Tests live in `symv/__root__.zig`. The AVX-512 and FMA4 cores' kernels
//! run only on such hardware (T43).

const coretypes = @import("../../../cpu/coretypes.zig");
const registry = @import("../../../registry.zig");
const fma = @import("../../../fma.zig");
const util = @import("../../../util.zig");
const kernel = @import("../driver.zig");

const Slot = util.Slot;
const Shape = kernel.Shape;

pub const Core = coretypes.X86CoreType;

/// Which microkernels a core's drivers include.
pub const Micro = enum { nehalem, sandy, haswell, skylakex, bulldozer };

/// The block shape of `slot` in `micro`'s microkernel, from its assembly or
/// intrinsics.
fn shape(comptime micro: Micro, comptime slot: Slot) Shape {
    return switch (micro) {
        // SSE, unfused; the U kernels store their sums
        .nehalem => switch (slot) {
            .dsymv_L => .{ .lanes = 2, .fused = false },
            .dsymv_U => .{ .lanes = 2, .fused = false, .store = true },
            .ssymv_L => .{ .lanes = 4, .fused = false },
            .ssymv_U => .{ .lanes = 4, .fused = false, .store = true },
            else => unreachable,
        },
        // AVX, unfused
        .sandy => switch (slot) {
            .dsymv_L, .dsymv_U => .{ .lanes = 4, .fused = false },
            .ssymv_L => .{ .lanes = 8, .fused = false, .narrow_odd = true },
            .ssymv_U => .{ .lanes = 8, .fused = false, .store = true },
            else => unreachable,
        },
        // FMA3; `ssymv_L` in 128 bits
        .haswell => switch (slot) {
            .dsymv_L, .dsymv_U, .ssymv_L => .{ .lanes = 4, .fused = true },
            .ssymv_U => .{ .lanes = 8, .fused = true, .store = true },
            else => unreachable,
        },
        // `dsymv_L`'s intrinsics, contracted as the core's flags say
        .skylakex => switch (slot) {
            .dsymv_L => .{ .lanes = 4, .fused = fma.fuse, .sum_y = true, .wide = true },
            else => shape(.haswell, slot),
        },
        // FMA4, 128 bits
        .bulldozer => switch (slot) {
            .dsymv_L => .{ .lanes = 2, .fused = true },
            .dsymv_U => .{ .lanes = 2, .fused = true, .store = true },
            .ssymv_L => .{ .lanes = 4, .fused = true },
            .ssymv_U => .{ .lanes = 4, .fused = true, .store = true },
            else => unreachable,
        },
    };
}

/// The ports of one microkernel family.
pub fn Ported(comptime micro: Micro) type {
    return struct {
        pub fn entry(comptime slot: Slot) util.SlotFn(slot) {
            return kernel.x86(slot, shape(micro, slot));
        }
    };
}

/// The ports and where they serve (`KERNEL.<CORE>`'s `?SYMV_?_KERNEL`).
pub const ports = [_]registry.Port(Core){
    .{ .module = Ported(.nehalem), .cores = &.{.nehalem}, .slots = &.{ .ssymv_U, .ssymv_L, .dsymv_U, .dsymv_L }, .written = .vector },
    .{ .module = Ported(.sandy), .cores = &.{.sandybridge}, .slots = &.{ .ssymv_U, .ssymv_L, .dsymv_U, .dsymv_L }, .written = .vector },
    .{ .module = Ported(.haswell), .cores = &.{ .haswell, .zen }, .slots = &.{ .ssymv_U, .ssymv_L, .dsymv_U, .dsymv_L }, .written = .vector },
    .{ .module = Ported(.skylakex), .cores = &.{ .skylakex, .cooperlake, .sapphirerapids }, .slots = &.{ .ssymv_U, .ssymv_L, .dsymv_U, .dsymv_L }, .written = .vector },
    .{ .module = Ported(.bulldozer), .cores = &.{ .bulldozer, .steamroller, .excavator }, .slots = &.{ .ssymv_U, .ssymv_L, .dsymv_U, .dsymv_L }, .written = .vector },
};

const P = registry.Ports(Core, "x86_64", &ports);
const A = registry.Kernels(Core, @import("asm/__root__.zig"));

const U = registry.Union(Core, .{ A, P });
pub const hasSlot = U.hasSlot;
pub const has = U.has;
pub const skips = U.skips;
pub const exportCore = U.exportCore;
pub const install = U.install;
