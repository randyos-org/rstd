//! x86_64 sum kernels. `ssum` and `dsum` are `kernel/arm/sum.c` on every core
//! (`kernel/x86_64/KERNEL`), which is not plain C: at unit stride its loop is
//! `kernel/simd/intrin.h`'s vectors, four sums of one each in the element's
//! type (`dsum` from AVX2 up); a strided `x` adds in one chain. Those slots
//! take `vector.zig`'s `sum` instead (T192, issue #545): `f32` summed in
//! `f64` and rounded once, as rstd's sums are, and four sums strided. Zen 3's
//! timings against upstream's: `tools/blasbench/results/results-sum.md`.
//! `csum` and `zsum` are `zsum_sse*.S`, which peel (`vector.peels`; `asm/`'s
//! `skipped`). Tests live in `sum/__root__.zig`.

const coretypes = @import("../../../cpu/coretypes.zig");
const registry = @import("../../../registry.zig");
const util = @import("../../../util.zig");
const vector = @import("../../../vector.zig");

const Slot = util.Slot;

pub const Core = coretypes.X86CoreType;

/// `vector.zig`'s real `sum`, four vectors a step.
pub const Vec = struct {
    pub fn entry(comptime slot: Slot) util.SlotFn(slot) {
        return switch (slot) {
            .ssum_k => &vector.Real(f32, 4).sum,
            .dsum_k => &vector.Real(f64, 4).sum,
            else => @compileError(@tagName(slot) ++ " is not a real sum slot"),
        };
    }
};

/// The ports and where they serve (`KERNEL`'s `SSUMKERNEL` and `DSUMKERNEL`,
/// which no `KERNEL.<CORE>` changes): every core.
pub const ports = [_]registry.Port(Core){
    .{ .module = Vec, .cores = &.{ .prescott, .atom, .core2, .penryn, .dunnington, .nehalem, .opteron, .opteron_sse3, .barcelona, .nano, .sandybridge, .bobcat, .bulldozer, .piledriver, .haswell, .steamroller, .excavator, .zen, .skylakex, .cooperlake, .sapphirerapids }, .slots = &.{ .ssum_k, .dsum_k }, .written = .vector },
};

const P = registry.Ports(Core, "x86_64", &ports);
const A = registry.Kernels(Core, @import("asm/__root__.zig"));

const U = registry.Union(Core, .{ A, P });
pub const hasSlot = U.hasSlot;
pub const has = U.has;
pub const skips = U.skips;
pub const exportCore = U.exportCore;
pub const install = U.install;
