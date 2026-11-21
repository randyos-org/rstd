//! Which converted file has which kernel, its cores in its `keys`; converted from OpenBLAS `v0.3.34-794-gdbf2c24c9` by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const registry = @import("../../../../registry.zig");

pub const Entry = registry.Entry;

pub const entries = [_]Entry{};

/// The kernels upstream builds here that are not converted, with the cores
/// it builds them for.
pub const skipped = [_]registry.Skipped{
    .{ .slot = .csum_k, .cores = &.{
        .atom,
        .barcelona,
        .bobcat,
        .bulldozer,
        .cooperlake,
        .core2,
        .dunnington,
        .excavator,
        .generic,
        .haswell,
        .nano,
        .nehalem,
        .opteron,
        .opteron_sse3,
        .penryn,
        .piledriver,
        .prescott,
        .sandybridge,
        .sapphirerapids,
        .skylakex,
        .steamroller,
        .zen,
    } },
    .{ .slot = .zsum_k, .cores = &.{
        .atom,
        .barcelona,
        .bobcat,
        .bulldozer,
        .cooperlake,
        .core2,
        .dunnington,
        .excavator,
        .generic,
        .haswell,
        .nano,
        .nehalem,
        .opteron,
        .opteron_sse3,
        .penryn,
        .piledriver,
        .prescott,
        .sandybridge,
        .sapphirerapids,
        .skylakex,
        .steamroller,
        .zen,
    } },
};
