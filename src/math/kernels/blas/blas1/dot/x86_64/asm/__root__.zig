//! Which converted file has which kernel, its cores in its `keys`; converted from OpenBLAS `v0.3.34-794-gdbf2c24c9` by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const registry = @import("../../../../registry.zig");

pub const Entry = registry.Entry;

pub const entries = [_]Entry{};

/// The kernels upstream builds here that are not converted, with the cores
/// it builds them for.
pub const skipped = [_]registry.Skipped{
    .{ .slot = .ddot_k, .cores = &.{
        .atom,
        .barcelona,
        .bobcat,
        .bulldozer,
        .core2,
        .dunnington,
        .nano,
        .opteron,
        .opteron_sse3,
        .penryn,
        .prescott,
    } },
    .{ .slot = .zdotu_k, .cores = &.{
        .atom,
        .barcelona,
        .bobcat,
        .core2,
        .dunnington,
        .nano,
        .nehalem,
        .opteron,
        .opteron_sse3,
        .penryn,
        .prescott,
    } },
    .{ .slot = .zdotc_k, .cores = &.{
        .atom,
        .barcelona,
        .bobcat,
        .core2,
        .dunnington,
        .nano,
        .nehalem,
        .opteron,
        .opteron_sse3,
        .penryn,
        .prescott,
    } },
    .{ .slot = .cdotu_k, .cores = &.{
        .atom,
        .barcelona,
        .bobcat,
        .core2,
        .dunnington,
        .nano,
        .nehalem,
        .opteron,
        .opteron_sse3,
        .penryn,
        .prescott,
    } },
    .{ .slot = .cdotc_k, .cores = &.{
        .atom,
        .barcelona,
        .bobcat,
        .core2,
        .dunnington,
        .nano,
        .nehalem,
        .opteron,
        .opteron_sse3,
        .penryn,
        .prescott,
    } },
};
