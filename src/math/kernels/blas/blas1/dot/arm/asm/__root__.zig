//! Which converted file has which kernel, its cores in its `keys`; converted from OpenBLAS `v0.3.34-794-gdbf2c24c9` by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const registry = @import("../../../../registry.zig");

pub const Entry = registry.Entry;

pub const entries = [_]Entry{
    .{ .slot = .ddot_k, .module = @import("ddot_vfp.zig") },
    .{ .slot = .zdotu_k, .module = @import("zdotu_vfp.zig") },
    .{ .slot = .zdotc_k, .module = @import("zdotc_vfp.zig") },
    .{ .slot = .cdotu_k, .module = @import("cdotu_vfp.zig") },
    .{ .slot = .cdotc_k, .module = @import("cdotc_vfp.zig") },
};
