//! Which converted file has which kernel, its cores in its `keys`; converted from OpenBLAS `v0.3.34-794-gdbf2c24c9` by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const registry = @import("../../../../registry.zig");

pub const Entry = registry.Entry;

pub const entries = [_]Entry{
    .{ .slot = .sscal_k, .module = @import("sscal.zig") },
    .{ .slot = .dscal_k, .module = @import("dscal.zig") },
    .{ .slot = .cscal_k, .module = @import("cscal.zig") },
    .{ .slot = .zscal_k, .module = @import("zscal.zig") },
};
