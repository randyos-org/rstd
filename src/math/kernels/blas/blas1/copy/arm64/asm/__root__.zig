//! Which converted file has which kernel, its cores in its `keys`; converted from OpenBLAS `v0.3.34-794-gdbf2c24c9` by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const registry = @import("../../../../registry.zig");

pub const Entry = registry.Entry;

pub const entries = [_]Entry{
    .{ .slot = .scopy_k, .module = @import("scopy.zig") },
    .{ .slot = .dcopy_k, .module = @import("dcopy.zig") },
    .{ .slot = .ccopy_k, .module = @import("ccopy.zig") },
    .{ .slot = .zcopy_k, .module = @import("zcopy.zig") },
};
