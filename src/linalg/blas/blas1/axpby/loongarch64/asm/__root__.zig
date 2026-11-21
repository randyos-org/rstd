//! Which converted file has which kernel, its cores in its `keys`; converted from OpenBLAS `v0.3.34-794-gdbf2c24c9` by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const registry = @import("../../../../registry.zig");

pub const Entry = registry.Entry;

pub const entries = [_]Entry{
    .{ .slot = .saxpby_k, .module = @import("saxpby_lsx/__root__.zig") },
    .{ .slot = .daxpby_k, .module = @import("daxpby_lsx/__root__.zig") },
    .{ .slot = .caxpby_k, .module = @import("caxpby_lsx.zig") },
    .{ .slot = .zaxpby_k, .module = @import("zaxpby_lsx/__root__.zig") },
    .{ .slot = .saxpby_k, .module = @import("saxpby_lasx/__root__.zig") },
    .{ .slot = .daxpby_k, .module = @import("daxpby_lasx/__root__.zig") },
    .{ .slot = .caxpby_k, .module = @import("caxpby_lasx/__root__.zig") },
    .{ .slot = .zaxpby_k, .module = @import("zaxpby_lasx.zig") },
};
