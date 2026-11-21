//! Which converted file has which kernel, its cores in its `keys`; converted from OpenBLAS `v0.3.34-794-gdbf2c24c9` by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const registry = @import("../../../../registry.zig");

pub const Entry = registry.Entry;

pub const entries = [_]Entry{
    .{ .slot = .saxpy_k, .module = @import("saxpy_loongson3a.zig") },
    .{ .slot = .daxpy_k, .module = @import("daxpy_loongson3a_simd.zig") },
};
