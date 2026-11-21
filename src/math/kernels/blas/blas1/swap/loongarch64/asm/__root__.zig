//! Which converted file has which kernel, its cores in its `keys`; converted from OpenBLAS `v0.3.34-794-gdbf2c24c9` by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const registry = @import("../../../../registry.zig");

pub const Entry = registry.Entry;

pub const entries = [_]Entry{
    .{ .slot = .sswap_k, .module = @import("sswap_lsx.zig") },
    .{ .slot = .dswap_k, .module = @import("dswap_lsx.zig") },
    .{ .slot = .cswap_k, .module = @import("cswap_lsx.zig") },
    .{ .slot = .zswap_k, .module = @import("zswap_lsx.zig") },
    .{ .slot = .sswap_k, .module = @import("sswap_lasx.zig") },
    .{ .slot = .dswap_k, .module = @import("dswap_lasx.zig") },
    .{ .slot = .cswap_k, .module = @import("cswap_lasx.zig") },
    .{ .slot = .zswap_k, .module = @import("zswap_lasx.zig") },
    .{ .slot = .sswap_k, .module = @import("sswap.zig") },
    .{ .slot = .dswap_k, .module = @import("dswap.zig") },
};
