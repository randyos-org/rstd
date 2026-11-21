//! Which converted file has which kernel, its cores in its `keys`; converted from OpenBLAS `v0.3.34-794-gdbf2c24c9` by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const registry = @import("../../../../registry.zig");

pub const Entry = registry.Entry;

pub const entries = [_]Entry{
    .{ .slot = .srot_k, .module = @import("srot_lsx.zig") },
    .{ .slot = .drot_k, .module = @import("drot_lsx.zig") },
    .{ .slot = .csrot_k, .module = @import("csrot_lsx.zig") },
    .{ .slot = .zdrot_k, .module = @import("zdrot_lsx.zig") },
    .{ .slot = .srot_k, .module = @import("srot_lasx.zig") },
    .{ .slot = .drot_k, .module = @import("drot_lasx.zig") },
    .{ .slot = .csrot_k, .module = @import("csrot_lasx/__root__.zig") },
    .{ .slot = .zdrot_k, .module = @import("zdrot_lasx.zig") },
};
