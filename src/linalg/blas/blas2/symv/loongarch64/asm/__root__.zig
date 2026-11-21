//! Which converted file has which kernel, its cores in its `keys`; converted from OpenBLAS `v0.3.34-794-gdbf2c24c9` by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const registry = @import("../../../../registry.zig");

pub const Entry = registry.Entry;

pub const entries = [_]Entry{
    .{ .slot = .ssymv_U, .module = @import("ssymv_U_lsx/__root__.zig") },
    .{ .slot = .ssymv_L, .module = @import("ssymv_L_lsx/__root__.zig") },
    .{ .slot = .dsymv_U, .module = @import("dsymv_U_lsx/__root__.zig") },
    .{ .slot = .dsymv_L, .module = @import("dsymv_L_lsx/__root__.zig") },
    .{ .slot = .ssymv_U, .module = @import("ssymv_U_lasx/__root__.zig") },
    .{ .slot = .ssymv_L, .module = @import("ssymv_L_lasx/__root__.zig") },
    .{ .slot = .dsymv_U, .module = @import("dsymv_U_lasx/__root__.zig") },
    .{ .slot = .dsymv_L, .module = @import("dsymv_L_lasx/__root__.zig") },
};
