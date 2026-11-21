//! Which converted file has which kernel, its cores in its `keys`; converted from OpenBLAS `v0.3.34-794-gdbf2c24c9` by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const registry = @import("../../../../registry.zig");

pub const Entry = registry.Entry;

pub const entries = [_]Entry{
    .{ .slot = .sgemv_n, .module = @import("sgemv_n/__root__.zig") },
    .{ .slot = .sgemv_t, .module = @import("sgemv_t/__root__.zig") },
    .{ .slot = .dgemv_n, .module = @import("dgemv_n/__root__.zig") },
    .{ .slot = .dgemv_t, .module = @import("dgemv_t/__root__.zig") },
    .{ .slot = .cgemv_n, .module = @import("cgemv_n/__root__.zig") },
    .{ .slot = .cgemv_t, .module = @import("cgemv_t/__root__.zig") },
    .{ .slot = .cgemv_r, .module = @import("cgemv_r/__root__.zig") },
    .{ .slot = .cgemv_c, .module = @import("cgemv_c/__root__.zig") },
    .{ .slot = .cgemv_o, .module = @import("cgemv_o/__root__.zig") },
    .{ .slot = .cgemv_u, .module = @import("cgemv_u/__root__.zig") },
    .{ .slot = .cgemv_s, .module = @import("cgemv_s/__root__.zig") },
    .{ .slot = .cgemv_d, .module = @import("cgemv_d/__root__.zig") },
    .{ .slot = .zgemv_n, .module = @import("zgemv_n/__root__.zig") },
    .{ .slot = .zgemv_t, .module = @import("zgemv_t/__root__.zig") },
    .{ .slot = .zgemv_r, .module = @import("zgemv_r/__root__.zig") },
    .{ .slot = .zgemv_c, .module = @import("zgemv_c/__root__.zig") },
    .{ .slot = .zgemv_o, .module = @import("zgemv_o/__root__.zig") },
    .{ .slot = .zgemv_u, .module = @import("zgemv_u/__root__.zig") },
    .{ .slot = .zgemv_s, .module = @import("zgemv_s/__root__.zig") },
    .{ .slot = .zgemv_d, .module = @import("zgemv_d/__root__.zig") },
};
