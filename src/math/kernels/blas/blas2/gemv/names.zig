//! The gemv blocks' names in a core's `uk/` module and the driver's own, as
//! upstream's microkernels name them. A form's `n` and `t` blocks live in
//! separate modules under the same names, so they share an entry.

const util = @import("../../util.zig");

const Slot = util.Slot;

/// `k4x8` is `sgemv_n`'s alone, `add_y` the complex `n` microkernel's.
pub const Names = struct {
    k4x8: ?[:0]const u8 = null,
    k4x4: [:0]const u8,
    k4x2: [:0]const u8,
    k4x1: [:0]const u8,
    add_y: ?[:0]const u8 = null,
};

/// The block names of `slot`, one of `driver.slots`.
pub fn names(comptime slot: Slot) Names {
    return switch (slot) {
        .sgemv_n, .sgemv_t => .{ .k4x8 = "sgemv_kernel_4x8", .k4x4 = "sgemv_kernel_4x4", .k4x2 = "sgemv_kernel_4x2", .k4x1 = "sgemv_kernel_4x1" },
        .dgemv_n, .dgemv_t => .{ .k4x4 = "dgemv_kernel_4x4", .k4x2 = "dgemv_kernel_4x2", .k4x1 = "dgemv_kernel_4x1" },
        .cgemv_n, .cgemv_t => .{ .k4x4 = "cgemv_kernel_4x4", .k4x2 = "cgemv_kernel_4x2", .k4x1 = "cgemv_kernel_4x1", .add_y = "add_y_asm" },
        .cgemv_r, .cgemv_c => .{ .k4x4 = "cgemv_kernel_4x4_c", .k4x2 = "cgemv_kernel_4x2_c", .k4x1 = "cgemv_kernel_4x1_c", .add_y = "add_y_c_asm" },
        .cgemv_o, .cgemv_u => .{ .k4x4 = "cgemv_kernel_4x4_x", .k4x2 = "cgemv_kernel_4x2_x", .k4x1 = "cgemv_kernel_4x1_x", .add_y = "add_y_x_asm" },
        .cgemv_s, .cgemv_d => .{ .k4x4 = "cgemv_kernel_4x4_cx", .k4x2 = "cgemv_kernel_4x2_cx", .k4x1 = "cgemv_kernel_4x1_cx", .add_y = "add_y_cx_asm" },
        .zgemv_n, .zgemv_t => .{ .k4x4 = "zgemv_kernel_4x4", .k4x2 = "zgemv_kernel_4x2", .k4x1 = "zgemv_kernel_4x1", .add_y = "add_y_asm" },
        .zgemv_r, .zgemv_c => .{ .k4x4 = "zgemv_kernel_4x4_c", .k4x2 = "zgemv_kernel_4x2_c", .k4x1 = "zgemv_kernel_4x1_c", .add_y = "add_y_c_asm" },
        .zgemv_o, .zgemv_u => .{ .k4x4 = "zgemv_kernel_4x4_x", .k4x2 = "zgemv_kernel_4x2_x", .k4x1 = "zgemv_kernel_4x1_x", .add_y = "add_y_x_asm" },
        .zgemv_s, .zgemv_d => .{ .k4x4 = "zgemv_kernel_4x4_cx", .k4x2 = "zgemv_kernel_4x2_cx", .k4x1 = "zgemv_kernel_4x1_cx", .add_y = "add_y_cx_asm" },
        else => @compileError(@tagName(slot) ++ " is not a gemv kernel"),
    };
}
