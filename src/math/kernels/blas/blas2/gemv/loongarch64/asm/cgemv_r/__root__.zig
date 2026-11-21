//! `cgemv_r` for la64_generic, from OpenBLAS `kernel/loongarch64/zgemv_n.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const text01 = @import("text01.zig").text;

/// Each function's core and variant.
pub const keys = .{
    .{ .core = .la64_generic, .variant = "" },
};

/// The function's text, in its own file.
const text = text01;

pub const Body = struct {
    pub fn entry() callconv(.naked) void {
        asm volatile (text);
    }
};
