//! `zgemv_u` for la464, from OpenBLAS `kernel/loongarch64/zgemv_t_4_lasx.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const text01 = @import("text01.zig").text;
const text02 = @import("text02.zig").text;

/// Each function's core and variant.
pub const keys = .{
    .{ .core = .la464, .variant = "" },
};

/// The function's text, in sections cut at its labels.
const text = text01 ++ text02;

pub const Body = struct {
    pub fn entry() callconv(.naked) void {
        asm volatile (text);
    }
};
