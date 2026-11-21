//! `dsymv_L` for core2 (sysv), from OpenBLAS `kernel/x86_64/symv_L_sse2.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const text01 = @import("text01.zig").text;
const text02 = @import("text02.zig").text;

/// Each function's core and variant.
pub const keys = .{
    .{ .core = .core2, .variant = "sysv" },
};

/// The function's text, in sections cut at its labels.
const text = text01 ++ text02;

pub const Body = struct {
    pub fn entry() callconv(.naked) void {
        asm volatile (text);
    }
};
