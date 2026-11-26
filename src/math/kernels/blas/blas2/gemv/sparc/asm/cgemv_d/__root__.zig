//! `cgemv_d` for sparc, sparcv7 (32); sparc, sparcv7 (64), from OpenBLAS `kernel/sparc/zgemv_t.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const text01 = @import("text01.zig").text;
const text02 = @import("text02.zig").text;
const text03 = @import("text03.zig").text;
const text04 = @import("text04.zig").text;

/// Each function's core and variant, and its parts (`__bk`'s arguments).
pub const keys = .{
    .{ .core = .sparc, .variant = "32", .parts = "1, 0" },
    .{ .core = .sparcv7, .variant = "32", .parts = "1, 0" },
    .{ .core = .sparc, .variant = "64", .parts = "0, 1" },
    .{ .core = .sparcv7, .variant = "64", .parts = "0, 1" },
};

/// The functions' one text, a macro whose parameters select the parts, in sections cut at its labels.
const text = text01 ++ text02 ++ text03 ++ text04;

pub fn Body(comptime args: []const u8) type {
    return struct {
        pub fn entry() callconv(.naked) void {
            asm volatile (text ++ (" __bk " ++ args ++ "\n .purgem __bk\n"));
        }
    };
}
