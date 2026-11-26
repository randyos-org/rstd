//! `caxpyc` for atom, barcelona, bobcat, core2, nano, nehalem, penryn, prescott (sysv), from OpenBLAS `kernel/x86_64/zaxpy_sse.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const text01 = @import("text01.zig").text;
const text02 = @import("text02.zig").text;
const text03 = @import("text03.zig").text;
const text04 = @import("text04.zig").text;
const text05 = @import("text05.zig").text;
const text06 = @import("text06.zig").text;
const text07 = @import("text07.zig").text;

/// Each function's core and variant, and its parts (`__bk`'s arguments).
pub const keys = .{
    .{ .core = .atom, .variant = "sysv", .parts = "1, 0, 0, 0, 0, 1" },
    .{ .core = .barcelona, .variant = "sysv", .parts = "1, 0, 1, 0, 0, 0" },
    .{ .core = .bobcat, .variant = "sysv", .parts = "1, 0, 1, 0, 0, 0" },
    .{ .core = .core2, .variant = "sysv", .parts = "1, 0, 0, 0, 0, 0" },
    .{ .core = .nano, .variant = "sysv", .parts = "1, 0, 0, 1, 0, 0" },
    .{ .core = .nehalem, .variant = "sysv", .parts = "1, 0, 0, 0, 1, 0" },
    .{ .core = .penryn, .variant = "sysv", .parts = "1, 0, 0, 0, 0, 0" },
    .{ .core = .prescott, .variant = "sysv", .parts = "0, 1, 0, 0, 0, 0" },
};

/// The functions' one text, a macro whose parameters select the parts, in sections cut at its labels.
const text = text01 ++ text02 ++ text03 ++ text04 ++ text05 ++ text06 ++ text07;

pub fn Body(comptime args: []const u8) type {
    return struct {
        pub fn entry() callconv(.naked) void {
            asm volatile (text ++ (" __bk " ++ args ++ "\n .purgem __bk\n"));
        }
    };
}
