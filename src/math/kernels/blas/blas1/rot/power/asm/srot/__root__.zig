//! `srot` for cell, power3, power4, power5, power6 (32); cell, power3, power4, power5, power6 (32le); cell, power3, power4, power5, power6 (64); cell, power3, power4, power5, power6 (64be), from OpenBLAS `kernel/power/rot.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const text01 = @import("text01.zig").text;

/// Each function's core and variant, and its parts (`__bk`'s arguments).
pub const keys = .{
    .{ .core = .cell, .variant = "32", .parts = "1, 0, 1, 0, 0, 1" },
    .{ .core = .power3, .variant = "32", .parts = "1, 0, 0, 1, 0, 1" },
    .{ .core = .power4, .variant = "32", .parts = "1, 0, 0, 1, 0, 1" },
    .{ .core = .power5, .variant = "32", .parts = "1, 0, 0, 1, 0, 1" },
    .{ .core = .power6, .variant = "32", .parts = "1, 0, 0, 0, 1, 0" },
    .{ .core = .cell, .variant = "32le", .parts = "1, 0, 1, 0, 0, 1" },
    .{ .core = .power3, .variant = "32le", .parts = "1, 0, 0, 1, 0, 1" },
    .{ .core = .power4, .variant = "32le", .parts = "1, 0, 0, 1, 0, 1" },
    .{ .core = .power5, .variant = "32le", .parts = "1, 0, 0, 1, 0, 1" },
    .{ .core = .power6, .variant = "32le", .parts = "1, 0, 0, 0, 1, 0" },
    .{ .core = .cell, .variant = "64", .parts = "0, 1, 1, 0, 0, 1" },
    .{ .core = .power3, .variant = "64", .parts = "0, 1, 0, 1, 0, 1" },
    .{ .core = .power4, .variant = "64", .parts = "0, 1, 0, 1, 0, 1" },
    .{ .core = .power5, .variant = "64", .parts = "0, 1, 0, 1, 0, 1" },
    .{ .core = .power6, .variant = "64", .parts = "0, 1, 0, 0, 1, 0" },
    .{ .core = .cell, .variant = "64be", .parts = "0, 1, 1, 0, 0, 1" },
    .{ .core = .power3, .variant = "64be", .parts = "0, 1, 0, 1, 0, 1" },
    .{ .core = .power4, .variant = "64be", .parts = "0, 1, 0, 1, 0, 1" },
    .{ .core = .power5, .variant = "64be", .parts = "0, 1, 0, 1, 0, 1" },
    .{ .core = .power6, .variant = "64be", .parts = "0, 1, 0, 0, 1, 0" },
};

/// The functions' one text, a macro whose parameters select the parts, in its own file.
const text = text01;

pub fn Body(comptime args: []const u8) type {
    return struct {
        pub fn entry() callconv(.naked) void {
            asm volatile (text ++ (" __bk " ++ args ++ "\n .purgem __bk\n"));
        }
    };
}
