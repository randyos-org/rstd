//! `cgemv_n` for cell, power3, power4, power5, power6, ppc440fp2, ppc970, ppcg4 (32); cell, power3, power4, power5, power6, ppc440fp2, ppc970, ppcg4 (32le); cell, power3, power4, power5, power6, ppc970 (64); cell, power3, power4, power5, power6, ppc970 (64be), from OpenBLAS `kernel/power/zgemv_n.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const text01 = @import("text01.zig").text;
const text02 = @import("text02.zig").text;
const text03 = @import("text03.zig").text;
const text04 = @import("text04.zig").text;
const text05 = @import("text05.zig").text;
const text06 = @import("text06.zig").text;
const text07 = @import("text07.zig").text;
const text08 = @import("text08.zig").text;
const text09 = @import("text09.zig").text;
const text10 = @import("text10.zig").text;
const text11 = @import("text11.zig").text;
const text12 = @import("text12.zig").text;

/// Each function's core and variant, and its parts (`__bk`'s arguments).
pub const keys = .{
    .{ .core = .cell, .variant = "32", .parts = "1, 0, 1, 0, 1, 0, 0, 0, 1, 0" },
    .{ .core = .power3, .variant = "32", .parts = "1, 0, 0, 1, 1, 0, 0, 0, 1, 0" },
    .{ .core = .power4, .variant = "32", .parts = "1, 0, 0, 1, 1, 0, 0, 0, 0, 1" },
    .{ .core = .power5, .variant = "32", .parts = "1, 0, 0, 0, 0, 1, 0, 1, 0, 1" },
    .{ .core = .power6, .variant = "32", .parts = "1, 0, 0, 0, 0, 0, 1, 1, 0, 0" },
    .{ .core = .ppc440fp2, .variant = "32", .parts = "1, 0, 0, 1, 1, 0, 0, 0, 0, 1" },
    .{ .core = .ppc970, .variant = "32", .parts = "1, 0, 1, 0, 1, 0, 0, 0, 1, 0" },
    .{ .core = .ppcg4, .variant = "32", .parts = "1, 0, 0, 1, 1, 0, 0, 0, 1, 0" },
    .{ .core = .cell, .variant = "32le", .parts = "1, 0, 1, 0, 1, 0, 0, 0, 1, 0" },
    .{ .core = .power3, .variant = "32le", .parts = "1, 0, 0, 1, 1, 0, 0, 0, 1, 0" },
    .{ .core = .power4, .variant = "32le", .parts = "1, 0, 0, 1, 1, 0, 0, 0, 0, 1" },
    .{ .core = .power5, .variant = "32le", .parts = "1, 0, 0, 0, 0, 1, 0, 1, 0, 1" },
    .{ .core = .power6, .variant = "32le", .parts = "1, 0, 0, 0, 0, 0, 1, 1, 0, 0" },
    .{ .core = .ppc440fp2, .variant = "32le", .parts = "1, 0, 0, 1, 1, 0, 0, 0, 0, 1" },
    .{ .core = .ppc970, .variant = "32le", .parts = "1, 0, 1, 0, 1, 0, 0, 0, 1, 0" },
    .{ .core = .ppcg4, .variant = "32le", .parts = "1, 0, 0, 1, 1, 0, 0, 0, 1, 0" },
    .{ .core = .cell, .variant = "64", .parts = "0, 1, 1, 0, 1, 0, 0, 0, 1, 0" },
    .{ .core = .power3, .variant = "64", .parts = "0, 1, 0, 1, 1, 0, 0, 0, 1, 0" },
    .{ .core = .power4, .variant = "64", .parts = "0, 1, 0, 1, 1, 0, 0, 0, 0, 1" },
    .{ .core = .power5, .variant = "64", .parts = "0, 1, 0, 0, 0, 1, 0, 1, 0, 1" },
    .{ .core = .power6, .variant = "64", .parts = "0, 1, 0, 0, 0, 0, 1, 1, 0, 0" },
    .{ .core = .ppc970, .variant = "64", .parts = "0, 1, 1, 0, 1, 0, 0, 0, 1, 0" },
    .{ .core = .cell, .variant = "64be", .parts = "0, 1, 1, 0, 1, 0, 0, 0, 1, 0" },
    .{ .core = .power3, .variant = "64be", .parts = "0, 1, 0, 1, 1, 0, 0, 0, 1, 0" },
    .{ .core = .power4, .variant = "64be", .parts = "0, 1, 0, 1, 1, 0, 0, 0, 0, 1" },
    .{ .core = .power5, .variant = "64be", .parts = "0, 1, 0, 0, 0, 1, 0, 1, 0, 1" },
    .{ .core = .power6, .variant = "64be", .parts = "0, 1, 0, 0, 0, 0, 1, 1, 0, 0" },
    .{ .core = .ppc970, .variant = "64be", .parts = "0, 1, 1, 0, 1, 0, 0, 0, 1, 0" },
};

/// The functions' one text, a macro whose parameters select the parts, in sections cut at its labels.
const text = text01 ++ text02 ++ text03 ++ text04 ++ text05 ++ text06 ++ text07 ++ text08 ++ text09 ++ text10 ++ text11 ++ text12;

pub fn Body(comptime args: []const u8) type {
    return struct {
        pub fn entry() callconv(.naked) void {
            asm volatile (text ++ (" __bk " ++ args ++ "\n .purgem __bk\n"));
        }
    };
}
