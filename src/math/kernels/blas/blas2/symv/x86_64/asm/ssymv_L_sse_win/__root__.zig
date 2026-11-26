//! `ssymv_L` for atom, barcelona, bobcat, core2, dunnington, nano, opteron, opteron_sse3, penryn, piledriver, prescott (win), from OpenBLAS `kernel/x86_64/symv_L_sse.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const text01 = @import("text01.zig").text;
const text02 = @import("text02.zig").text;
const text03 = @import("text03.zig").text;

/// Each function's core and variant, and its parts (`__bk`'s arguments).
pub const keys = .{
    .{ .core = .atom, .variant = "win", .parts = "1, 0, 1, 0, 0, 0, 1, 1, 0" },
    .{ .core = .barcelona, .variant = "win", .parts = "1, 0, 0, 1, 0, 0, 0, 1, 0" },
    .{ .core = .bobcat, .variant = "win", .parts = "1, 0, 0, 1, 0, 0, 0, 1, 0" },
    .{ .core = .core2, .variant = "win", .parts = "1, 0, 1, 0, 0, 0, 0, 1, 0" },
    .{ .core = .dunnington, .variant = "win", .parts = "1, 0, 1, 0, 0, 0, 0, 1, 0" },
    .{ .core = .nano, .variant = "win", .parts = "1, 0, 0, 0, 1, 0, 0, 1, 0" },
    .{ .core = .opteron, .variant = "win", .parts = "0, 1, 0, 0, 0, 0, 0, 0, 1" },
    .{ .core = .opteron_sse3, .variant = "win", .parts = "0, 1, 0, 0, 0, 0, 0, 1, 0" },
    .{ .core = .penryn, .variant = "win", .parts = "1, 0, 1, 0, 0, 0, 0, 1, 0" },
    .{ .core = .piledriver, .variant = "win", .parts = "1, 0, 0, 1, 0, 0, 0, 1, 0" },
    .{ .core = .prescott, .variant = "win", .parts = "1, 0, 0, 0, 0, 1, 0, 1, 0" },
};

/// The functions' one text, a macro whose parameters select the parts, in sections cut at its labels.
const text = text01 ++ text02 ++ text03;

pub fn Body(comptime args: []const u8) type {
    return struct {
        pub fn entry() callconv(.naked) void {
            asm volatile (text ++ (" __bk " ++ args ++ "\n .purgem __bk\n"));
        }
    };
}
