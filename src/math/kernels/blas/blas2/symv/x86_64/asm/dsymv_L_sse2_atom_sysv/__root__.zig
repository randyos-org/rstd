//! `dsymv_L` for atom, barcelona, bobcat, dunnington, nano, penryn, piledriver, prescott (sysv), from OpenBLAS `kernel/x86_64/symv_L_sse2.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const text01 = @import("text01.zig").text;
const text02 = @import("text02.zig").text;

/// Each function's core and variant, and its parts (`__bk`'s arguments).
pub const keys = .{
    .{ .core = .atom, .variant = "sysv", .parts = "1, 0, 0, 0, 1" },
    .{ .core = .barcelona, .variant = "sysv", .parts = "0, 1, 0, 0, 0" },
    .{ .core = .bobcat, .variant = "sysv", .parts = "0, 1, 0, 0, 0" },
    .{ .core = .dunnington, .variant = "sysv", .parts = "1, 0, 0, 0, 0" },
    .{ .core = .nano, .variant = "sysv", .parts = "0, 0, 1, 0, 0" },
    .{ .core = .penryn, .variant = "sysv", .parts = "1, 0, 0, 0, 0" },
    .{ .core = .piledriver, .variant = "sysv", .parts = "0, 1, 0, 0, 0" },
    .{ .core = .prescott, .variant = "sysv", .parts = "0, 0, 0, 1, 0" },
};

/// The functions' one text, a macro whose parameters select the parts, in sections cut at its labels.
const text = text01 ++ text02;

pub fn Body(comptime args: []const u8) type {
    return struct {
        pub fn entry() callconv(.naked) void {
            asm volatile (text ++ (" __bk " ++ args ++ "\n .purgem __bk\n"));
        }
    };
}
