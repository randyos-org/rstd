//! `zsymv_U` for atom, barcelona, bobcat, bulldozer, cooperlake, dunnington, excavator, haswell, nano, nehalem, penryn, piledriver, prescott, sandybridge, sapphirerapids, skylakex, steamroller, zen (win), from OpenBLAS `kernel/x86_64/zsymv_U_sse2.S` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const text01 = @import("text01.zig").text;
const text02 = @import("text02.zig").text;

/// Each function's core and variant, and its parts (`__bk`'s arguments).
pub const keys = .{
    .{ .core = .atom, .variant = "win", .parts = "1, 0, 0, 0" },
    .{ .core = .barcelona, .variant = "win", .parts = "0, 1, 0, 0" },
    .{ .core = .bobcat, .variant = "win", .parts = "0, 1, 0, 0" },
    .{ .core = .bulldozer, .variant = "win", .parts = "0, 1, 0, 0" },
    .{ .core = .cooperlake, .variant = "win", .parts = "1, 0, 0, 0" },
    .{ .core = .dunnington, .variant = "win", .parts = "1, 0, 0, 0" },
    .{ .core = .excavator, .variant = "win", .parts = "0, 1, 0, 0" },
    .{ .core = .haswell, .variant = "win", .parts = "1, 0, 0, 0" },
    .{ .core = .nano, .variant = "win", .parts = "0, 0, 1, 0" },
    .{ .core = .nehalem, .variant = "win", .parts = "1, 0, 0, 0" },
    .{ .core = .penryn, .variant = "win", .parts = "1, 0, 0, 0" },
    .{ .core = .piledriver, .variant = "win", .parts = "0, 1, 0, 0" },
    .{ .core = .prescott, .variant = "win", .parts = "0, 0, 0, 1" },
    .{ .core = .sandybridge, .variant = "win", .parts = "1, 0, 0, 0" },
    .{ .core = .sapphirerapids, .variant = "win", .parts = "1, 0, 0, 0" },
    .{ .core = .skylakex, .variant = "win", .parts = "1, 0, 0, 0" },
    .{ .core = .steamroller, .variant = "win", .parts = "0, 1, 0, 0" },
    .{ .core = .zen, .variant = "win", .parts = "1, 0, 0, 0" },
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
