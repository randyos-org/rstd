//! `load_window <name> <base> <size> <image> <stamp>`: fails when an ELF image's `PT_LOAD` span is not inside the window a rig's firmware
//! maps; `<stamp>` gets the span, so the build caches a pass.

const std = @import("std");
const Io = std.Io;
const elf = std.elf;

pub fn main(init: std.process.Init) !void {
    const arena = init.arena.allocator();
    const io = init.io;
    const args = try init.minimal.args.toSlice(arena);
    if (args.len != 6) {
        std.debug.print("usage: load_window <name> <base> <size> <image> <stamp>\n", .{});
        std.process.exit(2);
    }
    const name = args[1];
    const base = try std.fmt.parseInt(u64, args[2], 0);
    const size = try std.fmt.parseInt(u64, args[3], 0);
    const data = try Io.Dir.cwd().readFileAlloc(io, args[4], arena, .unlimited);
    var reader: Io.Reader = .fixed(data);
    const header = try elf.Header.read(&reader);
    var lo: u64 = std.math.maxInt(u64);
    var hi: u64 = 0;
    var it = header.iterateProgramHeadersBuffer(data);
    while (try it.next()) |ph| {
        if (ph.type != .LOAD) continue;
        lo = @min(lo, ph.vaddr);
        hi = @max(hi, ph.vaddr + ph.memsz);
    }
    if (hi == 0) return error.NoLoadSegments;
    if (lo < base or hi > base + size) {
        std.debug.print("{s}: loaded span 0x{x} (0x{x}..0x{x}) is past the firmware's window, 0x{x} at 0x{x}; split its part in build/qemu/machines.zig\n", .{ name, hi - lo, lo, hi, size, base });
        std.process.exit(1);
    }
    const file = try Io.Dir.cwd().createFile(io, args[5], .{});
    defer file.close(io);
    var buf: [32]u8 = undefined;
    try file.writeStreamingAll(io, try std.fmt.bufPrint(&buf, "0x{x}\n", .{hi - lo}));
}
