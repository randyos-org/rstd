//! `load_window <name> <base> <size> <image> <stamp> [memory]`: fails when an ELF image's `PT_LOAD` span, and its `__stack_top` where
//! the linker script puts the boot stack past the image, is not inside the window a rig's firmware maps (or with `memory`, the rig's
//! RAM: issue #575); `<stamp>` gets the span, so the build caches a pass.

const std = @import("std");
const Io = std.Io;
const elf = std.elf;

pub fn main(init: std.process.Init) !void {
    const arena = init.arena.allocator();
    const io = init.io;
    const args = try init.minimal.args.toSlice(arena);
    if (args.len != 6 and args.len != 7) {
        std.debug.print("usage: load_window <name> <base> <size> <image> <stamp> [memory]\n", .{});
        std.process.exit(2);
    }
    const name = args[1];
    const base = try std.fmt.parseInt(u64, args[2], 0);
    const size = try std.fmt.parseInt(u64, args[3], 0);
    const memory = args.len == 7;
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
    if (try symbol(header, data, "__stack_top")) |top| hi = @max(hi, top);
    if (lo < base or hi > base + size) {
        if (memory) {
            std.debug.print("{s}: the image and its stack end at 0x{x} ({d} MiB), past the rig's memory, 0x{x} at 0x{x}; raise `memory` in build/qemu/machines.zig\n", .{ name, hi, (hi - base) >> 20, size, base });
        } else {
            std.debug.print("{s}: loaded span 0x{x} (0x{x}..0x{x}) is past the firmware's window, 0x{x} at 0x{x}; split its part in build/qemu/machines.zig\n", .{ name, hi - lo, lo, hi, size, base });
        }
        std.process.exit(1);
    }
    const file = try Io.Dir.cwd().createFile(io, args[5], .{});
    defer file.close(io);
    var buf: [32]u8 = undefined;
    try file.writeStreamingAll(io, try std.mem.print(&buf, "0x{x}\n", .{hi - lo}));
}

/// The value of the symbol `want` in the image's `.symtab`, read from its
/// bytes (`Elf32_Sym`, `Elf64_Sym`); null without it.
fn symbol(header: elf.Header, data: []const u8, want: []const u8) !?u64 {
    var symtab: ?elf.Elf64_Shdr = null;
    var it = header.iterateSectionHeadersBuffer(data);
    while (try it.next()) |sh| {
        if (sh.sh_type == elf.SHT_SYMTAB) symtab = sh;
    }
    const st = symtab orelse return null;
    // the symbols' names: the section `sh_link` names
    var strtab: ?elf.Elf64_Shdr = null;
    var index: usize = 0;
    var it2 = header.iterateSectionHeadersBuffer(data);
    while (try it2.next()) |sh| : (index += 1) {
        if (index == st.sh_link) strtab = sh;
    }
    const names = strtab orelse return null;
    const entsize: usize = if (header.is_64) 24 else 16;
    for (0..st.sh_size / entsize) |i| {
        const at = st.sh_offset + i * entsize;
        if (at + entsize > data.len) return error.EndOfStream;
        const e = data[at..][0..entsize];
        const start = names.sh_offset + std.mem.readInt(u32, e[0..4], header.endian);
        if (start >= data.len) continue;
        const rest = data[start..];
        const len = std.mem.findScalar(u8, rest, 0) orelse continue;
        if (!std.mem.eql(u8, rest[0..len], want)) continue;
        return if (header.is_64) std.mem.readInt(u64, e[8..16], header.endian) else std.mem.readInt(u32, e[4..8], header.endian);
    }
    return null;
}
