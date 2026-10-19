//! Writes a Linux initramfs (a `newc` cpio archive) for the Linux qemu rigs:
//! the directories the rig's init mounts over, `/dev/console` (so the kernel
//! can give init its stdio), and the named executables.
//!
//!   initramfs <out.cpio> <name>=<file>...
//!
//! puts each `<file>` at `/<name>`, mode 0755.

const std = @import("std");
const Io = std.Io;

pub fn main(init: std.process.Init) !void {
    const arena = init.arena.allocator();
    const io = init.io;
    const args = try init.minimal.args.toSlice(arena);
    if (args.len < 3) {
        std.debug.print("usage: initramfs <out.cpio> <name>=<file>...\n", .{});
        std.process.exit(2);
    }
    var out: std.ArrayList(u8) = .empty;
    var ino: u32 = 1;
    for ([_][]const u8{ "dev", "proc", "sys", "tmp" }) |d| try entry(arena, &out, &ino, d, 0o040755, 0, 0, "");
    // c 5 1, as the kernel's own default initramfs has it
    try entry(arena, &out, &ino, "dev/console", 0o020600, 5, 1, "");
    for (args[2..]) |arg| {
        const eq = std.mem.indexOfScalar(u8, arg, '=') orelse return error.BadArgument;
        const data = try Io.Dir.cwd().readFileAlloc(io, arg[eq + 1 ..], arena, .unlimited);
        try entry(arena, &out, &ino, arg[0..eq], 0o100755, 0, 0, data);
    }
    try entry(arena, &out, &ino, "TRAILER!!!", 0, 0, 0, "");
    const file = try Io.Dir.cwd().createFile(io, args[1], .{});
    defer file.close(io);
    try file.writeStreamingAll(io, out.items);
}

/// One `newc` record: header, name, data, each padded to 4 bytes.
fn entry(gpa: std.mem.Allocator, out: *std.ArrayList(u8), ino: *u32, name: []const u8, mode: u32, rmajor: u32, rminor: u32, data: []const u8) !void {
    const fields = [_]u32{
        ino.*, mode, 0, 0, if (mode & 0o170000 == 0o040000) 2 else 1, 0, @intCast(data.len), 0, 0, rmajor, rminor, @intCast(name.len + 1), 0,
    };
    ino.* += 1;
    try out.appendSlice(gpa, "070701");
    for (fields) |f| try out.print(gpa, "{x:0>8}", .{f});
    try out.appendSlice(gpa, name);
    try out.append(gpa, 0);
    try pad(gpa, out);
    try out.appendSlice(gpa, data);
    try pad(gpa, out);
}

fn pad(gpa: std.mem.Allocator, out: *std.ArrayList(u8)) !void {
    while (out.items.len % 4 != 0) try out.append(gpa, 0);
}
