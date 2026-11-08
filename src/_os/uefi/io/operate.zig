const std = @import("std");
const Io = std.Io;

const state = @import("state.zig");
const console = @import("console.zig");
const fdtable = @import("fdtable.zig");

/// Backs `Io.File.writeStreamingAll`/`.readStreamingAll`; the streaming
/// counterpart to `dir.zig`/`file.zig`'s positional ops. stdin/stdout/stderr
/// (fds 0/1/2, see `fdtable.zig`) go through the console; real opened files
/// (fd >= 3) go through `fdtable.get`. Device ioctl and networking have no
/// UEFI backing (networking planned, ioctl not).
pub fn operate(userdata: ?*anyopaque, operation: Io.Operation) Io.Cancelable!Io.Operation.Result {
    _ = userdata;
    return switch (operation) {
        .file_write_streaming => |o| .{ .file_write_streaming = writeStreaming(o) },
        .file_read_streaming => |o| .{ .file_read_streaming = readStreaming(o) },
        .device_io_control => unreachable, // no UEFI ioctl-equivalent planned
        .net_receive => .{
            // TODO: back with EFI network protocols (see stubs.zig)
            .net_receive = .{ error.NetworkDown, 0 },
        },
        .net_send => .{
            // TODO: back with EFI network protocols (see stubs.zig)
            .net_send = .{ error.NetworkDown, 0 },
        },
        .net_accept => .{
            // TODO: back with EFI network protocols (see stubs.zig)
            .net_accept = error.NetworkDown,
        },
        .net_read => .{
            // TODO: back with EFI network protocols (see stubs.zig)
            .net_read = error.NetworkDown,
        },
        .net_write => .{
            // TODO: back with EFI network protocols (see stubs.zig)
            .net_write = error.NetworkDown,
        },
    };
}

fn writeStreaming(o: Io.Operation.FileWriteStreaming) Io.Operation.FileWriteStreaming.Error!usize {
    if (o.file.handle == fdtable.stdout_fd or o.file.handle == fdtable.stderr_fd) {
        // only one console tracked today (state.con_out); stdout and stderr
        // both land there until a separate std_err protocol is wired up
        return console.writeSplat(state.con_out, o.header, o.data, o.splat);
    }

    const f = fdtable.get(o.file.handle) orelse return error.NotOpenForWriting;
    var n: usize = 0;
    if (o.header.len != 0) n += f.write(o.header) catch return error.InputOutput;
    if (o.data.len != 0) {
        for (o.data[0 .. o.data.len - 1]) |chunk| {
            n += f.write(chunk) catch return error.InputOutput;
        }
        const last = o.data[o.data.len - 1];
        for (0..o.splat) |_| n += f.write(last) catch return error.InputOutput;
    }
    return n;
}

fn readStreaming(o: Io.Operation.FileReadStreaming) Io.Operation.FileReadStreaming.Error!usize {
    if (o.file.handle == fdtable.stdin_fd) {
        // TODO: no ConIn (Simple Text Input Protocol) backing yet
        return error.InputOutput;
    }

    const f = fdtable.get(o.file.handle) orelse return error.NotOpenForReading;
    return fdtable.readInto(f, o.data);
}
