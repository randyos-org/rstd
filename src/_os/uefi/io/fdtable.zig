//! A fixed table from small fds to open UEFI Simple File System files; backs
//! `Io.File.Handle` (`i32`). fds 0, 1 and 2 are the console (`operate.zig`),
//! never stored here; opened files start at `first_real_fd`.

const std = @import("std");
const uefi = std.os.uefi;
const Io = std.Io;

pub const fd_t = i32;
pub const no_handle: fd_t = -1;
pub const stdin_fd: fd_t = std.c.STDIN_FILENO;
pub const stdout_fd: fd_t = std.c.STDOUT_FILENO;
pub const stderr_fd: fd_t = std.c.STDERR_FILENO;
const first_real_fd: fd_t = 3;

const max_open_files = 16;
var slots: [max_open_files]?*uefi.protocol.File = @splat(null);

/// Store `file` in the first free slot, returning its newly assigned fd.
/// `null` if the table is full.
pub fn alloc(file: *uefi.protocol.File) ?fd_t {
    for (&slots, 0..) |*slot, i| {
        if (slot.* == null) {
            slot.* = file;
            return first_real_fd + @as(fd_t, @intCast(i));
        }
    }
    return null;
}

/// Look up the real file behind `fd`, if `fd` is a real (>= first_real_fd),
/// currently-open handle.
pub fn get(fd: fd_t) ?*uefi.protocol.File {
    return slots[slotIndex(fd) orelse return null];
}

/// Free `fd`'s slot, if it held one. Does not close the underlying file;
/// callers close before freeing (see `file.zig`'s `fileClose`).
pub fn free(fd: fd_t) void {
    slots[slotIndex(fd) orelse return] = null;
}

/// `fd`'s slot, if `fd` is a real one in the table's range.
fn slotIndex(fd: fd_t) ?usize {
    if (fd < first_real_fd) return null;
    const index: usize = @intCast(@as(fd_t, fd) - @as(fd_t, first_real_fd));
    if (index >= slots.len) return null;
    return index;
}

/// Reads `f` into each of `data` in turn, up to a short read (EOF); the bytes read.
pub fn readInto(f: *uefi.protocol.File, data: []const []u8) error{InputOutput}!usize {
    var total: usize = 0;
    for (data) |buffer| {
        if (buffer.len == 0) continue;
        const n = f.read(buffer) catch return error.InputOutput;
        total += n;
        if (n < buffer.len) break;
    }
    return total;
}

/// Drop all entries without closing them; call right after
/// `exitBootServices()`, where the backing protocols are gone regardless
/// (see `state.stop()`).
pub fn reset() void {
    slots = @splat(null);
}

pub fn stdout() Io.File {
    return .{
        .handle = stdout_fd,
        .flags = .{ .nonblocking = false },
    };
}

pub fn stderr() Io.File {
    return .{
        .handle = stderr_fd,
        .flags = .{ .nonblocking = false },
    };
}

pub fn stdin() Io.File {
    return .{
        .handle = stdin_fd,
        .flags = .{ .nonblocking = false },
    };
}
