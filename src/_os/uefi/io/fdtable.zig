//! Small fixed-size table mapping small integer file descriptors to real
//! open UEFI Simple File System file instances. Backs `Io.File.Handle`
//! (== `std.c.fd_t`, `i32` on this target -- see the upstream toolchain
//! patch noted in `state.zig`). fds 0/1/2 are reserved for stdin/stdout/
//! stderr (handled directly by `operate.zig`'s console path -- never
//! stored here, since they're not `uefi.protocol.File` instances); real
//! opened files get fds starting at `first_real_fd`.

const std = @import("std");
const uefi = std.os.uefi;
const Io = std.Io;

pub const stdin_fd: i32 = std.c.STDIN_FILENO;
pub const stdout_fd: i32 = std.c.STDOUT_FILENO;
pub const stderr_fd: i32 = std.c.STDERR_FILENO;
const first_real_fd: i32 = 3;

const max_open_files = 16;
var slots: [max_open_files]?*uefi.protocol.File = @splat(null);

/// Store `file` in the first free slot, returning its newly assigned fd.
/// `null` if the table is full.
pub fn alloc(file: *uefi.protocol.File) ?i32 {
    for (&slots, 0..) |*slot, i| {
        if (slot.* == null) {
            slot.* = file;
            return first_real_fd + @as(i32, @intCast(i));
        }
    }
    return null;
}

/// Look up the real file behind `fd`, if `fd` is a real (>= first_real_fd),
/// currently-open handle.
pub fn get(fd: i32) ?*uefi.protocol.File {
    if (fd < first_real_fd) return null;
    const index: usize = @intCast(fd - first_real_fd);
    if (index >= slots.len) return null;
    return slots[index];
}

/// Free `fd`'s slot, if it held one. Does not close the underlying file --
/// callers close before freeing (see `file.zig`'s `fileClose`).
pub fn free(fd: i32) void {
    if (fd < first_real_fd) return;
    const index: usize = @intCast(fd - first_real_fd);
    if (index >= slots.len) return;
    slots[index] = null;
}

/// Drop all entries without closing them -- call right after
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
