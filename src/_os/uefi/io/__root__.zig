const std = @import("std");
const Io = std.Io;

pub const time = @import("time.zig");

const state = @import("state.zig");
const dir = @import("dir.zig");
const fdtable = @import("fdtable.zig");

pub const init: ?fn () void = state.init;
pub const stop: ?fn () void = state.stop;
pub const io_inst: Io = state.io_inst;
const PATH_MAX = 4096; // std.os.linux.PATH_MAX
const NAME_MAX = 255; // std.os.linux.NAME_MAX

// std.Options.OperatingSystem.Io overrides
pub const debug_io = io_inst;

pub const File = struct {
    pub const Handle = fdtable.fd_t;
    pub const Flags = struct { nonblocking: bool };
    pub const Permissions = enum(u0) {
        default_file = 0,

        pub const default_dir: @This() = .default_file;
        pub const executable_file: @This() = .default_file;
        pub const has_executable_bit = false;
    };
    pub const stdout = fdtable.stdout;
    pub const stderr = fdtable.stderr;
    pub const stdin = fdtable.stdin;
};

pub const Dir = struct {
    pub const Handle = fdtable.fd_t;
    pub const max_path_bytes = PATH_MAX;
    pub const max_name_bytes = NAME_MAX;
    pub const cwd = dir.openRootDir;

    pub const Reader = struct {
        /// Used by `dirRead`, unimplemented on this target, so
        /// this only has to be a valid size, not a meaningful one.
        pub const min_buffer_len = std.mem.alignForward(usize, NAME_MAX, @alignOf(usize));
    };
};

pub const Operation = struct {
    /// No UEFI ioctl-equivalent planned; matches `operate.operate`'s
    /// `.device_io_control => unreachable`.
    pub const DeviceIoControl = noreturn;
};
