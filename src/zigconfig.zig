const std = @import("std");
const io = @import("io.zig");

pub const os = struct {
    pub const PATH_MAX = io.PATH_MAX;
    pub const NAME_MAX = io.NAME_MAX;
    pub const heap = struct {
        // pub const page_allocator: std.mem.Allocator = std.heap.page_allocator;
    };
};

const default_opts = std.Options{};
pub const std_options_debug_io = io.io_inst orelse default_opts.debug_io;
pub const std_options_cwd: ?fn () std.Io.Dir = io.cwd;
