const sysinfo = @import("builtin");
const std = @import("std");

// const io = @import("io.zig");

pub const std_os_options: std.Options.OperatingSystem = .{
    .Io = switch (sysinfo.os.tag) {
        .uefi => @import("_os/uefi/io/__root__.zig"),
        else => null,
    },
};

// const default_opts = std.Options{};
// pub const std_options_debug_io = default_opts.debug_io;
// pub const std_options_cwd: ?fn () std.Io.Dir = io.cwd;
pub const log_module_name = "root";
