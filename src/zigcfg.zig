const sysinfo = @import("builtin");
const std = @import("std");
const fma = @import("linalg/blas/fma.zig");

// const io = @import("io.zig");

pub const std_os_options: std.Options.OperatingSystem = .{
    .Io = switch (sysinfo.target.os.tag) {
        .uefi => @import("_os/uefi/io/__root__.zig"),
        else => null,
    },
};

// const default_opts = std.Options{};
// pub const std_options_debug_io = default_opts.debug_io;
// pub const std_options_cwd: ?fn () std.Io.Dir = io.cwd;
pub const log_module_name = "root";

/// A root's `std_float_options` type, the same declaration for either compiler:
/// `pub const std_float_options: rstd.zigcfg.FloatOptions = .{ .fuse_mul_add_bits = 0 };`.
/// `std.Options.Float` where std has it (this project's fork), else a copy std
/// ignores and rstd reads (`float_options`). Defined in `fma.zig`.
pub const FloatOptions = fma.FloatOptions;

/// Whether the compiler reads `std_float_options` itself, for `@Matrix`.
pub const has_float_options = fma.has_float_options;

/// The build's float options, the root's or the defaults.
pub const float_options = fma.float_options;
