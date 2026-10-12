const sysinfo = @import("builtin");
const std = @import("std");

const impl = switch (sysinfo.os.tag) {
    .uefi => @import("_os/uefi/io/__root__.zig"),
    else => void,
};

/// Comptime-known `Io` value; only present for targets whose impl needs one
/// as `std_options_debug_io` (see `_os/uefi/io/__root__.zig`). `null` for
/// impls that don't provide one, so callers fall back to std's own default.
pub const io_inst: ?std.Io = if (@hasDecl(impl, "io_inst")) impl.io_inst else null;
pub const init: ?fn () void = if (@hasDecl(impl, "init")) impl.init else null;
pub const stop: ?fn () void = if (@hasDecl(impl, "stop")) impl.stop else null;
pub const stdout = std.Io.File.stdout;
pub const stderr = std.Io.File.stderr;
pub const stdin = std.Io.File.stdin;
pub const cwd = std.Io.Dir.cwd;

pub fn printout(local_io: std.Io, msg: []const u8) !void {
    try stdout().writeStreamingAll(local_io, msg);
}

pub fn printerr(local_io: std.Io, msg: []const u8) !void {
    try stderr().writeStreamingAll(local_io, msg);
}

pub const fileExists = @import("buildutils/io.zig").fileExists;

pub fn getTimeSeconds() f64 {
    const ns = std.Io.Clock.now(.boot, std.Options.debug_io).toNanoseconds();
    if (ns == std.math.maxInt(@TypeOf(ns))) return -1;
    return @as(f64, @floatFromInt(ns)) / std.time.ns_per_s;
}
