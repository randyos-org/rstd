const sysinfo = @import("builtin");
const std = @import("std");

const impl = switch (sysinfo.os.tag) {
    .linux, .macos, .freebsd, .openbsd, .netbsd, .dragonfly => std.Io,
    .uefi => @import("_os/uefi/io/__root__.zig"),
    else => void,
};

pub const time = impl.time;
pub const ioFactory = impl.ioFactory;
/// Comptime-known `Io` value; only present for targets whose impl needs one
/// as `std_options_debug_io` (see `_os/uefi/io/__root__.zig`).
pub const io_inst: std.Io = if (@hasDecl(impl, "io_inst")) impl.io_inst else {};
pub const init: ?fn () void = if (@hasDecl(impl, "init")) impl.init else null;
pub const stop: ?fn () void = if (@hasDecl(impl, "stop")) impl.stop else null;
pub const stdout = if (@hasDecl(impl, "stdout")) impl.stdout else std.Io.File.stdout;
pub const stderr = if (@hasDecl(impl, "stderr")) impl.stderr else std.Io.File.stderr;
pub const stdin = if (@hasDecl(impl, "stdin")) impl.stdin else std.Io.File.stdin;
pub const cwdFn: ?fn () std.Io.Dir = if (@hasDecl(impl, "cwdFn")) impl.cwdFn else null;

pub const cwd = cwdFn orelse std.Io.Dir.cwd;

pub fn printout(local_io: std.Io, msg: []const u8) !void {
    try stdout().writeStreamingAll(local_io, msg);
}

pub fn printerr(local_io: std.Io, msg: []const u8) !void {
    try stderr().writeStreamingAll(local_io, msg);
}

pub fn fileExists(local_io: std.Io, path: []const u8) bool {
    std.Io.Dir.cwd().access(local_io, path, .{}) catch {
        return false;
    };
    return true;
}
