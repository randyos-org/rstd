//! File helpers for build scripts; `src/io.zig` re-exports them.

const std = @import("std");

pub fn fileExists(local_io: std.Io, path: []const u8) bool {
    std.Io.Dir.cwd().access(local_io, path, .{}) catch {
        return false;
    };
    return true;
}

/// The size of the file at `path`, or null when it cannot be opened.
pub fn fileSize(local_io: std.Io, path: []const u8) ?u64 {
    const file = std.Io.Dir.cwd().openFile(local_io, path, .{}) catch return null;
    defer file.close(local_io);
    return file.length(local_io) catch null;
}
