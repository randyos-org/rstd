//! File helpers for build scripts; `src/io.zig` re-exports them.

const std = @import("std");

pub fn fileExists(local_io: std.Io, path: []const u8) bool {
    std.Io.Dir.cwd().access(local_io, path, .{}) catch {
        return false;
    };
    return true;
}
