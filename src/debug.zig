const sysinfo = @import("builtin");
const std = @import("std");

pub const debug = sysinfo.optimize == .debug;

pub fn assertMsg(condition: bool, msg: []const u8) void {
    if (debug and !condition) {
        std.debug.panic("assertion failed: {s}", .{msg});
    }
}
