//! Project-wide log format: std.log callback prefixing each line with
//! seconds-since-boot + colored level/scope, e.g. `12.3456 [info](kmain): ...`.
//! That's why this exists instead of std.log.defaultLog.
//!
//! Goes through std.Options.debug_io (same stderr lock as std.debug.print/
//! panic), like std's default logger -- needs no per-platform wiring beyond
//! the get_time hook. logToTerm writes directly to an rstd Terminal instead,
//! bypassing debug_io (e.g. mirror log lines to a 2nd console).

// const sysinfo = @import("builtin");
const build_options = @import("build_options");

const std = @import("std");
const log = std.log.scoped(.logging);

const Terminal = @import("Terminal.zig");
// const ansi = @import("ansi.zig");

/// default logging terminal for logToTerm callers; logFn doesn't use this,
/// writes to std.Options.debug_io instead
pub var log_term: ?*Terminal = null;

/// platform hook for "seconds since boot" timestamp
pub var get_time: ?*const fn () f64 = null;

/// fmt: timestamp, start_color, level, end_color, scope, then caller fmt/args
// const log_prefix_fmt = "{d:.4} [{s}{s}{s}]({s}): ";

/// guards re-entrant logging: terminals parsing their own output (e.g.
/// framebuffer ANSI/SGR) can call back into std.log mid-line; without this,
/// first color code recurses forever.
var logging_in_progress: bool = false;

/// std.log callback (via std_options.logFn). Writes through debug_io like
/// std.log.defaultLog -- shared stderr lock so lines don't tear against
/// debug.print/panic -- but with this project's line format.
pub fn logFn(
    comptime level: std.log.Level,
    comptime scope: @TypeOf(.EnumLiteral),
    comptime format: []const u8,
    args: anytype,
) void {
    if (logging_in_progress) return;
    logging_in_progress = true;
    defer logging_in_progress = false;

    const io = std.Options.debug_io;
    const prev = io.swapCancelProtection(.blocked);
    defer _ = io.swapCancelProtection(prev);
    var buffer: [64]u8 = undefined;
    const stderr = std.debug.lockStderr(&buffer).terminal();
    defer std.debug.unlockStderr();
    logToTerminal(stderr, level, @tagName(scope), format, args, true) catch @panic("log print failed");
}

/// Formats + writes one log line (scope filter, timestamp, colored level)
/// to any std.Io.Terminal. Shared by logFn and logToTerm; color degrades to
/// plain text when mode is .no_color.
pub fn logToTerminal(
    t: std.Io.Terminal,
    comptime level: std.log.Level,
    scope_name: []const u8,
    comptime format: []const u8,
    args: anytype,
    use_color: bool,
) std.Io.Writer.Error!void {
    for (build_options.logger_scopes_ignore) |ignore| {
        if (std.mem.eql(u8, ignore, scope_name)) {
            return;
        }
    }
    const color: std.Io.Terminal.Color = comptime switch (level) {
        .debug => .green,
        .info => .cyan,
        .warn => .yellow,
        .err => .red,
    };
    const time: f64 = if (get_time) |f| f() else -1;

    try t.writer.print("{d:.4} [", .{time});
    if (use_color) t.setColor(color) catch {};
    try t.writer.writeAll(level.asText());
    if (use_color) t.setColor(.reset) catch {};
    try t.writer.print("]({s}): " ++ format ++ "\n", .{scope_name} ++ args);

    // flush explicitly: buffer whole line for one puts/render, but still
    // want it on screen now -- e.g. kernel panics before idle loop
    t.writer.flush() catch @panic("log flush failed");
}
