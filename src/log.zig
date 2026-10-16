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
const root = @import("root");
const equal = @import("zigstd.zig").equal;

fn moduleName(comptime module: [:0]const u8) []const u8 {
    if (@hasDecl(root, "log_module_name") and equal(u8, module, "root"))
        return root.log_module_name;
    return module;
}

/// Like `std.log.scoped`, but the scope is the calling module's own
/// qualified name, and each call takes the caller's own `@src()` explicitly.
///
/// `Scope` must be `@This()`
pub fn Logger(comptime Scope: type) type {
    const scope_name = @typeName(Scope);
    return struct {
        /// Log an error message. This log level is intended to be used
        /// when something has gone wrong. This might be recoverable or might
        /// be followed by the program exiting.
        pub inline fn err(
            comptime src: std.lang.SourceLocation,
            comptime format: []const u8,
            args: anytype,
        ) void {
            @branchHint(.cold);
            logImpl(.err, scope_name, format, args, src);
        }

        /// Log a warning message. This log level is intended to be used if
        /// it is uncertain whether something has gone wrong or not, but the
        /// circumstances would be worth investigating.
        pub inline fn warn(
            comptime src: std.lang.SourceLocation,
            comptime format: []const u8,
            args: anytype,
        ) void {
            logImpl(.warn, scope_name, format, args, src);
        }

        /// Log an info message. This log level is intended to be used for
        /// general messages about the state of the program.
        pub inline fn info(
            comptime src: std.lang.SourceLocation,
            comptime format: []const u8,
            args: anytype,
        ) void {
            logImpl(.info, scope_name, format, args, src);
        }

        /// Log a debug message. This log level is intended to be used for
        /// messages which are only useful for debugging.
        pub inline fn debug(
            comptime src: std.lang.SourceLocation,
            comptime format: []const u8,
            args: anytype,
        ) void {
            logImpl(.debug, scope_name, format, args, src);
        }
    };
}

const log = Logger(@This());

/// platform hook for "seconds since boot" timestamp
pub var get_time: ?*const fn () f64 = null;

/// guards re-entrant logging: terminals parsing their own output (e.g.
/// framebuffer ANSI/SGR) can call back into std.log mid-line; without this,
/// first color code recurses forever.
var logging_in_progress = std.atomic.Value(bool).init(false);

/// std.log callback (via std_options.logFn). Writes through debug_io like
/// std.log.defaultLog -- shared stderr lock so lines don't tear against
/// debug.print/panic -- but with this project's line format.
pub fn logFn(
    comptime level: std.log.Level,
    comptime scope: @EnumLiteral(),
    comptime format: []const u8,
    args: anytype,
) void {
    logImpl(level, @tagName(scope), format, args, null);
}

fn levelEnabled(comptime level: std.log.Level, comptime scope_name: []const u8) bool {
    for (build_options.logger_scope_names, build_options.logger_scope_name_levels) |name, scope_level| {
        if (equal(u8, name, scope_name)) return @backingInt(level) <= scope_level;
    }
    return @backingInt(level) <= @backingInt(std.options.log_level);
}

/// Shared by `logFn` (std's own `@EnumLiteral()`-scoped entry point) and
/// `logger()`'s module-qualified wrapper type -- the only difference between
/// them is where `scope_name` comes from.
fn logImpl(
    comptime level: std.log.Level,
    comptime scope_name: []const u8,
    comptime format: []const u8,
    args: anytype,
    comptime src: ?std.lang.SourceLocation,
) void {
    // module:typename when src is available (real Logger() calls always
    // have one), else just typename (logFn's std.log.scoped(...) path,
    // which has no source location to pull a module from)
    const full_scope_name = comptime if (src) |sl| moduleName(sl.module) ++ ":" ++ scope_name else scope_name;
    if (comptime !levelEnabled(level, full_scope_name)) return;
    if (logging_in_progress.swap(true, .acq_rel)) return;
    defer logging_in_progress.store(false, .release);

    const io = std.Options.debug_io;
    const prev = io.swapCancelProtection(.blocked);
    defer _ = io.swapCancelProtection(prev);
    var buffer: [64]u8 = undefined;
    const stderr = std.debug.lockStderr(&buffer).terminal();
    defer std.debug.unlockStderr();
    // a line that cannot be written is dropped, as std's own logFn does:
    // logging must not end the program (a GUI app's console handle, a closed pipe)
    logToTerminal(stderr, level, full_scope_name, format, args, src) catch return;
}

/// Formats + writes one log line (scope filter, timestamp, colored level)
/// to any std.Io.Terminal. Shared by logFn and logToTerm; color degrades to
/// plain text when mode is .no_color.
///
/// `scope_name` here is already fully-qualified (module:typename)
pub fn logToTerminal(
    t: std.Io.Terminal,
    comptime level: std.log.Level,
    comptime scope_name: []const u8,
    comptime format: []const u8,
    args: anytype,
    comptime src: ?std.lang.SourceLocation,
    // use_color: bool,
) std.Io.Writer.Error!void {
    const color: std.Io.Terminal.Color = comptime switch (level) {
        .debug => .green,
        .info => .cyan,
        .warn => .yellow,
        .err => .red,
    };
    const time: f64 = if (get_time) |f| f() else -1;

    const scope_label = comptime if (src) |sl|
        scope_name ++ ":" ++ std.fmt.comptimePrint("{d}", .{sl.line})
    else
        scope_name;

    try t.writer.print("{d:.4} [", .{time});
    t.setColor(color) catch {};
    try t.writer.writeAll(level.asText());
    t.setColor(.reset) catch {};
    try t.writer.print("]({s}): " ++ format ++ "\n", .{scope_label} ++ args);

    // flush explicitly: buffer whole line for one puts/render, but still
    // want it on screen now -- e.g. kernel panics before idle loop
    try t.writer.flush();
}
