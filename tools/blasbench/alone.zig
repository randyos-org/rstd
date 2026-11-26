//! The alone run (issue #690): the libraries again in a process where rstd
//! never runs, each on its own threads with its own defaults (OpenBLAS's
//! workers spin between calls, which the race turns off). Its OpenBLAS rates
//! are the race's `ob alone` column.

const std = @import("std");
const merge = @import("merge.zig");
const table = @import("libs/table.zig");

const Io = std.Io;
const Allocator = std.mem.Allocator;

/// OpenBLAS's rates from the alone run, by `key`; null without one.
pub var ob_rates: ?std.StringHashMapUnmanaged(f64) = null;

/// This bench again with `args` and `--accelerate-alone=1`, in the
/// environment it started with: its output echoed, OpenBLAS's rates kept.
pub fn run(io: Io, arena: Allocator, environ: *const std.process.Environ.Map, args: []const [:0]const u8) !void {
    var argv: std.ArrayList([]const u8) = .empty;
    for (args) |a| try argv.append(arena, a);
    try argv.append(arena, "--accelerate-alone=1");
    std.debug.print("\n{s}: the libraries first in a process without rstd, on their own threads\n", .{args[0]});
    var child = try std.process.spawn(io, .{ .argv = argv.items, .environ_map = environ, .stdin = .ignore, .stdout = .inherit, .stderr = .pipe });
    var text: std.ArrayList(u8) = .empty;
    var buf: [4096]u8 = undefined;
    while (true) {
        const n = child.stderr.?.readStreaming(io, &.{&buf}) catch |e| switch (e) {
            error.EndOfStream => break,
            else => |x| return x,
        };
        std.debug.print("{s}", .{buf[0..n]});
        try text.appendSlice(arena, buf[0..n]);
    }
    switch (try child.wait(io)) {
        .exited => |code| if (code != 0) return error.AloneRunFailed,
        else => return error.AloneRunFailed,
    }
    ob_rates = try rates(arena, text.items);
    std.debug.print("\nthe race: rstd with OpenBLAS's workers idle between calls (`OPENBLAS_THREAD_TIMEOUT=4` unless set, which rstd's own pool then ignores); `ob alone` is OpenBLAS's rate above\n", .{});
}

/// OpenBLAS's rate in each row of `text`'s tables, by `key`.
fn rates(arena: Allocator, text: []const u8) !std.StringHashMapUnmanaged(f64) {
    var m: std.StringHashMapUnmanaged(f64) = .empty;
    for (try merge.parse(arena, text, &table.merge_labels)) |item| switch (item) {
        .text => {},
        .table => |t| {
            const ob = for (t.header, 0..) |h, j| {
                if (std.mem.eql(u8, h, "openblas")) break j;
            } else continue;
            const first = merge.firstRate(t.header, &table.merge_labels) orelse continue;
            for (t.keys.items) |k| {
                const row = t.rows.get(k).?;
                if (row.len <= ob) continue;
                const v = merge.num(row[ob]) orelse continue;
                var kb: [512]u8 = undefined;
                const rk = try key(&kb, t.section, t.header[0..first], row);
                try m.put(arena, try arena.dupe(u8, rk), v);
            }
        },
    };
    return m;
}

/// A row's key in both runs: its section, its table's lead header, and its
/// lead cells before `rstd runs` (which the alone run prints as `-`).
pub fn key(buf: []u8, section: []const u8, lead: []const []const u8, row: []const []const u8) ![]const u8 {
    var w: Io.Writer = .fixed(buf);
    try w.print("{s}#", .{section});
    var k = lead.len;
    for (lead, 0..) |h, j| {
        if (std.mem.eql(u8, h, "rstd runs")) {
            k = j;
            break;
        }
        try w.print("{s}|", .{h});
    }
    try w.writeAll("#");
    for (row[0..@min(k, row.len)]) |c| try w.print("{s}|", .{c});
    return w.buffered();
}
