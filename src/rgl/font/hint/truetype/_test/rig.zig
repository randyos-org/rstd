//! A small engine for the interpreter's tests: its buffers, a glyph zone of
//! eight points and a twilight zone of four (fontations
//! `skrifa/src/outline/glyf/hint/engine/mod.rs` `MockEngine`).

const std = @import("std");
const Engine = @import("../Engine.zig");
const definitions = @import("../definitions.zig");
const zone = @import("../zone.zig");

const Rig = @This();

stack: [64]i32 = @splat(-7),
fdefs: [8]definitions.Definition = @splat(.{}),
idefs: [4]definitions.Definition = @splat(.{}),
cvt: [8]i32 = .{ 0, 64, 128, 192, 256, 320, 384, 448 },
cvt_copy: [8]i32 = @splat(0),
store: [8]i32 = @splat(0),
store_copy: [8]i32 = @splat(0),
orus: [8]zone.Point = undefined,
org: [8]zone.Point = undefined,
cur: [8]zone.Point = undefined,
tags: [8]u8 = @splat(zone.on_curve),
contours: [1]u16 = .{7},
tw_org: [4]zone.Point = @splat(.{ 0, 0 }),
tw_cur: [4]zone.Point = @splat(.{ 0, 0 }),
tw_tags: [4]u8 = @splat(0),
e: Engine = undefined,

/// Point `i` at (100 + 10 i, -50 i) font units, scaled by 0.5 to 26.6.
pub fn init(r: *Rig) *Engine {
    for (0..8) |i| {
        const n: i32 = @intCast(i);
        r.orus[i] = .{ 100 + 10 * n, -50 * n };
        r.org[i] = .{ 50 + 5 * n, -25 * n };
        r.cur[i] = r.org[i];
    }
    r.e = .{
        .metrics = .{ .ppem = 12, .scale = 0x8000, .point_size = 12 * 64 },
        .glyph_count = 100,
        .stack = &r.stack,
        .functions = .{ .defs = &r.fdefs },
        .instructions = .{ .defs = &r.idefs },
        .cvt = .{ .base = &r.cvt, .copy = &r.cvt_copy },
        .storage = .{ .base = &r.store, .copy = &r.store_copy },
        .zones = .{
            .{ .org = &r.tw_org, .cur = &r.tw_cur, .tags = &r.tw_tags },
            .{ .orus = &r.orus, .org = &r.org, .cur = &r.cur, .tags = &r.tags, .contours = &r.contours },
        },
    };
    return &r.e;
}

/// Runs `code` as `program`.
pub fn run(e: *Engine, program: Engine.Program, code: []const u8) Engine.Error!void {
    e.code[@backingInt(program)] = code;
    try e.run(program);
}

/// The stack after a run, bottom first.
pub fn values(e: *const Engine) []const i32 {
    return e.stack[0..e.top];
}

/// Runs `code` as the font program and expects `want` on the stack.
pub fn expectStack(e: *Engine, code: []const u8, want: []const i32) !void {
    try run(e, .font, code);
    try std.testing.expectEqualSlices(i32, want, values(e));
}
