//! A font's glyphs in font units, as `FT_Load_Glyph` with `FT_LOAD_NO_SCALE`
//! gives them to the autohinter, TrueType or CFF; and their advances, as
//! `FT_Get_Advance` with `FT_LOAD_NO_SCALE` (`hmtx`).

const std = @import("std");
const Sfnt = @import("../../sfnt/Sfnt.zig");
const truetype = @import("../truetype/__root__.zig");
const cff = @import("../cff/__root__.zig");
const Glyph = @import("../Glyph.zig");

const Allocator = std.mem.Allocator;
const Source = @This();

pub const Error = truetype.Loader.Error || cff.LoadError;
pub const InitError = truetype.Instance.InitError || cff.Instance.Error;

sfnt: *const Sfnt,
outlines: union(enum) {
    truetype: truetype.Instance,
    cff: cff.Instance,
},

pub fn init(allocator: Allocator, sfnt: *const Sfnt) InitError!Source {
    if (sfnt.cff) |table| return .{ .sfnt = sfnt, .outlines = .{ .cff = .{ .sfnt = sfnt, .cff = table, .scale = 1 << 16 } } };
    return .{ .sfnt = sfnt, .outlines = .{ .truetype = try .initUnscaled(allocator, sfnt) } };
}

pub fn deinit(s: *Source) void {
    switch (s.outlines) {
        .truetype => |*i| i.deinit(),
        .cff => {},
    }
}

/// The glyph's outline in font units; the caller owns it.
pub fn load(s: *Source, allocator: Allocator, gid: u16) Error!Glyph {
    return switch (s.outlines) {
        .truetype => |*i| truetype.Loader.load(allocator, i, gid, false),
        .cff => |*i| cff.load(allocator, i, @fromBackingInt(gid), false),
    };
}

pub fn advance(s: *const Source, gid: u16) i32 {
    const hmtx = s.sfnt.hmtx orelse return 0;
    return hmtx.getMetrics(@fromBackingInt(gid)).advance_width;
}

pub fn glyphCount(s: *const Source) u16 {
    return if (s.sfnt.maxp) |m| m.numGlyphs() else 0;
}

pub fn unitsPerEm(s: *const Source) u16 {
    return if (s.sfnt.head) |h| h.unitsPerEm() else 1000;
}
