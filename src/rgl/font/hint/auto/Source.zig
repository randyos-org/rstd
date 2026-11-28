//! A font's glyphs in font units, as `FT_Load_Glyph` with `FT_LOAD_NO_SCALE`
//! gives them to the autohinter, TrueType or CFF, at a variable TrueType
//! font's instance if one is set; and their advances, as `FT_Get_Advance`
//! with `FT_LOAD_NO_SCALE` gives them (`hmtx`, varied by `HVAR`, or by
//! loading the glyph without it: `src/base/ftadvanc.c`).

const std = @import("std");
const Sfnt = @import("../../sfnt/Sfnt.zig");
const truetype = @import("../truetype/__root__.zig");
const cff = @import("../cff/__root__.zig");
const vary = @import("../vary/__root__.zig");
const Glyph = @import("../Glyph.zig");

const Allocator = std.mem.Allocator;
const Source = @This();

pub const Error = truetype.Loader.Error || cff.LoadError;
pub const InitError = truetype.Instance.InitError || cff.Instance.Error;

sfnt: *const Sfnt,
/// The TrueType instance's; CFF fonts take none.
blend: ?*const vary.Blend = null,
outlines: union(enum) {
    truetype: truetype.Instance,
    cff: cff.Instance,
},

pub fn init(allocator: Allocator, sfnt: *const Sfnt, blend: ?*const vary.Blend) InitError!Source {
    if (sfnt.cff) |table| return .{ .sfnt = sfnt, .outlines = .{ .cff = .{ .sfnt = sfnt, .cff = table, .scale = 1 << 16 } } };
    return .{ .sfnt = sfnt, .blend = blend, .outlines = .{ .truetype = try .initUnscaled(allocator, sfnt, blend) } };
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

fn hmtxAdvance(s: *const Source, gid: u16) i32 {
    const hmtx = s.sfnt.hmtx orelse return 0;
    return hmtx.getMetrics(@fromBackingInt(gid)).advance_width;
}

/// `FT_Get_Advance`; null if the glyph fails to load, which leaves
/// FreeType's value as it was.
pub fn advance(s: *Source, allocator: Allocator, gid: u16) ?i32 {
    const b = s.blend orelse return s.hmtxAdvance(gid);
    if (b.hvar) |*h| return @as(u16, @truncate(@as(u32, @bitCast(s.hmtxAdvance(gid) +% h.advanceDelta(gid, b.coords)))));
    var g = s.load(allocator, gid) catch return null;
    defer g.deinit(allocator);
    return g.advance;
}

pub fn glyphCount(s: *const Source) u16 {
    return if (s.sfnt.maxp) |m| m.numGlyphs() else 0;
}

pub fn unitsPerEm(s: *const Source) u16 {
    return if (s.sfnt.head) |h| h.unitsPerEm() else 1000;
}
