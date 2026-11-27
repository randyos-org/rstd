//! FreeType's outlines for the hinting tests: `refs/<font>.bin.gz`, one per
//! test font (`build/testfonts.zig`), written by `tools/hintrefs/hintrefs.c`
//! from nixpkgs' FreeType (the version each file records). Regenerated with
//! `zig build hint-refs -Dtest-fonts`, never edited (issue #805).
//!
//! A file is gzipped: "RSTDHNT1"; then, as LEB128 varints (zigzag where
//! signed), FreeType's version, units per em, a CFF byte, the `fpgm`,
//! `prep` and `cvt ` lengths, the default load's mode byte and the light
//! target's (255: neither), the sizes, and the glyphs: codepoint and glyph
//! index, then for each size and mode an outline.

const std = @import("std");
const Allocator = std.mem.Allocator;

/// FreeType's load flags for each outline, as `hintrefs.c` has them.
pub const Mode = enum(u8) {
    /// `FT_LOAD_NO_HINTING`
    unhinted,
    /// `FT_LOAD_NO_AUTOHINT`: the bytecode at interpreter v40, or CFF hints
    native,
    /// `FT_LOAD_NO_AUTOHINT` at interpreter v35; TrueType only
    native_v35,
    /// `FT_LOAD_NO_AUTOHINT | FT_LOAD_TARGET_LIGHT`
    native_light,
    /// `FT_LOAD_FORCE_AUTOHINT`
    autohint,
    /// `FT_LOAD_FORCE_AUTOHINT | FT_LOAD_TARGET_LIGHT`
    autohint_light,
};

/// A `FT_CURVE_TAG` of a point.
pub const Tag = enum(u8) { conic = 0, on = 1, cubic = 2, _ };

/// One loaded glyph, in 26.6 pixels.
pub const Outline = struct {
    /// FreeType's error, 0 if the glyph loaded.
    err: i32,
    advance: i32 = 0,
    lsb_delta: i32 = 0,
    rsb_delta: i32 = 0,
    /// Each contour's last point.
    ends: []const u16 = &.{},
    tags: []const Tag = &.{},
    points: []const [2]i32 = &.{},
};

pub const Glyph = struct {
    codepoint: u21,
    /// 0 if the font has no glyph for `codepoint`; `outlines` is empty then.
    gid: u16,
    /// By size, then by mode, as `Refs.modes` lists them.
    outlines: []const Outline,
};

pub const Refs = struct {
    freetype: [3]u32,
    units_per_em: u16,
    cff: bool,
    fpgm_len: u32,
    prep_len: u32,
    cvt_len: u32,
    /// The mode FreeType's default load matched on every glyph and size.
    default_mode: ?Mode,
    /// The same for `FT_LOAD_TARGET_LIGHT`.
    default_light_mode: ?Mode,
    /// In pixels per em.
    sizes: []const u16,
    modes: []const Mode,
    glyphs: []const Glyph,

    /// `glyph`'s outline at `sizes[size]` in `mode`; null if the font has
    /// no such glyph or the file no such mode.
    pub fn outline(r: Refs, glyph: Glyph, size: usize, mode: Mode) ?Outline {
        if (glyph.gid == 0) return null;
        const m = std.mem.findScalar(Mode, r.modes, mode) orelse return null;
        return glyph.outlines[size * r.modes.len + m];
    }
};

pub const ParseError = error{ BadMagic, Truncated, BadValue } || Allocator.Error || std.compress.flate.Decompress.Error;

const ttf_modes = [_]Mode{ .unhinted, .native, .native_v35, .native_light, .autohint, .autohint_light };
const cff_modes = [_]Mode{ .unhinted, .native, .native_light, .autohint, .autohint_light };

/// Parses a `refs/<font>.bin.gz` file into `arena`.
pub fn parse(arena: Allocator, gz: []const u8) ParseError!Refs {
    var in: std.Io.Reader = .fixed(gz);
    var window: [std.compress.flate.max_window_len]u8 = undefined;
    var dec: std.compress.flate.Decompress = .init(&in, .gzip, &window);
    const raw = dec.reader.allocRemaining(arena, .unlimited) catch |e| switch (e) {
        error.OutOfMemory => return error.OutOfMemory,
        error.ReadFailed => return dec.err orelse error.Truncated,
        error.StreamTooLong => unreachable,
    };
    var c: Cursor = .{ .bytes = raw };
    if (!std.mem.eql(u8, try c.take(8), "RSTDHNT1")) return error.BadMagic;
    var r: Refs = undefined;
    for (&r.freetype) |*v| v.* = try c.int(u32);
    r.units_per_em = try c.int(u16);
    r.cff = switch (try c.byte()) {
        0 => false,
        1 => true,
        else => return error.BadValue,
    };
    r.fpgm_len = try c.int(u32);
    r.prep_len = try c.int(u32);
    r.cvt_len = try c.int(u32);
    r.default_mode = try c.mode();
    r.default_light_mode = try c.mode();
    const sizes = try arena.alloc(u16, try c.int(usize));
    for (sizes) |*s| s.* = try c.int(u16);
    r.sizes = sizes;
    r.modes = if (r.cff) &cff_modes else &ttf_modes;
    const glyphs = try arena.alloc(Glyph, try c.int(usize));
    for (glyphs) |*g| {
        g.codepoint = try c.int(u21);
        g.gid = try c.int(u16);
        if (g.gid == 0) {
            g.outlines = &.{};
            continue;
        }
        const outlines = try arena.alloc(Outline, sizes.len * r.modes.len);
        for (outlines) |*o| o.* = try c.outline(arena);
        g.outlines = outlines;
    }
    if (c.at != raw.len) return error.BadValue;
    r.glyphs = glyphs;
    return r;
}

const Cursor = struct {
    bytes: []const u8,
    at: usize = 0,

    fn take(c: *Cursor, n: usize) error{Truncated}![]const u8 {
        if (c.bytes.len - c.at < n) return error.Truncated;
        defer c.at += n;
        return c.bytes[c.at..][0..n];
    }

    fn byte(c: *Cursor) error{Truncated}!u8 {
        return (try c.take(1))[0];
    }

    fn varint(c: *Cursor) error{ Truncated, BadValue }!u64 {
        var v: u64 = 0;
        var shift: u7 = 0;
        while (true) : (shift += 7) {
            if (shift > 63) return error.BadValue;
            const b = try c.byte();
            v |= @as(u64, b & 0x7f) << @intCast(shift);
            if (b < 0x80) return v;
        }
    }

    fn int(c: *Cursor, comptime T: type) error{ Truncated, BadValue }!T {
        return std.math.cast(T, try c.varint()) orelse error.BadValue;
    }

    fn signed(c: *Cursor) error{ Truncated, BadValue }!i32 {
        const v = try c.varint();
        const s: i64 = @bitCast((v >> 1) ^ (0 -% (v & 1)));
        return std.math.cast(i32, s) orelse error.BadValue;
    }

    fn mode(c: *Cursor) error{ Truncated, BadValue }!?Mode {
        const b = try c.byte();
        if (b == 255) return null;
        return std.enums.fromInt(Mode, b) orelse error.BadValue;
    }

    fn outline(c: *Cursor, arena: Allocator) (error{ Truncated, BadValue } || Allocator.Error)!Outline {
        const err = try c.signed();
        if (err != 0) return .{ .err = err };
        var o: Outline = .{ .err = 0, .advance = try c.signed(), .lsb_delta = try c.signed(), .rsb_delta = try c.signed() };
        const ends = try arena.alloc(u16, try c.int(usize));
        var end: u32 = 0;
        for (ends, 0..) |*e, i| {
            end += try c.int(u16);
            if (i > 0 and end == ends[i - 1]) return error.BadValue;
            e.* = std.math.cast(u16, end) orelse return error.BadValue;
        }
        const n = try c.int(u16);
        const tags = try arena.alloc(Tag, n);
        for (tags, try c.take(n)) |*t, b| t.* = @fromBackingInt(b);
        const points = try arena.alloc([2]i32, n);
        var x: i32 = 0;
        var y: i32 = 0;
        for (points) |*p| {
            x +%= try c.signed();
            y +%= try c.signed();
            p.* = .{ x, y };
        }
        o.ends = ends;
        o.tags = tags;
        o.points = points;
        return o;
    }
};
