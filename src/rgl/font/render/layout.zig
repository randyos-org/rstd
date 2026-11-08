const std = @import("std");
const zigstd = @import("../../../zigstd.zig");
const Allocator = zigstd.Allocator;

const Sfnt = @import("../sfnt/Sfnt.zig");
const types = @import("../sfnt/types.zig");
const GlyphIndex = types.GlyphIndex;
const test_font = @import("_test/font.zig");

/// One glyph placed during layout: pixel bounds already include the glyph's
/// own bearing and the pen position at the time it was laid out, so these
/// can be blitted directly once shifted into an atlas/screen origin.
pub const PositionedGlyph = struct {
    glyph: GlyphIndex,
    codepoint: u21,
    pixel_x1: i32,
    pixel_x2: i32,
    pixel_y1: i32,
    pixel_y2: i32,
};

/// The result of laying a run of text out, optionally wrapped across
/// multiple lines. Caller owns `glyphs`.
pub const TextLayout = struct {
    glyphs: []PositionedGlyph,
    min_x: i32 = 0,
    max_x: i32 = 0,
    min_y: i32 = 0,
    max_y: i32 = 0,

    pub const empty: TextLayout = .{ .glyphs = &.{} };

    pub fn deinit(self: *TextLayout, gpa: Allocator) void {
        gpa.free(self.glyphs);
        self.* = undefined;
    }

    pub fn width(self: TextLayout) u32 {
        return @intCast(self.max_x - self.min_x);
    }

    pub fn height(self: TextLayout) u32 {
        return @intCast(self.max_y - self.min_y);
    }
};

pub const Options = struct {
    /// null: lay everything out on one unbounded line (matches
    /// buildTextAtlas's historical behavior). Otherwise, wrap at the last
    /// word boundary once a line's bounding box would reach this pixel
    /// width, force-breaking mid-word if a single word is wider than this
    /// on its own.
    wrap_width_px: ?u32 = null,
};

pub const LayoutError = Allocator.Error || error{ InvalidUtf8, MissingHead };

const Box = struct {
    min_x: i32 = 0,
    max_x: i32 = 0,
    min_y: i32 = 0,
    max_y: i32 = 0,

    fn width(self: Box) i32 {
        return self.max_x - self.min_x;
    }

    fn merge(a: Box, b: Box) Box {
        return .{
            .min_x = @min(a.min_x, b.min_x),
            .max_x = @max(a.max_x, b.max_x),
            .min_y = @min(a.min_y, b.min_y),
            .max_y = @max(a.max_y, b.max_y),
        };
    }
};

fn pxFromFunit(funit: i64, scale: f32) i32 {
    const funit_f: f32 = @floatFromInt(funit);
    return @round(funit_f * scale);
}

const LayoutState = enum { in_word, between_word };

const RollbackData = struct {
    byte_pos: usize = 0,
    glyphs_len: usize = 0,
    bounds: Box = .{},
    start_x: i64 = 0,
};

const Helper = struct {
    gpa: Allocator,
    sfnt: Sfnt,
    scale: f32,
    line_height_funit: i32,
    wrap_width_px: ?u32,
    utf8: std.unicode.Utf8Iterator,

    cursor_x: i64 = 0,
    cursor_y: i64 = 0,
    last_char_byte_pos: usize = 0,
    /// Last glyph actually placed, for kerning the next one against it.
    /// Reset to null on every line break; kerning never crosses one.
    prev_glyph: ?GlyphIndex = null,
    glyphs: std.ArrayList(PositionedGlyph) = .empty,
    bounds: Box,
    layout_state: LayoutState = .between_word,
    rollback_data: RollbackData = .{},

    fn init(gpa: Allocator, sfnt: Sfnt, text: []const u8, units_per_em: u16, pixel_size: f32, wrap_width_px: ?u32) LayoutError!Helper {
        const view = std.unicode.Utf8View.init(text) catch return error.InvalidUtf8;
        const scale = pixel_size / @as(f32, @floatFromInt(units_per_em));
        const vmetrics = sfnt.verticalMetrics();

        return .{
            .gpa = gpa,
            .sfnt = sfnt,
            .scale = scale,
            .line_height_funit = vmetrics.line_height,
            .wrap_width_px = wrap_width_px,
            .utf8 = view.iterator(),
            .bounds = .{
                .min_y = pxFromFunit(vmetrics.descender, scale),
                .max_y = pxFromFunit(vmetrics.ascender, scale),
            },
        };
    }

    fn deinit(self: *Helper) void {
        self.glyphs.deinit(self.gpa);
    }

    fn step(self: *Helper) !bool {
        const start_byte = self.utf8.i;
        const cp = self.utf8.nextCodepoint() orelse return false;
        self.last_char_byte_pos = start_byte;
        self.updateRollbackData(cp, start_byte);

        if (cp == '\n') {
            self.advanceLine();
            return true;
        }

        const glyph = self.sfnt.getGlyphIndex(cp);
        if (glyph) |g| {
            if (self.prev_glyph) |prev| self.cursor_x += self.sfnt.glyphKernAdvance(prev, g);
        }

        const glyph_bounds = if (glyph) |g| self.calcCharBounds(g) else null;

        if (glyph_bounds == null) {
            const advance = if (glyph) |g| self.sfnt.getGlyphMetrics(g).advance_width else 0;
            self.advanceNoGlyphChar(advance);
            self.prev_glyph = glyph;
            return true;
        }

        const new_bounds = self.bounds.merge(glyph_bounds.?);
        if (self.wrap_width_px) |wrap_width| {
            if (new_bounds.width() >= wrap_width) {
                self.doTextWrapping();
                return true;
            }
        }

        self.cursor_x += self.sfnt.getGlyphMetrics(glyph.?).advance_width;
        self.prev_glyph = glyph;

        try self.glyphs.append(self.gpa, .{
            .glyph = glyph.?,
            .codepoint = cp,
            .pixel_x1 = glyph_bounds.?.min_x,
            .pixel_x2 = glyph_bounds.?.max_x,
            .pixel_y1 = glyph_bounds.?.min_y,
            .pixel_y2 = glyph_bounds.?.max_y,
        });

        self.bounds = new_bounds;
        return true;
    }

    fn advanceLine(self: *Helper) void {
        self.cursor_y -= self.line_height_funit;
        self.cursor_x = 0;
        self.bounds.min_y -= pxFromFunit(self.line_height_funit, self.scale);
        // Moving up a line puts us back at the start of the line, not
        // wherever the current word started.
        self.rollback_data.start_x = 0;
        // Kerning never crosses a line break.
        self.prev_glyph = null;
    }

    fn updateRollbackData(self: *Helper, cp: u21, start_byte: usize) void {
        const is_whitespace = cp < 128 and std.ascii.isWhitespace(@intCast(cp));
        if (is_whitespace) {
            self.layout_state = .between_word;
            return;
        }
        if (self.layout_state == .in_word) return;

        self.layout_state = .in_word;
        self.rollback_data = .{
            .byte_pos = start_byte,
            .glyphs_len = self.glyphs.items.len,
            .bounds = self.bounds,
            .start_x = self.cursor_x,
        };
    }

    fn advanceNoGlyphChar(self: *Helper, advance_width: u16) void {
        self.cursor_x += advance_width;
        self.bounds.max_x += pxFromFunit(advance_width, self.scale);
        if (self.wrap_width_px) |wrap_width| {
            // -1 to stay BELOW the wrap width, or future checks get
            // confused about why the bounding box is already >= it.
            self.bounds.max_x = @min(@as(i32, @intCast(wrap_width)) - 1, self.bounds.max_x);
        }
    }

    fn doTextWrapping(self: *Helper) void {
        const word_at_line_start = self.rollback_data.start_x == 0;
        if (word_at_line_start) {
            // The word itself is wider than the wrap width; force a
            // mid-word break by re-decoding the current character on the
            // next line instead of looping on it forever.
            self.utf8.i = self.last_char_byte_pos;
        } else {
            self.rollback();
        }
        self.advanceLine();
    }

    fn rollback(self: *Helper) void {
        self.utf8.i = self.rollback_data.byte_pos;
        self.glyphs.shrinkRetainingCapacity(self.rollback_data.glyphs_len);
        self.bounds = self.rollback_data.bounds;
        self.cursor_x = self.rollback_data.start_x;
    }

    fn calcCharBounds(self: *Helper, glyph: GlyphIndex) ?Box {
        const bbox = self.sfnt.getGlyphBoundingBox(glyph) orelse return null;

        const x1_funit = self.cursor_x + bbox.x_min;
        const x2_funit = x1_funit + (bbox.x_max - bbox.x_min);
        const y1_funit = self.cursor_y + bbox.y_min;
        const y2_funit = y1_funit + (bbox.y_max - bbox.y_min);

        const x1_px = pxFromFunit(x1_funit, self.scale);
        const y1_px = pxFromFunit(y1_funit, self.scale);
        // Convert x1 and the width separately (not x2 directly) so a
        // glyph's pixel size stays consistent regardless of how far the
        // cursor has drifted from the origin.
        const x2_px = x1_px + pxFromFunit(x2_funit - x1_funit, self.scale);
        const y2_px = y1_px + pxFromFunit(y2_funit - y1_funit, self.scale);

        return .{ .min_x = x1_px, .max_x = x2_px, .min_y = y1_px, .max_y = y2_px };
    }
};

/// Lays `text` out left-to-right (and top-to-bottom across lines, when
/// wrapping), returning every glyph's pixel placement plus the overall
/// bounding box. With `options.wrap_width_px == null` this reproduces the
/// single-line behavior `render/atlas.zig`'s buildTextAtlas has always had;
/// otherwise it wraps at word boundaries, force-breaking a single word that
/// is wider than the wrap width on its own.
pub fn layoutText(gpa: Allocator, sfnt: Sfnt, text: []const u8, pixel_size: f32, options: Options) LayoutError!TextLayout {
    const head = sfnt.head orelse return error.MissingHead;

    var helper = try Helper.init(gpa, sfnt, text, head.unitsPerEm(), pixel_size, options.wrap_width_px);
    errdefer helper.deinit();

    while (try helper.step()) {}

    return .{
        .glyphs = try helper.glyphs.toOwnedSlice(gpa),
        .min_x = helper.bounds.min_x,
        .max_x = helper.bounds.max_x,
        .min_y = helper.bounds.min_y,
        .max_y = helper.bounds.max_y,
    };
}

test "layoutText lays out a single unwrapped line matching pen advance + glyph bounds" {
    const alloc = zigstd.testAllocator;
    const font_bytes = try test_font.build(alloc);
    defer alloc.free(font_bytes);

    const sfnt = try Sfnt.init(alloc, font_bytes, .{ .check = true });

    var layout = try layoutText(alloc, sfnt, "AB", 1000, .{});
    defer layout.deinit(alloc);

    // scale = 1000/1000 = 1.0: glyph box is 0..200 in x, glyph advance 250.
    try zigstd.testExpectEqual(@as(usize, 2), layout.glyphs.len);
    try zigstd.testExpectEqual(@as(i32, 0), layout.glyphs[0].pixel_x1);
    try zigstd.testExpectEqual(@as(i32, 200), layout.glyphs[0].pixel_x2);
    try zigstd.testExpectEqual(@as(i32, 250), layout.glyphs[1].pixel_x1);
    try zigstd.testExpectEqual(@as(i32, 450), layout.glyphs[1].pixel_x2);
}

test "layoutText wraps at a word boundary, rolling the second word onto a new line" {
    const alloc = zigstd.testAllocator;
    const font_bytes = try test_font.build(alloc);
    defer alloc.free(font_bytes);

    const sfnt = try Sfnt.init(alloc, font_bytes, .{ .check = true });

    // "AB CD": word "AB" spans x in [0,450), word "CD" would push the line
    // to [500,950) if kept on the same line. Wrap width 500 forces "CD"
    // onto its own line.
    var layout = try layoutText(alloc, sfnt, "AB CD", 1000, .{ .wrap_width_px = 500 });
    defer layout.deinit(alloc);

    try zigstd.testExpectEqual(@as(usize, 4), layout.glyphs.len);
    // "AB" stays on line 0 (y unshifted).
    try zigstd.testExpectEqual(@as(i32, 0), layout.glyphs[0].pixel_y1);
    try zigstd.testExpectEqual(@as(i32, 0), layout.glyphs[1].pixel_y1);
    // "CD" rolled onto line 1: x restarts at 0, y is shifted down by one line.
    try zigstd.testExpectEqual(@as(i32, 0), layout.glyphs[2].pixel_x1);
    try zigstd.testExpectEqual(@as(i32, 200), layout.glyphs[2].pixel_x2);
    try zigstd.testExpectTrue(layout.glyphs[2].pixel_y1 < layout.glyphs[0].pixel_y1);
    try zigstd.testExpectEqual(layout.glyphs[2].pixel_y1, layout.glyphs[3].pixel_y1);
}

test "layoutText force-breaks a single word wider than the wrap width, without looping forever" {
    const alloc = zigstd.testAllocator;
    const font_bytes = try test_font.build(alloc);
    defer alloc.free(font_bytes);

    const sfnt = try Sfnt.init(alloc, font_bytes, .{ .check = true });

    // Every letter is its own "word" character; each glyph alone (0..200)
    // is already narrower than wrap_width_px=150 is not possible (glyph
    // itself would never fit), so use a wrap width that fits exactly one
    // glyph (250, the advance) to force a break after every letter.
    var layout = try layoutText(alloc, sfnt, "AAAA", 1000, .{ .wrap_width_px = 250 });
    defer layout.deinit(alloc);

    try zigstd.testExpectEqual(@as(usize, 4), layout.glyphs.len);
    // Every glyph landed at the start of its own line.
    for (layout.glyphs) |g| {
        try zigstd.testExpectEqual(@as(i32, 0), g.pixel_x1);
    }
    // And each successive glyph is on a strictly lower line than the last.
    try zigstd.testExpectTrue(layout.glyphs[1].pixel_y1 < layout.glyphs[0].pixel_y1);
    try zigstd.testExpectTrue(layout.glyphs[2].pixel_y1 < layout.glyphs[1].pixel_y1);
    try zigstd.testExpectTrue(layout.glyphs[3].pixel_y1 < layout.glyphs[2].pixel_y1);
}

test "layoutText starts a new line on an explicit newline" {
    const alloc = zigstd.testAllocator;
    const font_bytes = try test_font.build(alloc);
    defer alloc.free(font_bytes);

    const sfnt = try Sfnt.init(alloc, font_bytes, .{ .check = true });

    var layout = try layoutText(alloc, sfnt, "A\nB", 1000, .{});
    defer layout.deinit(alloc);

    try zigstd.testExpectEqual(@as(usize, 2), layout.glyphs.len);
    try zigstd.testExpectEqual(@as(i32, 0), layout.glyphs[1].pixel_x1);
    try zigstd.testExpectTrue(layout.glyphs[1].pixel_y1 < layout.glyphs[0].pixel_y1);
}
