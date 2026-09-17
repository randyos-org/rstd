const std = @import("std");
const zigstd = @import("../../../../zigstd.zig");
const Allocator = zigstd.Allocator;

const Sfnt = @import("../Sfnt.zig");
const types = @import("../stream_types.zig");
const GlyphIndex = types.GlyphIndex;

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
    /// Reset to null on every line break - kerning never crosses one.
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
            // The word itself is wider than the wrap width - force a
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

// ---
// Tests build a minimal synthetic SFNT byte buffer directly (same low-level
// approach Sfnt.zig's own tests use) rather than embedding a real font:
// glyph 0 (.notdef) is empty and unused, glyph 1 is a visible 200x300-funit
// box mapped from every uppercase ASCII letter, and glyph 2 is an invisible
// (no outline) 150-funit-wide space. units_per_em = 1000.

const bytes = @import("../../../../bytes.zig");
const general = @import("../tables/general.zig");
const TableTag = general.TableTag;

fn appendTable(buf: []u8, table_dir_off: usize, table_idx: usize, tag: TableTag, data: []const u8, off: usize) void {
    const rec_off = table_dir_off + table_idx * 16;
    bytes.setU32BE(buf[rec_off..], tag.toInt());
    bytes.setU32BE(buf[rec_off + 4 ..], 0); // checksum, unused
    bytes.setU32BE(buf[rec_off + 8 ..], @intCast(off));
    bytes.setU32BE(buf[rec_off + 12 ..], @intCast(data.len));
    @memcpy(buf[off..][0..data.len], data);
}

fn buildTestFont(alloc: Allocator) ![]u8 {
    const num_tables = 7;
    const table_dir_off = 12;
    const tables_start = table_dir_off + num_tables * 16;

    var head_buf: [54]u8 = @splat(0);
    bytes.setU16BE(head_buf[18..], 1000); // units_per_em
    bytes.setI16BE(head_buf[50..], 1); // index_to_loc_format: long

    var hhea_buf: [36]u8 = @splat(0);
    bytes.setI16BE(hhea_buf[4..], 800); // ascent
    bytes.setI16BE(hhea_buf[6..], -200); // descent
    bytes.setI16BE(hhea_buf[8..], 0); // line_gap
    bytes.setU16BE(hhea_buf[34..], 3); // num_of_long_hor_metrics

    var maxp_buf: [6]u8 = @splat(0);
    bytes.setU32BE(maxp_buf[0..], 0x00005000); // version 0.5
    bytes.setU16BE(maxp_buf[4..], 3); // num_glyphs

    // hmtx: glyph0 (unused/.notdef), glyph1 (box, advance 250), glyph2 (space, advance 150)
    var hmtx_buf: [12]u8 = @splat(0);
    bytes.setU16BE(hmtx_buf[4..], 250);
    bytes.setU16BE(hmtx_buf[8..], 150);

    // loca (long): glyph0 empty [0,0), glyph1 box [0,34), glyph2 empty [34,34)
    var loca_buf: [16]u8 = @splat(0);
    bytes.setU32BE(loca_buf[0..], 0);
    bytes.setU32BE(loca_buf[4..], 0);
    bytes.setU32BE(loca_buf[8..], 34);
    bytes.setU32BE(loca_buf[12..], 34);

    // glyf: one simple glyph, 1 contour, 4 on-curve points forming a
    // 200x300 rectangle, using the "long" x/y coordinate encoding (flag
    // 0x01) throughout so deltas aren't limited to a single byte.
    var glyf_buf: [34]u8 = @splat(0);
    bytes.setI16BE(glyf_buf[0..], 1); // number_of_contours
    bytes.setI16BE(glyf_buf[2..], 0); // x_min
    bytes.setI16BE(glyf_buf[4..], 0); // y_min
    bytes.setI16BE(glyf_buf[6..], 200); // x_max
    bytes.setI16BE(glyf_buf[8..], 300); // y_max
    bytes.setU16BE(glyf_buf[10..], 3); // endPtsOfContours[0] = 3 (4 points)
    bytes.setU16BE(glyf_buf[12..], 0); // instructionLength
    const flags_start = 14;
    @memset(glyf_buf[flags_start..][0..4], 0x01); // on_curve, x/y long
    const x_start = flags_start + 4;
    bytes.setI16BE(glyf_buf[x_start + 0 ..], 0); // p0 dx
    bytes.setI16BE(glyf_buf[x_start + 2 ..], 200); // p1 dx
    bytes.setI16BE(glyf_buf[x_start + 4 ..], 0); // p2 dx
    bytes.setI16BE(glyf_buf[x_start + 6 ..], -200); // p3 dx
    const y_start = x_start + 8;
    bytes.setI16BE(glyf_buf[y_start + 0 ..], 0); // p0 dy
    bytes.setI16BE(glyf_buf[y_start + 2 ..], 0); // p1 dy
    bytes.setI16BE(glyf_buf[y_start + 4 ..], 300); // p2 dy
    bytes.setI16BE(glyf_buf[y_start + 6 ..], 0); // p3 dy

    // cmap: format 0, one Windows-Unicode-BMP record. 'A'-'Z' -> glyph 1,
    // ' ' -> glyph 2, everything else unmapped (.notdef).
    var cmap_buf: [4 + 8 + 262]u8 = @splat(0);
    bytes.setU16BE(cmap_buf[0..], 0); // version
    bytes.setU16BE(cmap_buf[2..], 1); // num_tables
    bytes.setU16BE(cmap_buf[4..], 3); // platform_id: microsoft
    bytes.setU16BE(cmap_buf[6..], 1); // encoding_id: unicode_bmp
    bytes.setU32BE(cmap_buf[8..], 12); // subtable offset
    const sub = cmap_buf[12..];
    bytes.setU16BE(sub, 0); // format
    bytes.setU16BE(sub[2..], 262); // length
    bytes.setU16BE(sub[4..], 0); // language
    var c: u8 = 'A';
    while (c <= 'Z') : (c += 1) sub[6 + @as(usize, c)] = 1;
    sub[6 + @as(usize, ' ')] = 2;

    var buf = try alloc.alloc(u8, tables_start + head_buf.len + hhea_buf.len + maxp_buf.len + hmtx_buf.len + loca_buf.len + cmap_buf.len + glyf_buf.len);
    errdefer alloc.free(buf);
    @memset(buf, 0);

    bytes.setU32BE(buf[0..], 0x00010000); // scaler
    bytes.setU16BE(buf[4..], num_tables);

    var off: usize = tables_start;
    appendTable(buf, table_dir_off, 0, .head, &head_buf, off);
    off += head_buf.len;
    appendTable(buf, table_dir_off, 1, .hhea, &hhea_buf, off);
    off += hhea_buf.len;
    appendTable(buf, table_dir_off, 2, .maxp, &maxp_buf, off);
    off += maxp_buf.len;
    appendTable(buf, table_dir_off, 3, .hmtx, &hmtx_buf, off);
    off += hmtx_buf.len;
    appendTable(buf, table_dir_off, 4, .loca, &loca_buf, off);
    off += loca_buf.len;
    appendTable(buf, table_dir_off, 5, .cmap, &cmap_buf, off);
    off += cmap_buf.len;
    appendTable(buf, table_dir_off, 6, .glyf, &glyf_buf, off);
    off += glyf_buf.len;

    return buf;
}

test "layoutText lays out a single unwrapped line matching pen advance + glyph bounds" {
    const alloc = zigstd.testAllocator;
    const font_bytes = try buildTestFont(alloc);
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
    const font_bytes = try buildTestFont(alloc);
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
    const font_bytes = try buildTestFont(alloc);
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
    const font_bytes = try buildTestFont(alloc);
    defer alloc.free(font_bytes);

    const sfnt = try Sfnt.init(alloc, font_bytes, .{ .check = true });

    var layout = try layoutText(alloc, sfnt, "A\nB", 1000, .{});
    defer layout.deinit(alloc);

    try zigstd.testExpectEqual(@as(usize, 2), layout.glyphs.len);
    try zigstd.testExpectEqual(@as(i32, 0), layout.glyphs[1].pixel_x1);
    try zigstd.testExpectTrue(layout.glyphs[1].pixel_y1 < layout.glyphs[0].pixel_y1);
}
