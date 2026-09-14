const bytes = @import("../../../../bytes.zig");
const zigstd = @import("../../../../zigstd.zig");
const LocaTable = @import("../loca.zig").LocaTable;

const Allocator = zigstd.Allocator;

const SimpleGlyph = @import("SimpleGlyph.zig");
const CompositeGlyph = @import("CompositeGlyph.zig");

// Glyph outlines. Unlike the fixed-size tables, glyf isn't decoded up front:
// a font can have thousands of glyphs and callers only ever need a handful
// per frame, so GlyfTable just borrows the raw table bytes and getGlyph()
// decodes one glyph's outline on demand, given the byte range LocaTable
// hands back for that glyph index.
//
// getGlyph() decodes the stream encoding fully but does not walk composite
// trees or flatten anything into renderable contours - that belongs to a
// separate consumer (shape.zig) that recurses through GlyfTable+LocaTable
// together.
const Self = @This();

data: []const u8,

pub fn init(table_data: []const u8) Self {
    return .{ .data = table_data };
}

pub const GlyphHeader = extern struct {
    /// >= 0: a simple glyph with this many contours. < 0: a composite
    /// glyph (-1 is the spec-recommended value, but any negative works).
    number_of_contours: i16,
    x_min: i16,
    y_min: i16,
    x_max: i16,
    y_max: i16,
};

pub const Glyph = union(enum) {
    simple: SimpleGlyph,
    composite: CompositeGlyph,

    pub fn destroy(self: *Glyph, alloc: Allocator) void {
        switch (self.*) {
            .simple => |*s| s.destroy(alloc),
            .composite => |*c| c.destroy(alloc),
        }
    }
};

pub const Error = error{
    /// argument1/argument2 identify points to align rather than an (x,y)
    /// offset; this requires walking the parent's already-placed points
    /// and isn't supported yet (it's also rare - the spec requires
    /// ARGS_ARE_XY_VALUES on every composite glyph's first component).
    PointMatchingNotSupported,
};

/// Decodes the glyph occupying `range` within this table's bytes. Range
/// comes from LocaTable.getGlyphRange(glyph); a null range (no outline)
/// isn't representable here since there's nothing to decode.
pub fn getGlyph(self: Self, alloc: Allocator, range: LocaTable.GlyphRange, check: bool) !Glyph {
    if (check) try bytes.checkRange(self.data, range.start, range.end);
    const glyph_data = self.data[range.start..range.end];

    if (check) try bytes.checkSize(glyph_data, 0, GlyphHeader);
    const header = bytes.typeBEToNative(GlyphHeader, glyph_data);

    if (header.number_of_contours >= 0) {
        return .{ .simple = try SimpleGlyph.parse(alloc, glyph_data, header, check) };
    }
    return .{ .composite = try CompositeGlyph.parse(alloc, glyph_data, header, check) };
}

test "GlyfTable.GlyphHeader matches the stream layout (10 bytes, no padding)" {
    try zigstd.testExpectEqual(@as(usize, 10), @sizeOf(GlyphHeader));
    try zigstd.testExpectEqual(@as(usize, 2), @offsetOf(GlyphHeader, "x_min"));
    try zigstd.testExpectEqual(@as(usize, 8), @offsetOf(GlyphHeader, "y_max"));
}

// One contour, 4 points, exercising all four x/y coordinate encodings:
// point0 short+positive, point1 short+negative, point2 long (i16), point3
// "same as previous" (zero delta, no bytes consumed).
fn buildSimpleGlyphCodingVariants(alloc: Allocator) ![]u8 {
    const total = 10 + 2 + 2 + 0 + 4 + 4 + 4;
    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);

    bytes.setU16BE(buf, @bitCast(@as(i16, 1))); // number_of_contours
    bytes.setU16BE(buf[10..], 3); // endPtsOfContours[0] = 3 (4 points)
    bytes.setU16BE(buf[12..], 0); // instructionLength

    const flags_start = 14;
    buf[flags_start + 0] = 0x37; // on_curve, x_short+pos, y_short+pos
    buf[flags_start + 1] = 0x07; // on_curve, x_short+neg, y_short+neg
    buf[flags_start + 2] = 0x01; // on_curve, x long, y long
    buf[flags_start + 3] = 0x31; // on_curve, x same, y same

    const x_start = flags_start + 4;
    buf[x_start] = 10; // point0: +10
    buf[x_start + 1] = 5; // point1: -5
    bytes.setI16BE(buf[x_start + 2 ..], 1000); // point2: +1000 (2 bytes)
    // point3: no bytes

    const y_start = x_start + 4;
    buf[y_start] = 20; // point0: +20
    buf[y_start + 1] = 8; // point1: -8
    bytes.setI16BE(buf[y_start + 2 ..], 2000); // point2: +2000 (2 bytes)
    // point3: no bytes

    return buf;
}

test "GlyfTable.getGlyph (simple) accumulates absolute coordinates across every encoding variant" {
    const alloc = zigstd.testAllocator;
    const glyph_buf = try buildSimpleGlyphCodingVariants(alloc);
    defer alloc.free(glyph_buf);

    const table = Self.init(glyph_buf);

    var glyph = try table.getGlyph(alloc, .{ .start = 0, .end = @intCast(glyph_buf.len) }, true);
    defer glyph.destroy(alloc);

    const simple = glyph.simple;
    try zigstd.testExpectEqualSlices(u16, &[_]u16{3}, simple.end_pts_of_contours);
    try zigstd.testExpectEqualSlices(i16, &[_]i16{ 10, 5, 1005, 1005 }, simple.x_coordinates);
    try zigstd.testExpectEqualSlices(i16, &[_]i16{ 20, 12, 2012, 2012 }, simple.y_coordinates);
}

test "GlyfTable.getGlyph (simple) expands a repeated flag run to one entry per point" {
    const alloc = zigstd.testAllocator;
    // 5 points sharing one flag (on_curve, x_short+pos, y_short+pos),
    // written once plus a repeat count of 4.
    const total = 10 + 2 + 2 + 2 + 5 + 5;
    const buf = try alloc.alloc(u8, total);
    defer alloc.free(buf);
    @memset(buf, 0);

    bytes.setU16BE(buf, @bitCast(@as(i16, 1))); // number_of_contours
    bytes.setU16BE(buf[10..], 4); // endPtsOfContours[0] = 4 (5 points)
    bytes.setU16BE(buf[12..], 0); // instructionLength

    buf[14] = 0x37 | 0x08; // base flag with REPEAT_FLAG set
    buf[15] = 4; // repeat 4 additional times (5 total)

    const x_start = 16;
    for (0..5) |i| buf[x_start + i] = 10; // each point: +10
    const y_start = x_start + 5;
    for (0..5) |i| buf[y_start + i] = 5; // each point: +5

    const table = Self.init(buf);

    var glyph = try table.getGlyph(alloc, .{ .start = 0, .end = @intCast(buf.len) }, true);
    defer glyph.destroy(alloc);

    const simple = glyph.simple;
    try zigstd.testExpectEqual(@as(usize, 5), simple.flags.len);
    try zigstd.testExpectEqualSlices(i16, &[_]i16{ 10, 20, 30, 40, 50 }, simple.x_coordinates);
    try zigstd.testExpectEqualSlices(i16, &[_]i16{ 5, 10, 15, 20, 25 }, simple.y_coordinates);
}

test "GlyfTable.getGlyph (composite) decodes a single scaled component" {
    const alloc = zigstd.testAllocator;
    const total = 10 + 2 + 2 + 2 + 2 + 2;
    const buf = try alloc.alloc(u8, total);
    defer alloc.free(buf);
    @memset(buf, 0);

    bytes.setI16BE(buf, -1); // number_of_contours: composite

    // flags: ARG_1_AND_2_ARE_WORDS | ARGS_ARE_XY_VALUES | WE_HAVE_A_SCALE
    bytes.setU16BE(buf[10..], 0x0001 | 0x0002 | 0x0008);
    bytes.setU16BE(buf[12..], 5); // glyphIndex
    bytes.setI16BE(buf[14..], 100); // dx
    bytes.setI16BE(buf[16..], -50); // dy
    bytes.setI16BE(buf[18..], 0x2000); // scale = 0.5 (F2Dot14)

    const table = Self.init(buf);

    var glyph = try table.getGlyph(alloc, .{ .start = 0, .end = @intCast(buf.len) }, true);
    defer glyph.destroy(alloc);

    const composite = glyph.composite;
    try zigstd.testExpectEqual(@as(usize, 1), composite.components.len);
    const c = composite.components[0];
    try zigstd.testExpectEqual(@as(u16, 5), @backingInt(c.glyph_index));
    try zigstd.testExpectEqual(@as(f32, 100), c.dx);
    try zigstd.testExpectEqual(@as(f32, -50), c.dy);
    try zigstd.testExpectEqual(@as(f32, 0.5), c.xscale);
    try zigstd.testExpectEqual(@as(f32, 0.5), c.yscale);
    try zigstd.testExpectEqual(@as(f32, 0), c.scale01);
    try zigstd.testExpectEqual(@as(f32, 0), c.scale10);
    try zigstd.testExpectEqual(@as(usize, 0), composite.instructions.len);
}

test "GlyfTable.getGlyph (composite) chains MORE_COMPONENTS and reads trailing instructions" {
    const alloc = zigstd.testAllocator;
    // component 0: byte args, WE_HAVE_AN_X_AND_Y_SCALE, MORE_COMPONENTS set
    // component 1: word args, WE_HAVE_A_TWO_BY_TWO, WE_HAVE_INSTRUCTIONS
    const total = 10 + 10 + 16 + 2 + 2;
    const buf = try alloc.alloc(u8, total);
    defer alloc.free(buf);
    @memset(buf, 0);

    bytes.setI16BE(buf, -1); // number_of_contours: composite

    var off: usize = 10;
    // component 0
    bytes.setU16BE(buf[off..], 0x0002 | 0x0020 | 0x0040); // ARGS_ARE_XY_VALUES | MORE_COMPONENTS | X_AND_Y_SCALE
    bytes.setU16BE(buf[off + 2 ..], 7); // glyphIndex
    buf[off + 4] = 10; // dx (i8)
    buf[off + 5] = @bitCast(@as(i8, -20)); // dy (i8)
    bytes.setI16BE(buf[off + 6 ..], 0x4000); // xscale = 1.0
    bytes.setI16BE(buf[off + 8 ..], 0x6000); // yscale = 1.5
    off += 10;

    // component 1
    bytes.setU16BE(buf[off..], 0x0001 | 0x0002 | 0x0080 | 0x0100); // WORDS | XY_VALUES | TWO_BY_TWO | INSTRUCTIONS
    bytes.setU16BE(buf[off + 2 ..], 9); // glyphIndex
    bytes.setI16BE(buf[off + 4 ..], 0); // dx
    bytes.setI16BE(buf[off + 6 ..], 0); // dy
    bytes.setI16BE(buf[off + 8 ..], 0x4000); // xscale = 1.0
    bytes.setI16BE(buf[off + 10 ..], 0x1000); // scale01 = 0.25
    bytes.setI16BE(buf[off + 12 ..], @bitCast(@as(i16, -4096))); // scale10 = -0.25
    bytes.setI16BE(buf[off + 14 ..], 0x4000); // yscale = 1.0
    off += 16;

    bytes.setU16BE(buf[off..], 2); // numInstr
    buf[off + 2] = 0xAA;
    buf[off + 3] = 0xBB;

    const table = Self.init(buf);

    var glyph = try table.getGlyph(alloc, .{ .start = 0, .end = @intCast(buf.len) }, true);
    defer glyph.destroy(alloc);

    const composite = glyph.composite;
    try zigstd.testExpectEqual(@as(usize, 2), composite.components.len);

    const c0 = composite.components[0];
    try zigstd.testExpectEqual(@as(u16, 7), @backingInt(c0.glyph_index));
    try zigstd.testExpectEqual(@as(f32, 10), c0.dx);
    try zigstd.testExpectEqual(@as(f32, -20), c0.dy);
    try zigstd.testExpectEqual(@as(f32, 1.0), c0.xscale);
    try zigstd.testExpectEqual(@as(f32, 1.5), c0.yscale);

    const c1 = composite.components[1];
    try zigstd.testExpectEqual(@as(u16, 9), @backingInt(c1.glyph_index));
    try zigstd.testExpectEqual(@as(f32, 1.0), c1.xscale);
    try zigstd.testExpectEqual(@as(f32, 0.25), c1.scale01);
    try zigstd.testExpectEqual(@as(f32, -0.25), c1.scale10);
    try zigstd.testExpectEqual(@as(f32, 1.0), c1.yscale);

    try zigstd.testExpectEqualSlices(u8, &[_]u8{ 0xAA, 0xBB }, composite.instructions);
}

test "GlyfTable.getGlyph (composite) rejects point-matching components" {
    const alloc = zigstd.testAllocator;
    var buf: [14]u8 = @splat(0);
    bytes.setI16BE(&buf, -1); // number_of_contours: composite
    bytes.setU16BE(buf[10..], 0x0001); // ARG_1_AND_2_ARE_WORDS only: ARGS_ARE_XY_VALUES clear
    bytes.setU16BE(buf[12..], 3); // glyphIndex

    const table = Self.init(&buf);

    try zigstd.testExpectError(
        Error.PointMatchingNotSupported,
        table.getGlyph(alloc, .{ .start = 0, .end = buf.len }, true),
    );
}
