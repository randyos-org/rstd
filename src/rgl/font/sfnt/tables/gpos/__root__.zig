const bytes = @import("../../../../../bytes.zig");
const zigstd = @import("../../../../../zigstd.zig");
const types = @import("../../types.zig");

const gpos_table = @import("GposTable.zig");
pub const GposTable = gpos_table;
pub const Lookup = gpos_table.Lookup;

const pair_pos = @import("PairPos.zig");
pub const ValueFormat = pair_pos.ValueFormat;
pub const PairPosFormat1 = pair_pos.PairPosFormat1;
pub const PairPosFormat2 = pair_pos.PairPosFormat2;
pub const PairSet = pair_pos.PairSet;

pub const Coverage = @import("Coverage.zig");
pub const ClassDef = @import("ClassDef.zig");

/// GPOS header, one LookupList entry and one type-2 Lookup whose single
/// subtable starts at `subtable_start`.
fn writeGposLookup(buf: []u8, minor_version: u16, subtable_start: usize) void {
    const lookup_list_start = 10;
    const lookup_start = lookup_list_start + 4;

    bytes.setU16BE(buf[0..], 1); // majorVersion
    bytes.setU16BE(buf[2..], minor_version);
    bytes.setU16BE(buf[4..], 0); // scriptListOffset (unused)
    bytes.setU16BE(buf[6..], 0); // featureListOffset (unused)
    bytes.setU16BE(buf[8..], lookup_list_start); // lookupListOffset

    bytes.setU16BE(buf[lookup_list_start..], 1); // numLookups
    bytes.setU16BE(buf[lookup_list_start + 2 ..], lookup_start - lookup_list_start); // lookupOffsets[0]

    bytes.setU16BE(buf[lookup_start..], 2); // lookupType (Pair Adjustment)
    bytes.setU16BE(buf[lookup_start + 2 ..], 0); // lookupFlag
    bytes.setU16BE(buf[lookup_start + 4 ..], 1); // subtableCount
    bytes.setU16BE(buf[lookup_start + 6 ..], @intCast(subtable_start - lookup_start)); // subtableOffsets[0]
}

/// PairPosFormat1 kerning glyph 5 against glyph 10. `values` fills the one
/// PairValueRecord's valueRecord1 then valueRecord2, field by field.
fn buildGpos(alloc: zigstd.Allocator, vf1: u16, vf2: u16, minor_version: u16, values: []const i16) ![]u8 {
    const subtable_start = 10 + 4 + 8;
    const pair_pos_header_size = 10;
    const coverage_start = subtable_start + pair_pos_header_size + 2; // one pairSetOffset
    const pair_set_start = coverage_start + 6; // coverage fmt1: format+glyphCount+1 glyph
    const total = pair_set_start + 2 + 2 + 2 * values.len; // pairValueCount + secondGlyph + values

    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);
    writeGposLookup(buf, minor_version, subtable_start);

    bytes.setU16BE(buf[subtable_start..], 1); // format
    bytes.setU16BE(buf[subtable_start + 2 ..], coverage_start - subtable_start); // coverageOffset
    bytes.setU16BE(buf[subtable_start + 4 ..], vf1); // valueFormat1
    bytes.setU16BE(buf[subtable_start + 6 ..], vf2); // valueFormat2
    bytes.setU16BE(buf[subtable_start + 8 ..], 1); // pairSetCount
    bytes.setU16BE(buf[subtable_start + 10 ..], pair_set_start - subtable_start); // pairSetOffsets[0]

    bytes.setU16BE(buf[coverage_start..], 1); // Coverage format 1
    bytes.setU16BE(buf[coverage_start + 2 ..], 1); // glyphCount
    bytes.setU16BE(buf[coverage_start + 4 ..], 5); // glyphArray[0] = glyph 5

    bytes.setU16BE(buf[pair_set_start..], 1); // pairValueCount
    bytes.setU16BE(buf[pair_set_start + 2 ..], 10); // secondGlyph = 10
    for (values, 0..) |v, i| bytes.setI16BE(buf[pair_set_start + 4 + 2 * i ..], v);

    return buf;
}

/// PairPosFormat2 with 2x2 classes: glyph 5 is class 1 of ClassDef1,
/// glyph 10 class 1 of ClassDef2. Class2Records hold vf1 = xPlacement |
/// xAdvance and vf2 = xAdvance; (1, 1)'s first xAdvance is -40.
fn buildGposClass(alloc: zigstd.Allocator) ![]u8 {
    const subtable_start = 10 + 4 + 8;
    const record_size = 6;
    const coverage_start = subtable_start + 16 + 2 * 2 * record_size;
    const class_def1_start = coverage_start + 6;
    const class_def2_start = class_def1_start + 8;
    const total = class_def2_start + 8;

    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);
    writeGposLookup(buf, 0, subtable_start);

    bytes.setU16BE(buf[subtable_start..], 2); // format
    bytes.setU16BE(buf[subtable_start + 2 ..], coverage_start - subtable_start); // coverageOffset
    bytes.setU16BE(buf[subtable_start + 4 ..], 0x0005); // valueFormat1: xPlacement | xAdvance
    bytes.setU16BE(buf[subtable_start + 6 ..], 0x0004); // valueFormat2: xAdvance
    bytes.setU16BE(buf[subtable_start + 8 ..], class_def1_start - subtable_start);
    bytes.setU16BE(buf[subtable_start + 10 ..], class_def2_start - subtable_start);
    bytes.setU16BE(buf[subtable_start + 12 ..], 2); // class1Count
    bytes.setU16BE(buf[subtable_start + 14 ..], 2); // class2Count

    const records = subtable_start + 16;
    for (0..2) |c1| for (0..2) |c2| {
        const rec = records + record_size * (c1 * 2 + c2);
        bytes.setI16BE(buf[rec..], 7); // valueRecord1.xPlacement
        bytes.setI16BE(buf[rec + 2 ..], if (c1 == 1 and c2 == 1) -40 else -1); // valueRecord1.xAdvance
        bytes.setI16BE(buf[rec + 4 ..], 3); // valueRecord2.xAdvance
    };

    bytes.setU16BE(buf[coverage_start..], 1); // Coverage format 1
    bytes.setU16BE(buf[coverage_start + 2 ..], 1); // glyphCount
    bytes.setU16BE(buf[coverage_start + 4 ..], 5); // glyphArray[0] = glyph 5

    for ([_]usize{ class_def1_start, class_def2_start }, [_]u16{ 5, 10 }) |start, first_glyph| {
        bytes.setU16BE(buf[start..], 1); // ClassDef format 1
        bytes.setU16BE(buf[start + 2 ..], first_glyph); // startGlyphID
        bytes.setU16BE(buf[start + 4 ..], 1); // glyphCount
        bytes.setU16BE(buf[start + 6 ..], 1); // classValueArray[0]
    }

    return buf;
}

fn glyph(i: u16) types.GlyphIndex {
    return @fromBackingInt(i);
}

test "GPOS structural accessors walk header -> LookupList -> Lookup -> PairPosFormat1 -> Coverage/PairSet" {
    const alloc = zigstd.testAllocator;
    const buf = try buildGpos(alloc, 0x0004, 0, 0, &.{-25});
    defer alloc.free(buf);

    const gpos = try GposTable.init(buf, true);
    try zigstd.testExpectEqual(@as(u16, 1), gpos.majorVersion());
    try zigstd.testExpectEqual(@as(u16, 1), gpos.numLookups());

    const lookup = gpos.getLookup(0);
    try zigstd.testExpectEqual(@as(u16, 2), lookup.lookupType());
    try zigstd.testExpectEqual(@as(u16, 1), lookup.subtableCount());

    const pp: PairPosFormat1 = .{ .table_data = lookup.getSubtableData(0) };
    try zigstd.testExpectEqual(@as(u16, 1), pp.format());
    try zigstd.testExpectEqual(@as(u16, 1), pp.pairSetCount());
    try zigstd.testExpectTrue(pp.valueFormat1().x_advance);
    try zigstd.testExpectEqual(@as(usize, 2), pp.valueFormat1().size());

    const coverage = pp.getCoverage();
    try zigstd.testExpectEqual(@as(u16, 1), coverage.format());
    try zigstd.testExpectEqual(@as(u16, 1), coverage.glyphCount());
    try zigstd.testExpectEqual(@as(u16, 5), @backingInt(coverage.getGlyph(0)));

    const pair_set = pp.getPairSet(0);
    try zigstd.testExpectEqual(@as(u16, 1), pair_set.pairValueCount());
    try zigstd.testExpectEqual(@as(u16, 10), @backingInt(pair_set.getSecondGlyph(0)));

    const xa_off = pp.valueFormat1().xAdvanceOffset().?;
    const rec_off = pair_set.getPairValueRecordOffset(0);
    const x_advance = bytes.readIntBE(i16, pair_set.table_data[rec_off + @sizeOf(u16) + xa_off ..]);
    try zigstd.testExpectEqual(@as(i16, -25), x_advance);
}

test "glyphKernAdvance finds xAdvance past xPlacement in PairPosFormat1, GPOS 1.1" {
    const alloc = zigstd.testAllocator;
    // vf1 = xPlacement | xAdvance, vf2 = xPlacement
    const buf = try buildGpos(alloc, 0x0005, 0x0001, 1, &.{ 7, -25, 3 });
    defer alloc.free(buf);

    const gpos = try GposTable.init(buf, true);
    try zigstd.testExpectEqual(@as(i16, -25), gpos.glyphKernAdvance(glyph(5), glyph(10)));
    try zigstd.testExpectEqual(@as(i16, 0), gpos.glyphKernAdvance(glyph(5), glyph(11)));
    try zigstd.testExpectEqual(@as(i16, 0), gpos.glyphKernAdvance(glyph(6), glyph(10)));
}

test "glyphKernAdvance skips a pair adjustment with no first-glyph xAdvance" {
    const alloc = zigstd.testAllocator;
    const buf = try buildGpos(alloc, 0x0001, 0x0004, 0, &.{ 7, -25 });
    defer alloc.free(buf);

    const gpos = try GposTable.init(buf, true);
    try zigstd.testExpectEqual(@as(i16, 0), gpos.glyphKernAdvance(glyph(5), glyph(10)));
}

test "glyphKernAdvance indexes PairPosFormat2 class records by both value formats" {
    const alloc = zigstd.testAllocator;
    const buf = try buildGposClass(alloc);
    defer alloc.free(buf);

    const gpos = try GposTable.init(buf, true);
    try zigstd.testExpectEqual(@as(i16, -40), gpos.glyphKernAdvance(glyph(5), glyph(10)));
    // glyph 11 is unclassed (class 0)
    try zigstd.testExpectEqual(@as(i16, -1), gpos.glyphKernAdvance(glyph(5), glyph(11)));
}
