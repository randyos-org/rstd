const bytes = @import("../../../../bytes.zig");
const zigstd = @import("../../../../zigstd.zig");

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

fn buildGpos(alloc: zigstd.Allocator) ![]u8 {
    // header(10) -> lookupList at 10 -> [numLookups(2)=1, lookupOffset(2)]
    // -> Lookup at 16 -> [type(2)=2, flag(2)=0, count(2)=1, subOff(2)]
    // -> PairPosFormat1 at 26 -> [format=1, covOff, vf1=X_ADVANCE, vf2=0, pairSetCount=1, pairSetOff]
    const header_size = 10;
    const lookup_list_start = header_size;
    const lookup_start = lookup_list_start + 4;
    const subtable_start = lookup_start + 8;
    const pair_pos_header_size = 10;
    const coverage_start = subtable_start + pair_pos_header_size + 2; // one pairSetOffset
    const pair_set_start = coverage_start + 6; // coverage fmt1: format+glyphCount+1 glyph
    const total = pair_set_start + 2 + 4; // pairValueCount + one PairValueRecord (secondGlyph+xAdvance)

    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);

    bytes.setU16BE(buf[0..], 1); // majorVersion
    bytes.setU16BE(buf[2..], 0); // minorVersion
    bytes.setU16BE(buf[4..], 0); // scriptListOffset (unused)
    bytes.setU16BE(buf[6..], 0); // featureListOffset (unused)
    bytes.setU16BE(buf[8..], lookup_list_start); // lookupListOffset

    bytes.setU16BE(buf[lookup_list_start..], 1); // numLookups
    bytes.setU16BE(buf[lookup_list_start + 2 ..], @intCast(lookup_start - lookup_list_start)); // lookupOffsets[0]

    bytes.setU16BE(buf[lookup_start..], 2); // lookupType (Pair Adjustment)
    bytes.setU16BE(buf[lookup_start + 2 ..], 0); // lookupFlag
    bytes.setU16BE(buf[lookup_start + 4 ..], 1); // subtableCount
    bytes.setU16BE(buf[lookup_start + 6 ..], @intCast(subtable_start - lookup_start)); // subtableOffsets[0]

    bytes.setU16BE(buf[subtable_start..], 1); // format
    bytes.setU16BE(buf[subtable_start + 2 ..], @intCast(coverage_start - subtable_start)); // coverageOffset
    const vf1: ValueFormat = .{ .x_placement = false, .y_placement = false, .x_advance = true, .y_advance = false, .x_placement_device = false, .y_placement_device = false, .x_advance_device = false, .y_advance_device = false };
    bytes.setU16BE(buf[subtable_start + 4 ..], @bitCast(vf1)); // valueFormat1
    bytes.setU16BE(buf[subtable_start + 6 ..], 0); // valueFormat2
    bytes.setU16BE(buf[subtable_start + 8 ..], 1); // pairSetCount
    bytes.setU16BE(buf[subtable_start + 10 ..], @intCast(pair_set_start - subtable_start)); // pairSetOffsets[0]

    bytes.setU16BE(buf[coverage_start..], 1); // Coverage format 1
    bytes.setU16BE(buf[coverage_start + 2 ..], 1); // glyphCount
    bytes.setU16BE(buf[coverage_start + 4 ..], 5); // glyphArray[0] = glyph 5

    bytes.setU16BE(buf[pair_set_start..], 1); // pairValueCount
    bytes.setU16BE(buf[pair_set_start + 2 ..], 10); // secondGlyph = 10
    bytes.setI16BE(buf[pair_set_start + 4 ..], -25); // xAdvance

    return buf;
}

test "GPOS structural accessors walk header -> LookupList -> Lookup -> PairPosFormat1 -> Coverage/PairSet" {
    const alloc = zigstd.testAllocator;
    const buf = try buildGpos(alloc);
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
