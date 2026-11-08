const types = @import("../types.zig");
const bytes = @import("../../../../bytes.zig");
const zigstd = @import("../../../../zigstd.zig");

/// The 'head' table.
const StreamLayout = extern struct {
    version: types.Fixed align(2),
    font_revision: types.Fixed align(2),
    checksum_adjustment: u32 align(2),
    magic_number: u32 align(2),
    flags: u16,
    units_per_em: u16,
    created: types.LongDateTime align(2),
    modified: types.LongDateTime align(2),
    x_min: i16,
    y_min: i16,
    x_max: i16,
    y_max: i16,
    mac_style: u16,
    lowest_rec_ppem: u16,
    font_direction_hint: i16,
    index_to_loc_format: i16,
    glyph_data_format: i16,
};

/// Font header.
pub const HeadTable = struct {
    const Self = @This();

    table_data: []const u8,

    pub fn init(table_data: []const u8, check: bool) !Self {
        if (check) try bytes.checkSize(table_data, 0, StreamLayout);
        return .{ .table_data = table_data };
    }

    pub inline fn version(self: Self) types.Fixed {
        return types.Fixed.fromBytes(self.table_data[@offsetOf(StreamLayout, "version")..]);
    }

    pub inline fn fontRevision(self: Self) types.Fixed {
        return types.Fixed.fromBytes(self.table_data[@offsetOf(StreamLayout, "font_revision")..]);
    }

    pub inline fn checksumAdjustment(self: Self) u32 {
        return bytes.getU32BE(self.table_data[@offsetOf(StreamLayout, "checksum_adjustment")..]);
    }

    pub inline fn magicNumber(self: Self) u32 {
        return bytes.getU32BE(self.table_data[@offsetOf(StreamLayout, "magic_number")..]);
    }

    pub inline fn flags(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "flags")..]);
    }

    pub inline fn unitsPerEm(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "units_per_em")..]);
    }

    pub inline fn created(self: Self) types.LongDateTime {
        return bytes.readIntBE(types.LongDateTime, self.table_data[@offsetOf(StreamLayout, "created")..]);
    }

    pub inline fn modified(self: Self) types.LongDateTime {
        return bytes.readIntBE(types.LongDateTime, self.table_data[@offsetOf(StreamLayout, "modified")..]);
    }

    pub inline fn xMin(self: Self) i16 {
        return bytes.getI16BE(self.table_data[@offsetOf(StreamLayout, "x_min")..]);
    }

    pub inline fn yMin(self: Self) i16 {
        return bytes.getI16BE(self.table_data[@offsetOf(StreamLayout, "y_min")..]);
    }

    pub inline fn xMax(self: Self) i16 {
        return bytes.getI16BE(self.table_data[@offsetOf(StreamLayout, "x_max")..]);
    }

    pub inline fn yMax(self: Self) i16 {
        return bytes.getI16BE(self.table_data[@offsetOf(StreamLayout, "y_max")..]);
    }

    pub inline fn macStyle(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "mac_style")..]);
    }

    pub inline fn lowestRecPpem(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "lowest_rec_ppem")..]);
    }

    pub inline fn fontDirectionHint(self: Self) i16 {
        return bytes.getI16BE(self.table_data[@offsetOf(StreamLayout, "font_direction_hint")..]);
    }

    pub inline fn indexToLocFormat(self: Self) i16 {
        return bytes.getI16BE(self.table_data[@offsetOf(StreamLayout, "index_to_loc_format")..]);
    }

    pub inline fn glyphDataFormat(self: Self) i16 {
        return bytes.getI16BE(self.table_data[@offsetOf(StreamLayout, "glyph_data_format")..]);
    }

    /// When calculating the raw_font_checksum_adj0, the checksum_adjustment field is assumed to be 0.
    /// You may need to actually temporarily set it to zero since we don't necessarily know the packing scheme,
    /// so set it to zero, calculate the checksum for the entire font byte array with the zero checksum,
    /// and pass that value here.  Do not forget to then re-set the checksum_adjustment following the calculation.
    pub fn calcChecksumAdj(raw_font_checksum_adj0: u32) u32 {
        return 0xB1B0AFBA - raw_font_checksum_adj0;
    }
};

test "head StreamLayout matches the OpenType stream layout (54 bytes, no padding)" {
    try zigstd.testExpectEqual(@as(usize, 54), @sizeOf(StreamLayout));
    try zigstd.testExpectEqual(@as(usize, 0), @offsetOf(StreamLayout, "version"));
    try zigstd.testExpectEqual(@as(usize, 4), @offsetOf(StreamLayout, "font_revision"));
    try zigstd.testExpectEqual(@as(usize, 8), @offsetOf(StreamLayout, "checksum_adjustment"));
    try zigstd.testExpectEqual(@as(usize, 12), @offsetOf(StreamLayout, "magic_number"));
    try zigstd.testExpectEqual(@as(usize, 16), @offsetOf(StreamLayout, "flags"));
    try zigstd.testExpectEqual(@as(usize, 18), @offsetOf(StreamLayout, "units_per_em"));
    try zigstd.testExpectEqual(@as(usize, 20), @offsetOf(StreamLayout, "created"));
    try zigstd.testExpectEqual(@as(usize, 28), @offsetOf(StreamLayout, "modified"));
    try zigstd.testExpectEqual(@as(usize, 36), @offsetOf(StreamLayout, "x_min"));
    try zigstd.testExpectEqual(@as(usize, 50), @offsetOf(StreamLayout, "index_to_loc_format"));
    try zigstd.testExpectEqual(@as(usize, 52), @offsetOf(StreamLayout, "glyph_data_format"));
}

test "HeadTable getters read a real big-endian 'head' table correctly" {
    var buf: [54]u8 = @splat(0);
    bytes.setU32BE(buf[12..], 0x5F0F3CF5); // magic_number
    bytes.setU16BE(buf[16..], 0x0003); // flags
    bytes.setU16BE(buf[18..], 2048); // units_per_em
    bytes.setI16BE(buf[36..], -100); // x_min

    const head = try HeadTable.init(buf[0..], true);

    try zigstd.testExpectEqual(@as(u32, 0x5F0F3CF5), head.magicNumber());
    try zigstd.testExpectEqual(@as(u16, 3), head.flags());
    try zigstd.testExpectEqual(@as(u16, 2048), head.unitsPerEm());
    try zigstd.testExpectEqual(@as(i16, -100), head.xMin());
}
