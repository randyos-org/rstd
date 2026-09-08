const types = @import("../types.zig");
const bytes = @import("../../../bytes.zig");
const zigstd = @import("../../../zigstd.zig");

const Allocator = zigstd.Allocator;

pub const HeadTable = extern struct {
    const Self = @This();

    version: types.Fixed align(2),
    font_revision: types.Fixed align(2),
    checksum_adjustment: u32 align(2),
    magic_number: u32 align(2),
    flags: u16,
    units_per_em: u16,
    created: i64 align(2),
    modified: i64 align(2),
    x_min: i16,
    y_min: i16,
    x_max: i16,
    y_max: i16,
    mac_style: u16,
    lowest_rec_ppem: u16,
    font_direction_hint: i16,
    index_to_loc_format: i16,
    glyph_data_format: i16,

    pub fn calcTableChecksum(self: *Self, table_data: *[]const u8) u32 {
        var sum: u32 = self.version + self.font_revision;
        // assume checksum_adjustment is 0 for this calculation
        const nLongs = (@sizeOf(Self) + 3) / 4 - 3; // subtract 3 for the two Fixed fields and checksum_adjustment
        for (0..nLongs) |i| sum += @as(u32, table_data[12 + i * 4]);
        return sum;
    }

    /// When calculating the raw_font_checksum_adj0, the checksum_adjustment field is assumed to be 0.
    /// You may need to actually temporarily set it to zero since we don't necessarily know the packing scheme,
    /// so set it to zero, calculate the checksum for the entire font byte array with the zero checksum,
    /// and pass that value here.  Do not forget to then re-set the checksum_adjustment following the calculation.
    pub fn calcChecksumAdj(raw_font_checksum_adj0: u32) u32 {
        return 0xB1B0AFBA - raw_font_checksum_adj0;
    }

    pub fn create(alloc: Allocator, table_data: []const u8, check: bool) !*Self {
        if (check) try bytes.checkSize(table_data, 0, Self);
        return bytes.allocTypeBEToNative(alloc, Self, table_data);
    }

    pub fn destroy(self: *Self, alloc: Allocator) void {
        alloc.destroy(self);
    }
};

test "HeadTable matches the OpenType wire layout (54 bytes, no padding)" {
    try zigstd.testExpectEqual(@as(usize, 54), @sizeOf(HeadTable));
    try zigstd.testExpectEqual(@as(usize, 0), @offsetOf(HeadTable, "version"));
    try zigstd.testExpectEqual(@as(usize, 4), @offsetOf(HeadTable, "font_revision"));
    try zigstd.testExpectEqual(@as(usize, 8), @offsetOf(HeadTable, "checksum_adjustment"));
    try zigstd.testExpectEqual(@as(usize, 12), @offsetOf(HeadTable, "magic_number"));
    try zigstd.testExpectEqual(@as(usize, 16), @offsetOf(HeadTable, "flags"));
    try zigstd.testExpectEqual(@as(usize, 18), @offsetOf(HeadTable, "units_per_em"));
    try zigstd.testExpectEqual(@as(usize, 20), @offsetOf(HeadTable, "created"));
    try zigstd.testExpectEqual(@as(usize, 28), @offsetOf(HeadTable, "modified"));
    try zigstd.testExpectEqual(@as(usize, 36), @offsetOf(HeadTable, "x_min"));
    try zigstd.testExpectEqual(@as(usize, 50), @offsetOf(HeadTable, "index_to_loc_format"));
    try zigstd.testExpectEqual(@as(usize, 52), @offsetOf(HeadTable, "glyph_data_format"));
}

test "HeadTable.create reads a real big-endian 'head' table correctly" {
    const alloc = zigstd.testAllocator;
    var buf: [54]u8 = @splat(0);
    bytes.writeIntBE(u32, buf[12..16], 0x5F0F3CF5); // magic_number
    bytes.writeIntBE(u16, buf[16..18], 0x0003); // flags
    bytes.writeIntBE(u16, buf[18..20], 2048); // units_per_em
    bytes.writeIntBE(i16, buf[36..38], -100); // x_min

    const head = try HeadTable.create(alloc, buf[0..], true);
    defer head.destroy(alloc);

    try zigstd.testExpectEqual(@as(u32, 0x5F0F3CF5), head.magic_number);
    try zigstd.testExpectEqual(@as(u16, 3), head.flags);
    try zigstd.testExpectEqual(@as(u16, 2048), head.units_per_em);
    try zigstd.testExpectEqual(@as(i16, -100), head.x_min);
}
