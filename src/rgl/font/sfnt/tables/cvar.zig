const bytes = @import("../../../../bytes.zig");
const tv = @import("tuple_variation.zig");
const types = @import("../types.zig");
const zigstd = @import("../../../../zigstd.zig");

/// cvar's fixed header, majorVersion and minorVersion only; what
/// follows is the shared GlyphVariationData format (tuple_variation.zig).
const StreamLayout = extern struct {
    major_version: u16,
    minor_version: u16,
};

/// CVT variations: deltas for 'cvt ' control values, using the same Tuple
/// Variation Store format 'gvar' uses per-glyph, applied once for the whole
/// table. Not self-describing: axisCount comes from 'fvar', so init() takes
/// it (like LocaTable takes num_glyphs/long_format). Interpolating actual
/// instance values is deferred.
pub const CvarTable = struct {
    const Self = @This();

    table_data: []const u8,
    axis_count: u16,

    pub fn init(table_data: []const u8, axis_count: u16, check: bool) !Self {
        if (check) try bytes.checkSize(table_data, 0, StreamLayout);
        return .{ .table_data = table_data, .axis_count = axis_count };
    }

    pub inline fn majorVersion(self: Self) u16 {
        return types.field(StreamLayout, "major_version", self.table_data);
    }

    pub inline fn minorVersion(self: Self) u16 {
        return types.field(StreamLayout, "minor_version", self.table_data);
    }

    /// The tupleVariationCount/dataOffset/tupleVariationHeaders[] portion,
    /// in the shared GlyphVariationData format.
    pub fn getData(self: Self) tv.GlyphVariationData {
        return .{ .table_data = self.table_data[@sizeOf(StreamLayout)..], .axis_count = self.axis_count };
    }
};

test "cvar StreamLayout matches the OpenType stream layout" {
    try zigstd.testExpectEqual(@as(usize, 4), @sizeOf(StreamLayout));
    try zigstd.testExpectEqual(@as(usize, 2), @offsetOf(StreamLayout, "minor_version"));
}

test "CvarTable.getData exposes the GlyphVariationData portion after the version fields" {
    const header_size = @sizeOf(StreamLayout);
    const gvd_header = 4;
    const tvh_size = 4; // fixed_size, embedded peak tuple mandatory in cvar but omitted here for simplicity

    var buf: [header_size + gvd_header + tvh_size]u8 = @splat(0);
    bytes.setU16BE(&buf, 1); // majorVersion
    bytes.setU16BE(buf[2..], 0); // minorVersion

    const tvc: tv.TupleVariationCount = .{ .count = 1, .shared_point_numbers = false };
    bytes.setU16BE(buf[header_size..], @bitCast(tvc));
    bytes.setU16BE(buf[header_size + 2 ..], gvd_header + tvh_size); // dataOffset

    bytes.setU16BE(buf[header_size + gvd_header ..], 6); // header0 variationDataSize

    const cvar = try CvarTable.init(&buf, 1, true);
    try zigstd.testExpectEqual(@as(u16, 1), cvar.majorVersion());

    const data = cvar.getData();
    try zigstd.testExpectEqual(@as(u12, 1), data.tupleVariationCount().count);
    try zigstd.testExpectEqual(@as(u16, 6), data.getTupleVariationHeader(0).variationDataSize());
}
