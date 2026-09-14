const types = @import("../stream_types.zig");
const bytes = @import("../../../bytes.zig");
const zigstd = @import("../../../zigstd.zig");

const HeaderLayout = extern struct {
    version: types.Version16Dot16 align(2),
    num_glyphs: u16,
};

const TtfLayout = extern struct {
    max_points: u16,
    max_contours: u16,
    max_component_points: u16,
    max_component_contours: u16,
    max_zones: u16,
    max_twilight_points: u16,
    max_storage: u16,
    max_function_defs: u16,
    max_instruction_defs: u16,
    max_stack_elements: u16,
    max_size_of_instructions: u16,
    max_component_elements: u16,
    max_component_depth: u16,
};

const ttf_start = @sizeOf(HeaderLayout);

pub const MaxpTable = struct {
    const Self = @This();

    table_data: []const u8,

    pub fn init(table_data: []const u8, check: bool) !Self {
        if (check) try bytes.checkSize(table_data, 0, HeaderLayout);
        const self: Self = .{ .table_data = table_data };

        // Version 0.5 fonts have no TTF-specific section; anything else is
        // expected to carry the full 32-byte table.
        if (check and @as(u32, @bitCast(self.version())) != 0x00005000) {
            try bytes.checkSize(table_data, ttf_start, TtfLayout);
        }

        return self;
    }

    pub inline fn version(self: Self) types.Version16Dot16 {
        return types.Version16Dot16.fromBytes(self.table_data[@offsetOf(HeaderLayout, "version")..]);
    }

    pub inline fn numGlyphs(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(HeaderLayout, "num_glyphs")..]);
    }

    /// Valid only when version() is 1.0 (0x00010000); version 0.5 fonts have
    /// no TTF-specific section, so calling these on one reads past the table.
    pub inline fn maxPoints(self: Self) u16 {
        return bytes.getU16BE(self.table_data[ttf_start + @offsetOf(TtfLayout, "max_points")..]);
    }

    pub inline fn maxContours(self: Self) u16 {
        return bytes.getU16BE(self.table_data[ttf_start + @offsetOf(TtfLayout, "max_contours")..]);
    }

    pub inline fn maxComponentPoints(self: Self) u16 {
        return bytes.getU16BE(self.table_data[ttf_start + @offsetOf(TtfLayout, "max_component_points")..]);
    }

    pub inline fn maxComponentContours(self: Self) u16 {
        return bytes.getU16BE(self.table_data[ttf_start + @offsetOf(TtfLayout, "max_component_contours")..]);
    }

    pub inline fn maxZones(self: Self) u16 {
        return bytes.getU16BE(self.table_data[ttf_start + @offsetOf(TtfLayout, "max_zones")..]);
    }

    pub inline fn maxTwilightPoints(self: Self) u16 {
        return bytes.getU16BE(self.table_data[ttf_start + @offsetOf(TtfLayout, "max_twilight_points")..]);
    }

    pub inline fn maxStorage(self: Self) u16 {
        return bytes.getU16BE(self.table_data[ttf_start + @offsetOf(TtfLayout, "max_storage")..]);
    }

    pub inline fn maxFunctionDefs(self: Self) u16 {
        return bytes.getU16BE(self.table_data[ttf_start + @offsetOf(TtfLayout, "max_function_defs")..]);
    }

    pub inline fn maxInstructionDefs(self: Self) u16 {
        return bytes.getU16BE(self.table_data[ttf_start + @offsetOf(TtfLayout, "max_instruction_defs")..]);
    }

    pub inline fn maxStackElements(self: Self) u16 {
        return bytes.getU16BE(self.table_data[ttf_start + @offsetOf(TtfLayout, "max_stack_elements")..]);
    }

    pub inline fn maxSizeOfInstructions(self: Self) u16 {
        return bytes.getU16BE(self.table_data[ttf_start + @offsetOf(TtfLayout, "max_size_of_instructions")..]);
    }

    pub inline fn maxComponentElements(self: Self) u16 {
        return bytes.getU16BE(self.table_data[ttf_start + @offsetOf(TtfLayout, "max_component_elements")..]);
    }

    pub inline fn maxComponentDepth(self: Self) u16 {
        return bytes.getU16BE(self.table_data[ttf_start + @offsetOf(TtfLayout, "max_component_depth")..]);
    }
};

test "maxp layouts match the OpenType stream layout (32 bytes total, no padding)" {
    try zigstd.testExpectEqual(@as(usize, 6), @sizeOf(HeaderLayout));
    try zigstd.testExpectEqual(@as(usize, 4), @offsetOf(HeaderLayout, "num_glyphs"));
    try zigstd.testExpectEqual(@as(usize, 26), @sizeOf(TtfLayout));
    try zigstd.testExpectEqual(@as(usize, 24), @offsetOf(TtfLayout, "max_component_depth"));
}

test "MaxpTable getters read a real big-endian 'maxp' table correctly" {
    var buf: [32]u8 = @splat(0);
    bytes.setU16BE(buf[4..], 500); // num_glyphs
    bytes.setU16BE(buf[30..], 2); // max_component_depth (ttf offset 24)

    const maxp = try MaxpTable.init(buf[0..], true);

    try zigstd.testExpectEqual(@as(u16, 500), maxp.numGlyphs());
    try zigstd.testExpectEqual(@as(u16, 2), maxp.maxComponentDepth());
}

test "MaxpTable.init validates only the header for version 0.5 (no TTF fields)" {
    var buf: [6]u8 = @splat(0);
    bytes.setU32BE(buf[0..], 0x00005000); // version 0.5
    bytes.setU16BE(buf[4..], 12); // num_glyphs

    const maxp = try MaxpTable.init(buf[0..], true);

    try zigstd.testExpectEqual(@as(u16, 12), maxp.numGlyphs());
}
