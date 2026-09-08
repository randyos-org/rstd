const types = @import("../types.zig");
const bytes = @import("../../../bytes.zig");
const zigstd = @import("../../../zigstd.zig");

const Allocator = zigstd.Allocator;

pub const MaxpTable = extern struct {
    const Self = @This();

    version: types.Fixed,
    num_glyphs: u16,
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

    pub fn create(alloc: Allocator, table_data: []const u8, check: bool) !*Self {
        if (check) try bytes.checkSize(table_data, 0, Self);
        return bytes.allocTypeBEToNative(alloc, Self, table_data);
    }

    pub fn destroy(self: *Self, alloc: Allocator) void {
        alloc.destroy(self);
    }
};

test "MaxpTable matches the OpenType wire layout (32 bytes, no padding)" {
    try zigstd.testExpectEqual(@as(usize, 32), @sizeOf(MaxpTable));
    try zigstd.testExpectEqual(@as(usize, 4), @offsetOf(MaxpTable, "num_glyphs"));
    try zigstd.testExpectEqual(@as(usize, 30), @offsetOf(MaxpTable, "max_component_depth"));
}

test "MaxpTable.create reads a real big-endian 'maxp' table correctly" {
    const alloc = zigstd.testAllocator;
    var buf: [32]u8 = @splat(0);
    bytes.writeIntBE(u16, buf[4..6], 500); // num_glyphs
    bytes.writeIntBE(u16, buf[30..32], 2); // max_component_depth

    const maxp = try MaxpTable.create(alloc, buf[0..], true);
    defer maxp.destroy(alloc);

    try zigstd.testExpectEqual(@as(u16, 500), maxp.num_glyphs);
    try zigstd.testExpectEqual(@as(u16, 2), maxp.max_component_depth);
}
