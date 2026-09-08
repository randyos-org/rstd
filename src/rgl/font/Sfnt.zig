const bytes = @import("../../bytes.zig");
const zigstd = @import("../../zigstd.zig");
const tables = @import("tables/__root__.zig");

const Allocator = zigstd.Allocator;
const typeBEToNative = bytes.typeBEToNative;
const allocTypeBEToNative = bytes.allocTypeBEToNative;
const OffsetTable = tables.general.OffsetTable;
const TableRecord = tables.general.TableRecord;
const TableTag = tables.general.TableTag;
const ScalerType = tables.general.ScalerType;

const Self = @This();

/// Naming table.  Required.
// name,
/// Character to glyph mapping.  Required.
cmap: ?*tables.CmapTable = null,
/// Font header. May be missing in bitmap-only fonts, otherwise required.
head: ?*tables.HeadTable = null,
/// Horizontal header. May be missing in bitmap-only fonts, otherwise required.
hhea: ?*tables.HheaTable = null,
/// Horizontal metrics. May be missing in bitmap-only fonts, otherwise required.
// hmtx,
/// Maximum profile. Optional per TTF, required per OpenType spec.
maxp: ?*tables.MaxpTable = null,
/// PostScript table. Required.
// post,

// ---
// Classic glyph-based TrueType tables. Only required if not CFF.

/// Glyph data. Only required if not CFF.
// glyf,
/// Index to location. Only required if glyf is present
// loca,

// ---
// CFF TABLES: only required if not glyf

/// Compact Font Format.
// CFF,

// ---
// OPTIONAL TTF TABLES:

/// kerning
// kern,

// ---
// Advanced typographic tables

/// Glyph positioning data
// GPOS,

pub const SfntError = error{
    MissingRequiredTable,
    IndexMapMissing,
} || TableTag.Error || ScalerType.Error || bytes.BytesError || Allocator.Error; // || CffData.InitError;

const CreateOptions = struct {
    check: bool = false,
};

pub fn create(alloc: Allocator, font_data: []const u8, options: CreateOptions) SfntError!Self {
    if (options.check) try bytes.checkSize(font_data, 0, OffsetTable);
    const offset_table = typeBEToNative(OffsetTable, font_data);
    if (options.check) _ = try ScalerType.fromInt(offset_table.scaler);

    const tables_start = @sizeOf(OffsetTable);
    const num_tables = offset_table.num_tables;
    if (options.check) try bytes.checkArray(font_data, tables_start, num_tables, TableRecord);
    const table_recs = try bytes.allocArrayBEToNative(alloc, TableRecord, num_tables, font_data[tables_start..]);
    defer alloc.free(table_recs);

    var self: Self = .{};
    errdefer self.destroy(alloc);
    for (table_recs) |record| {
        const table_data = try bytes.safeSlice(font_data, record.offset, record.length, options.check);
        const table_tag = try TableTag.fromTagStr(record.tag);
        switch (table_tag) {
            .cmap => self.cmap = try tables.CmapTable.create(alloc, table_data, options.check),
            .head => self.head = try tables.HeadTable.create(alloc, table_data, options.check),
            .hhea => self.hhea = try tables.HheaTable.create(alloc, table_data, options.check),
            .maxp => self.maxp = try tables.MaxpTable.create(alloc, table_data, options.check),
            .hmtx => {},
            .loca => {},
            .glyf => {},
            .kern => {},
            .GPOS => {},
            .CFF => {},
            else => {},
        }
    }

    return self;
}

pub fn destroy(self: *Self, alloc: Allocator) void {
    if (self.cmap) |cmap| cmap.destroy(alloc);
    if (self.head) |head| head.destroy(alloc);
    if (self.hhea) |hhea| hhea.destroy(alloc);
    if (self.maxp) |maxp| maxp.destroy(alloc);
    // if (self.hmtx) |hmtx| hmtx.destroy(alloc);
    // if (self.loca) |loca| loca.destroy(alloc);
    // if (self.glyf) |glyf| glyf.destroy(alloc);
    // if (self.kern) |kern| kern.destroy(alloc);
    // if (self.GPOS) |GPOS| GPOS.destroy(alloc);
    // if (self.cff) |cff| cff.destroy(alloc);
}

fn setU16BE(dst: []u8, offset: usize, v: u16) void {
    bytes.writeIntBE(u16, dst[offset..][0..2], v);
}

fn setU32BE(dst: []u8, offset: usize, v: u32) void {
    bytes.writeIntBE(u32, dst[offset..][0..4], v);
}

test "Sfnt.create rejects a font with a bad scaler tag" {
    var buf: [12]u8 = @splat(0);
    setU32BE(&buf, 0, 0xDEADBEEF);
    try zigstd.testExpectError(error.UnsupportedScalerType, create(zigstd.testAllocator, buf[0..], .{ .check = true }));
}

test "Sfnt.create leaks nothing when a later table's allocation fails" {
    const offset_table_size = 12;
    const rec_size = 16;
    const head_start = offset_table_size + 2 * rec_size; // 44
    const head_size = 54;
    const hhea_start = head_start + head_size; // 98
    const hhea_size = 36;
    const total = hhea_start + hhea_size;

    const buf = try zigstd.testAllocator.alloc(u8, total);
    defer zigstd.testAllocator.free(buf);
    @memset(buf, 0);

    setU32BE(buf, 0, 0x00010000); // scaler
    setU16BE(buf, 4, 2); // num_tables

    @memcpy(buf[offset_table_size..][0..4], "head");
    setU32BE(buf, offset_table_size + 8, head_start);
    setU32BE(buf, offset_table_size + 12, head_size);

    @memcpy(buf[offset_table_size + rec_size ..][0..4], "hhea");
    setU32BE(buf, offset_table_size + rec_size + 8, hhea_start);
    setU32BE(buf, offset_table_size + rec_size + 12, hhea_size);

    var found_failure_point = false;
    var fail_at: usize = 0;
    while (fail_at < 6) : (fail_at += 1) {
        var failing = zigstd.TestFailingAllocator.init(zigstd.testAllocator, .{ .fail_index = fail_at });
        var sfnt = create(failing.allocator(), buf, .{ .check = true }) catch {
            found_failure_point = true;
            continue;
        };
        sfnt.destroy(failing.allocator());
    }
    try zigstd.testExpectTrue(found_failure_point);
}
