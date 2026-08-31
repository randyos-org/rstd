const std = @import("std");

const bytes = @import("../../bytes.zig");
const zigstd = @import("../../zigstd.zig");
const tables = @import("tables/__root__.zig");

const Allocator = zigstd.Allocator;
const typeBEToNative = bytes.typeBEToNative;
const OffsetTable = tables.general.OffsetTable;
const TableRecord = tables.general.TableRecord;
const TableTag = tables.general.TableTag;

const Self = @This();

offset_table: OffsetTable,
/// Naming table.  Required.
// name,
/// Character to glyph mapping.  Required.
// cmap,
/// Font header. May be missing in bitmap-only fonts, otherwise required.
head: ?tables.general.HeadTable = null,
/// Horizontal header. May be missing in bitmap-only fonts, otherwise required.
hhea: ?tables.general.HheaTable = null,
/// Horizontal metrics. May be missing in bitmap-only fonts, otherwise required.
// hmtx,
/// Maximum profile. Optional per TTF, required per OpenType spec.
maxp: ?tables.general.MaxpTable = null,
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
    /// The font file unexpectedly ended when more data was expected.
    EndOfStream,
    MissingRequiredTable,
    IndexMapMissing,
}; // || CffData.InitError;

const offset_size = @sizeOf(OffsetTable);
const table_record_size = @sizeOf(TableRecord);

pub fn init(font_data: []const u8) SfntError!Self {
    if (font_data.len < offset_size) return error.EndOfStream;
    const offset_table = typeBEToNative(OffsetTable, &font_data);
    if (offset_size + offset_table.num_tables * table_record_size > font_data.len) return error.EndOfStream;

    var self: Self = .{
        .offset_table = offset_table,
    };
    for (0..offset_table.num_tables) |i| {
        const record_offset = offset_size + i * table_record_size;
        const record = typeBEToNative(TableRecord, &font_data[record_offset..]);
        if (record.offset + record.length > font_data.len) return error.EndOfStream;
        const table_tag = try TableTag.fromTagStr(record.tag);
        switch (table_tag) {
            .cmap => {},
            .head => self.head = typeBEToNative(tables.general.HeadTable, &font_data[record.offset..]),
            .hhea => self.hhea = typeBEToNative(tables.general.HheaTable, &font_data[record.offset..]),
            .hmtx => {},
            .loca => {},
            .glyf => {},
            .maxp => self.maxp = typeBEToNative(tables.general.MaxpTable, &font_data[record.offset..]),
            .kern => {},
            .GPOS => {},
            .cff => {},
            else => {},
        }
    }

    return .{
        .raw_bytes = font_data,
        .offset_table = offset_table,
    };
}
