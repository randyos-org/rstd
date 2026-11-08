const types = @import("../../types.zig");
const bytes = @import("../../../../../bytes.zig");
const zigstd = @import("../../../../../zigstd.zig");

const GlyphIndex = types.GlyphIndex;
const Offset32 = types.Offset32;

/// The format-14 header. length/num_var_selector_records need an explicit
/// align(2): the wire layout packs them right after the u16 format field with
/// no padding, but extern struct layout otherwise follows C alignment rules and
/// would insert 2 bytes before a natively-4-byte-aligned u32 field, silently
/// shifting every offset after it.
const StreamLayout = extern struct {
    format: u16,
    length: u32 align(2),
    num_var_selector_records: u32 align(2),
};

const records_start = @sizeOf(StreamLayout);

/// A u24 (a Unicode scalar value, or here also a variation selector, which
/// is itself just a Unicode scalar value used in that role) as it actually
/// sits in the stream. Zig doesn't allow u24 as an extern struct field
/// directly (only power-of-two-bit integers are extern-compatible), so every
/// 24-bit field below is stored as this raw byte triple instead (like
/// _encoding_id on EncRecord, prefer the type's accessor method over
/// touching the field directly) and converted to a real u24 only there.
const UnicodeValue = [3]u8;

/// One VariationSelectorRecord. default_uvs_offset needs align(1) since it
/// otherwise wants 4-byte alignment, which would insert a padding byte after
/// the 3-byte field before it; the wire format has no such gap.
pub const VariationSelectorRecord = extern struct {
    _var_selector: UnicodeValue,
    default_uvs_offset: Offset32 align(1),
    non_default_uvs_offset: Offset32 align(1),

    pub fn varSelector(self: @This()) u24 {
        return bytes.getU24BE(&self._var_selector);
    }
};

/// One UnicodeRange (Default UVS table entry). additional_count needs no
/// align() override since u8 has no alignment requirement to violate.
pub const UnicodeRange = extern struct {
    _start_unicode_value: UnicodeValue,
    additional_count: u8,

    pub fn startUnicodeValue(self: @This()) u24 {
        return bytes.getU24BE(&self._start_unicode_value);
    }
};

/// One UvsMapping (Non-Default UVS table entry). glyph_id can be a real
/// GlyphIndex field directly (u16-backed, so extern-compatible), but still
/// needs align(1) for the same reason default_uvs_offset does above.
pub const UvsMapping = extern struct {
    _unicode_value: UnicodeValue,
    glyph_id: GlyphIndex align(1),

    pub fn unicodeValue(self: @This()) u24 {
        return bytes.getU24BE(&self._unicode_value);
    }
};

const record_size = @sizeOf(VariationSelectorRecord);
const range_size = @sizeOf(UnicodeRange);
const mapping_size = @sizeOf(UvsMapping);

// Unicode Variation Sequences: maps (codepoint, variation selector) pairs to
// glyphs, rather than codepoints alone. Structurally unlike every other cmap
// format, so it doesn't fit the getGlyphIndex(code) shape the others share.
//
// Every field here (the variation selector records, and each record's
// default/non-default UVS tables (which live at their own byte offsets
// elsewhere in table_data, much like CmapTable's own encoding records point
// at subtables) is read directly from table_data on each call rather than
// being decoded into owned arrays up front.
pub const Format14 = struct {
    const Self = @This();

    table_data: []const u8,
    num_records: u32,

    pub fn init(table_data: []const u8, check: bool) !Self {
        if (check) try bytes.checkSize(table_data, 0, StreamLayout);
        const num_records = bytes.getU32BE(table_data[@offsetOf(StreamLayout, "num_var_selector_records")..]);
        if (check) try bytes.checkArray(table_data, records_start, num_records, VariationSelectorRecord);
        return .{ .table_data = table_data, .num_records = num_records };
    }

    pub inline fn format(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "format")..]);
    }

    pub inline fn length(self: Self) u32 {
        return bytes.getU32BE(self.table_data[@offsetOf(StreamLayout, "length")..]);
    }

    pub inline fn numRecords(self: Self) u32 {
        return self.num_records;
    }

    pub fn getRecord(self: Self, index: usize) VariationSelectorRecord {
        return bytes.typeBEToNative(VariationSelectorRecord, self.table_data[records_start + index * record_size ..]);
    }

    // Records are sorted numerically by var_selector per spec.
    pub fn findRecord(self: Self, selector: u24) ?VariationSelectorRecord {
        var lo: usize = 0;
        var hi: usize = self.num_records;
        while (lo < hi) {
            const mid = lo + (hi - lo) / 2;
            const record = self.getRecord(mid);
            const mid_selector = record.varSelector();
            if (selector < mid_selector) {
                hi = mid;
            } else if (selector > mid_selector) {
                lo = mid + 1;
            } else {
                return record;
            }
        }
        return null;
    }

    /// table_offset comes from a record's default_uvs_offset (0 means "no
    /// default UVS table"; callers only reach here when it's non-zero).
    /// Returns 0 if table_offset itself doesn't fit table_data, so a
    /// malformed offset just yields an empty range list rather than an
    /// out-of-bounds read.
    pub fn numDefaultUvsRanges(self: Self, table_offset: u32) u32 {
        if (table_offset + @sizeOf(u32) > self.table_data.len) return 0;
        return bytes.getU32BE(self.table_data[table_offset..]);
    }

    pub fn getDefaultUvsRange(self: Self, table_offset: u32, index: usize) UnicodeRange {
        const off = table_offset + @sizeOf(u32) + index * range_size;
        return bytes.typeBEToNative(UnicodeRange, self.table_data[off..]);
    }

    /// table_offset comes from a record's non_default_uvs_offset (0 means
    /// "no non-default UVS table"). Same malformed-offset handling as
    /// numDefaultUvsRanges().
    pub fn numNonDefaultUvsMappings(self: Self, table_offset: u32) u32 {
        if (table_offset + @sizeOf(u32) > self.table_data.len) return 0;
        return bytes.getU32BE(self.table_data[table_offset..]);
    }

    pub fn getNonDefaultUvsMapping(self: Self, table_offset: u32, index: usize) UvsMapping {
        const off = table_offset + @sizeOf(u32) + index * mapping_size;
        return bytes.typeBEToNative(UvsMapping, self.table_data[off..]);
    }

    pub const VariantResult = union(enum) {
        /// An explicit glyph from the non-default UVS table.
        glyph: GlyphIndex,
        /// Codepoint falls in a default UVS range: fall back to the font's
        /// normal cmap lookup for the base codepoint.
        use_default,
    };

    /// Looks up how `codepoint` should render with `selector` applied.
    /// Returns null if this subtable has no record for `selector` at all.
    pub fn getGlyphVariant(self: Self, codepoint: u32, selector: u32) ?VariantResult {
        // Both fields are 24-bit in the stream; anything wider than that can't
        // be a valid Unicode scalar value or variation selector anyway.
        if (codepoint > 0xFFFFFF or selector > 0xFFFFFF) return null;
        const cp: u24 = @intCast(codepoint);

        const record = self.findRecord(@intCast(selector)) orelse return null;

        if (record.non_default_uvs_offset != 0) {
            const table_offset = record.non_default_uvs_offset;
            const num_mappings = self.numNonDefaultUvsMappings(table_offset);
            var i: u32 = 0;
            while (i < num_mappings) : (i += 1) {
                const off = table_offset + @sizeOf(u32) + @as(usize, i) * mapping_size;
                if (off + mapping_size > self.table_data.len) break;
                const mapping = self.getNonDefaultUvsMapping(table_offset, i);
                if (mapping.unicodeValue() == cp) return .{ .glyph = mapping.glyph_id };
            }
        }

        if (record.default_uvs_offset != 0) {
            const table_offset = record.default_uvs_offset;
            const num_ranges = self.numDefaultUvsRanges(table_offset);
            var i: u32 = 0;
            while (i < num_ranges) : (i += 1) {
                const off = table_offset + @sizeOf(u32) + @as(usize, i) * range_size;
                if (off + range_size > self.table_data.len) break;
                const range = self.getDefaultUvsRange(table_offset, i);
                const start = range.startUnicodeValue();
                const range_end = @as(u32, start) + range.additional_count;
                if (cp >= start and cp <= range_end) return .use_default;
            }
        }

        return null;
    }
};

fn buildFormat14(
    alloc: zigstd.Allocator,
    selector: u24,
    default_ranges: []const struct { start_unicode_value: u24, additional_count: u8 },
    non_default_mappings: []const struct { unicode_value: u24, glyph_id: GlyphIndex },
) ![]u8 {
    const default_size: usize = if (default_ranges.len > 0) @sizeOf(u32) + default_ranges.len * range_size else 0;
    const non_default_size: usize = if (non_default_mappings.len > 0) @sizeOf(u32) + non_default_mappings.len * mapping_size else 0;

    const default_offset = records_start + record_size;
    const non_default_offset = default_offset + default_size;
    const total = non_default_offset + non_default_size;

    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);

    bytes.setU16BE(buf, 14); // format
    bytes.setU32BE(buf[@offsetOf(StreamLayout, "length")..], @intCast(total));
    bytes.setU32BE(buf[@offsetOf(StreamLayout, "num_var_selector_records")..], 1);

    const default_offset_pos = @offsetOf(VariationSelectorRecord, "default_uvs_offset");
    const non_default_offset_pos = @offsetOf(VariationSelectorRecord, "non_default_uvs_offset");
    bytes.setU24BE(buf[records_start..], selector);
    bytes.setU32BE(buf[records_start + default_offset_pos ..], if (default_ranges.len > 0) @intCast(default_offset) else 0);
    bytes.setU32BE(buf[records_start + non_default_offset_pos ..], if (non_default_mappings.len > 0) @intCast(non_default_offset) else 0);

    const additional_count_pos = @offsetOf(UnicodeRange, "additional_count");
    if (default_ranges.len > 0) {
        bytes.setU32BE(buf[default_offset..], @intCast(default_ranges.len));
        var off: usize = default_offset + @sizeOf(u32);
        for (default_ranges) |r| {
            bytes.setU24BE(buf[off..], r.start_unicode_value);
            buf[off + additional_count_pos] = r.additional_count;
            off += range_size;
        }
    }

    const glyph_id_pos = @offsetOf(UvsMapping, "glyph_id");
    if (non_default_mappings.len > 0) {
        bytes.setU32BE(buf[non_default_offset..], @intCast(non_default_mappings.len));
        var off: usize = non_default_offset + @sizeOf(u32);
        for (non_default_mappings) |m| {
            bytes.setU24BE(buf[off..], m.unicode_value);
            bytes.setU16BE(buf[off + glyph_id_pos ..], @backingInt(m.glyph_id));
            off += mapping_size;
        }
    }

    return buf;
}

test "Format14.getGlyphVariant finds an explicit non-default glyph mapping" {
    const alloc = zigstd.testAllocator;
    const buf = try buildFormat14(
        alloc,
        0xFE0F,
        &.{},
        &.{.{ .unicode_value = 0x2603, .glyph_id = @fromBackingInt(@intCast(42)) }},
    );
    defer alloc.free(buf);

    const fmt14 = try Format14.init(buf, true);

    const result = fmt14.getGlyphVariant(0x2603, 0xFE0F).?;
    try zigstd.testExpectEqual(@as(u16, 42), @backingInt(result.glyph));
}

test "Format14.getGlyphVariant falls back to the default cmap inside a default range" {
    const alloc = zigstd.testAllocator;
    const buf = try buildFormat14(
        alloc,
        0xFE0E,
        &.{.{ .start_unicode_value = 0x2600, .additional_count = 10 }},
        &.{},
    );
    defer alloc.free(buf);

    const fmt14 = try Format14.init(buf, true);

    try zigstd.testExpectEqual(Format14.VariantResult.use_default, fmt14.getGlyphVariant(0x2605, 0xFE0E).?);
    try zigstd.testExpectEqual(@as(?Format14.VariantResult, null), fmt14.getGlyphVariant(0x2611, 0xFE0E));
}

test "Format14.getGlyphVariant returns null for an unrecorded variation selector" {
    const alloc = zigstd.testAllocator;
    const buf = try buildFormat14(alloc, 0xFE0F, &.{}, &.{});
    defer alloc.free(buf);

    const fmt14 = try Format14.init(buf, true);

    try zigstd.testExpectEqual(@as(?Format14.VariantResult, null), fmt14.getGlyphVariant(0x2603, 0xFE0E));
}

test "Format14.getGlyphVariant returns null past the 24-bit range either argument can represent" {
    const alloc = zigstd.testAllocator;
    const buf = try buildFormat14(alloc, 0xFE0F, &.{}, &.{});
    defer alloc.free(buf);

    const fmt14 = try Format14.init(buf, true);

    try zigstd.testExpectEqual(@as(?Format14.VariantResult, null), fmt14.getGlyphVariant(0x1000000, 0xFE0F));
    try zigstd.testExpectEqual(@as(?Format14.VariantResult, null), fmt14.getGlyphVariant(0x2603, 0x1000000));
}

test "Format14 record/range/mapping layouts match the OpenType stream layout" {
    try zigstd.testExpectEqual(@as(usize, 10), records_start);
    try zigstd.testExpectEqual(@as(usize, 11), @sizeOf(VariationSelectorRecord));
    try zigstd.testExpectEqual(@as(usize, 3), @offsetOf(VariationSelectorRecord, "default_uvs_offset"));
    try zigstd.testExpectEqual(@as(usize, 7), @offsetOf(VariationSelectorRecord, "non_default_uvs_offset"));
    try zigstd.testExpectEqual(@as(usize, 4), @sizeOf(UnicodeRange));
    try zigstd.testExpectEqual(@as(usize, 3), @offsetOf(UnicodeRange, "additional_count"));
    try zigstd.testExpectEqual(@as(usize, 5), @sizeOf(UvsMapping));
    try zigstd.testExpectEqual(@as(usize, 3), @offsetOf(UvsMapping, "glyph_id"));
}
