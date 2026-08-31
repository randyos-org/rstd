const types = @import("types.zig");
const Fixed = types.Fixed;

const bytes = @import("../../bytes.zig");
const zigstd = @import("../../zigstd.zig");

const readIntBE = bytes.readIntBE;
const writeIntBE = bytes.writeIntBE;
const Allocator = zigstd.Allocator;

const OffsetTable = packed struct {
    scaler: u32,
    num_tables: u16,
    search_range: u16,
    entry_selector: u16,
    range_shift: u16,
};

const TableDirectoryEntry = extern struct {
    tag: [4]u8,
    check_sum: u32,
    offset: u32,
    length: u32,
};

const HeadTable = packed struct {
    version: Fixed,
    font_revision: Fixed,
    check_sum_adjustment: u32,
    magic_number: u32,
    flags: u16,
    units_per_em: u16,
    created: i64,
    modified: i64,
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

const MaxpTable = packed struct {
    version: Fixed,
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
    maxStackElements: u16,
    maxSizeOfInstructions: u16,
    maxComponentElements: u16,
    maxComponentDepth: u16,
};

const HheaTable = packed struct {
    version: Fixed,
    ascent: i16,
    descent: i16,
    line_gap: i16,
    advance_width_max: u16,
    min_left_side_bearing: i16,
    min_right_side_bearing: i16,
    x_max_extent: i16,
    caret_slope_rise: i16,
    caret_slope_run: i16,
    caret_offset: i16,
    reserved1: i16,
    reserved2: i16,
    reserved3: i16,
    reserved4: i16,
    metric_data_format: i16,
    num_of_long_hor_metrics: u16,
};

pub const HmtxTable = struct {
    hmtx_bytes: []const u8,

    pub const LongHorMetric = packed struct {
        advance_width: u16,
        left_side_bearing: i16,
    };

    // pub fn getMetrics(self: HmtxTable, num_hor_metrics: usize, glyph_index: usize) LongHorMetric {
    //     if (glyph_index < num_hor_metrics) {
    //         return self.loadHorMetric(glyph_index);
    //     } else {
    //         const last = self.loadHorMetric(num_hor_metrics - 1);
    //         const lsb_index = glyph_index - num_hor_metrics;
    //         const lsb_offs = num_hor_metrics * @bitSizeOf(LongHorMetric) / 8 + lsb_index * 2;
    //         const lsb = fixEndianness(std.mem.bytesToValue(i16, self.hmtx_bytes[lsb_offs..]));
    //         return .{
    //             .advance_width = last.advance_width,
    //             .left_side_bearing = lsb,
    //         };
    //     }
    // }

    // fn loadHorMetric(self: HmtxTable, idx: usize) LongHorMetric {
    //     const offs = idx * @bitSizeOf(LongHorMetric) / 8;
    //     return fixEndianness(std.mem.bytesToValue(LongHorMetric, self.hmtx_bytes[offs..]));
    // }
};

pub const CmapTable = struct {
    cmap_bytes: []const u8,

    const Index = packed struct {
        version: u16,
        num_subtables: u16,
    };

    const SubtableLookup = packed struct {
        platform_id: u16,
        platform_specific_id: u16,
        offset: u32,

        fn isUnicodeBmp(self: SubtableLookup) bool {
            return (self.platform_id == 0 and self.platform_specific_id == 3) or // unicode + bmp
                (self.platform_id == 3 and self.platform_specific_id == 1) // windows + unicode ucs 2
            ;
        }
    };

    const Subtable = packed struct {
        platform_id: u16,
        platform_specific_id: u16,
        offset: u32,
    };

    pub const SubtableFormat4 = struct {
        format: u16,
        length: u16,
        language: u16,
        seg_count_x2: u16,
        search_range: u16,
        entry_selector: u16,
        range_shift: u16,
        end_code: []const u16,
        reserved_pad: u16,
        start_code: []const u16,
        id_delta: []const u16,
        id_range_offset: []const u16,
        glyph_indices: []const u16,

        fn getGlyphIndex(self: SubtableFormat4, c: u16) u16 {
            // This won't make sense if you don't read the spec...
            var i: usize = 0;
            while (i < self.end_code.len) {
                if (self.end_code[i] >= c and self.start_code[i] <= c) {
                    break;
                }
                i += 1;
            }

            if (i >= self.end_code.len) return 0;

            const byte_offset_from_id_offset = self.id_range_offset[i];
            if (byte_offset_from_id_offset == 0) {
                return self.id_delta[i] +% c;
            } else {
                // We apply the pointer offset a little different than the spec
                // suggests. We made individual allocations when copying/byte
                // swapping the id_range_offset and glyph_indices out of the input
                // data. This means that we can't just do pointer addition
                //
                // Instead we look at the data as follows
                //
                // [ id range ] [glyph indices ]
                //     |--offs_bytes--|
                //     ^
                //     i
                //
                // To find the index into glyph indices, we just subtract i from
                // id_range.len, and subtract that from the offset
                const offs_from_loc = byte_offset_from_id_offset / 2 + (c - self.start_code[i]);
                const dist_to_end = self.id_range_offset.len - i;
                const glyph_index_index = offs_from_loc - dist_to_end;
                return self.glyph_indices[glyph_index_index] +% self.id_delta[i];
            }
        }
    };

    // fn readIndex(self: CmapTable) Index {
    //     return fixEndianness(std.mem.bytesToValue(Index, self.cmap_bytes[0 .. @bitSizeOf(Index) / 8]));
    // }

    // fn readSubtableLookup(self: CmapTable, idx: usize) SubtableLookup {
    //     const subtable_size = @bitSizeOf(SubtableLookup) / 8;
    //     const start = @bitSizeOf(Index) / 8 + idx * subtable_size;
    //     const end = start + subtable_size;

    //     return fixEndianness(std.mem.bytesToValue(SubtableLookup, self.cmap_bytes[start..end]));
    // }

    // fn readSubtableFormat(self: CmapTable, offset: usize) u16 {
    //     return fixEndianness(std.mem.bytesToValue(u16, self.cmap_bytes[offset .. offset + 2]));
    // }

    // fn readSubtableFormat4(self: CmapTable, alloc: Allocator, offset: usize) !SubtableFormat4 {
    //     var runtime_parser = RuntimeParser{ .data = self.cmap_bytes[offset..] };
    //     const format = runtime_parser.readVal(u16);
    //     const length = runtime_parser.readVal(u16);
    //     const language = runtime_parser.readVal(u16);
    //     const seg_count_x2 = runtime_parser.readVal(u16);
    //     const search_range = runtime_parser.readVal(u16);
    //     const entry_selector = runtime_parser.readVal(u16);
    //     const range_shift = runtime_parser.readVal(u16);

    //     const end_code: []const u16 = try runtime_parser.readArray(u16, alloc, seg_count_x2 / 2);
    //     const reserved_pad = runtime_parser.readVal(u16);
    //     const start_code: []const u16 = try runtime_parser.readArray(u16, alloc, seg_count_x2 / 2);
    //     const id_delta: []const u16 = try runtime_parser.readArray(u16, alloc, seg_count_x2 / 2);
    //     const id_range_offset: []const u16 = try runtime_parser.readArray(u16, alloc, seg_count_x2 / 2);
    //     const glyph_indices: []const u16 = try runtime_parser.readArray(u16, alloc, (runtime_parser.data.len - runtime_parser.idx) / 2);

    //     return .{
    //         .format = format,
    //         .length = length,
    //         .language = language,
    //         .seg_count_x2 = seg_count_x2,
    //         .search_range = search_range,
    //         .entry_selector = entry_selector,
    //         .range_shift = range_shift,
    //         .end_code = end_code,
    //         .reserved_pad = reserved_pad,
    //         .start_code = start_code,
    //         .id_delta = id_delta,
    //         .id_range_offset = id_range_offset,
    //         .glyph_indices = glyph_indices,
    //     };
    // }
};
