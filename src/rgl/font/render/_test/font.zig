//! A synthetic sfnt for the render tests: units_per_em 1000; glyph 0
//! (.notdef) empty, glyph 1 a 200x300-funit box (advance 250) mapped from
//! 'A' to 'Z' and 'É', glyph 2 an outline-less 150-funit space mapped from ' '.

const bytes = @import("../../../../bytes.zig");
const general = @import("../../sfnt/tables/general.zig");
const zigstd = @import("../../../../zigstd.zig");

const Allocator = zigstd.Allocator;
const TableTag = general.TableTag;

fn appendTable(buf: []u8, table_dir_off: usize, table_idx: usize, tag: TableTag, data: []const u8, off: usize) void {
    const rec_off = table_dir_off + table_idx * 16;
    bytes.setU32BE(buf[rec_off..], tag.toInt());
    bytes.setU32BE(buf[rec_off + 4 ..], 0); // checksum, unused
    bytes.setU32BE(buf[rec_off + 8 ..], @intCast(off));
    bytes.setU32BE(buf[rec_off + 12 ..], @intCast(data.len));
    @memcpy(buf[off..][0..data.len], data);
}

pub fn build(alloc: Allocator) ![]u8 {
    const num_tables = 7;
    const table_dir_off = 12;
    const tables_start = table_dir_off + num_tables * 16;

    var head_buf: [54]u8 = @splat(0);
    bytes.setU16BE(head_buf[18..], 1000); // units_per_em
    bytes.setI16BE(head_buf[50..], 1); // index_to_loc_format: long

    var hhea_buf: [36]u8 = @splat(0);
    bytes.setI16BE(hhea_buf[4..], 800); // ascent
    bytes.setI16BE(hhea_buf[6..], -200); // descent
    bytes.setI16BE(hhea_buf[8..], 0); // line_gap
    bytes.setU16BE(hhea_buf[34..], 3); // num_of_long_hor_metrics

    var maxp_buf: [6]u8 = @splat(0);
    bytes.setU32BE(maxp_buf[0..], 0x00005000); // version 0.5
    bytes.setU16BE(maxp_buf[4..], 3); // num_glyphs

    // hmtx: glyph0 (unused/.notdef), glyph1 (box, advance 250), glyph2 (space, advance 150)
    var hmtx_buf: [12]u8 = @splat(0);
    bytes.setU16BE(hmtx_buf[4..], 250);
    bytes.setU16BE(hmtx_buf[8..], 150);

    // loca (long): glyph0 empty [0,0), glyph1 box [0,34), glyph2 empty [34,34)
    var loca_buf: [16]u8 = @splat(0);
    bytes.setU32BE(loca_buf[0..], 0);
    bytes.setU32BE(loca_buf[4..], 0);
    bytes.setU32BE(loca_buf[8..], 34);
    bytes.setU32BE(loca_buf[12..], 34);

    // glyf: one simple glyph, 1 contour, 4 on-curve points forming a
    // 200x300 rectangle, using the "long" x/y coordinate encoding (flag
    // 0x01) throughout so deltas aren't limited to a single byte.
    var glyf_buf: [34]u8 = @splat(0);
    bytes.setI16BE(glyf_buf[0..], 1); // number_of_contours
    bytes.setI16BE(glyf_buf[2..], 0); // x_min
    bytes.setI16BE(glyf_buf[4..], 0); // y_min
    bytes.setI16BE(glyf_buf[6..], 200); // x_max
    bytes.setI16BE(glyf_buf[8..], 300); // y_max
    bytes.setU16BE(glyf_buf[10..], 3); // endPtsOfContours[0] = 3 (4 points)
    bytes.setU16BE(glyf_buf[12..], 0); // instructionLength
    const flags_start = 14;
    @memset(glyf_buf[flags_start..][0..4], 0x01); // on_curve, x/y long
    const x_start = flags_start + 4;
    bytes.setI16BE(glyf_buf[x_start + 0 ..], 0); // p0 dx
    bytes.setI16BE(glyf_buf[x_start + 2 ..], 200); // p1 dx
    bytes.setI16BE(glyf_buf[x_start + 4 ..], 0); // p2 dx
    bytes.setI16BE(glyf_buf[x_start + 6 ..], -200); // p3 dx
    const y_start = x_start + 8;
    bytes.setI16BE(glyf_buf[y_start + 0 ..], 0); // p0 dy
    bytes.setI16BE(glyf_buf[y_start + 2 ..], 0); // p1 dy
    bytes.setI16BE(glyf_buf[y_start + 4 ..], 300); // p2 dy
    bytes.setI16BE(glyf_buf[y_start + 6 ..], 0); // p3 dy

    // cmap: format 0, one Windows-Unicode-BMP record. 'A'-'Z' and 'É' -> glyph 1,
    // ' ' -> glyph 2, everything else unmapped (.notdef).
    var cmap_buf: [4 + 8 + 262]u8 = @splat(0);
    bytes.setU16BE(cmap_buf[0..], 0); // version
    bytes.setU16BE(cmap_buf[2..], 1); // num_tables
    bytes.setU16BE(cmap_buf[4..], 3); // platform_id: microsoft
    bytes.setU16BE(cmap_buf[6..], 1); // encoding_id: unicode_bmp
    bytes.setU32BE(cmap_buf[8..], 12); // subtable offset
    const sub = cmap_buf[12..];
    bytes.setU16BE(sub, 0); // format
    bytes.setU16BE(sub[2..], 262); // length
    bytes.setU16BE(sub[4..], 0); // language
    var c: u8 = 'A';
    while (c <= 'Z') : (c += 1) sub[6 + @as(usize, c)] = 1;
    sub[6 + @as(usize, ' ')] = 2;
    sub[6 + 0xC9] = 1; // 'É', two bytes in UTF-8

    var buf = try alloc.alloc(u8, tables_start + head_buf.len + hhea_buf.len + maxp_buf.len + hmtx_buf.len + loca_buf.len + cmap_buf.len + glyf_buf.len);
    errdefer alloc.free(buf);
    @memset(buf, 0);

    bytes.setU32BE(buf[0..], 0x00010000); // scaler
    bytes.setU16BE(buf[4..], num_tables);

    var off: usize = tables_start;
    appendTable(buf, table_dir_off, 0, .head, &head_buf, off);
    off += head_buf.len;
    appendTable(buf, table_dir_off, 1, .hhea, &hhea_buf, off);
    off += hhea_buf.len;
    appendTable(buf, table_dir_off, 2, .maxp, &maxp_buf, off);
    off += maxp_buf.len;
    appendTable(buf, table_dir_off, 3, .hmtx, &hmtx_buf, off);
    off += hmtx_buf.len;
    appendTable(buf, table_dir_off, 4, .loca, &loca_buf, off);
    off += loca_buf.len;
    appendTable(buf, table_dir_off, 5, .cmap, &cmap_buf, off);
    off += cmap_buf.len;
    appendTable(buf, table_dir_off, 6, .glyf, &glyf_buf, off);
    off += glyf_buf.len;

    return buf;
}
