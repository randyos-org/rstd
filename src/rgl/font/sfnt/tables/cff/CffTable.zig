const bytes = @import("../../../../../bytes.zig");
const types = @import("../../stream_types.zig");

const GlyphIndex = types.GlyphIndex;

const Index = @import("Index.zig");
const Dict = @import("Dict.zig");

// Index/Dict lookups are O(n) walks, so init() resolves CharStrings/Subrs
// indexes once up front instead of per glyph - same idea as LocaTable
// caching num_glyphs/long_format. Charset/Encoding are unused (cmap resolves
// GIDs directly); the Type 2 charstring interpreter is deferred.

/// Schema for the CFF header, used only to compute field offsets/sizes at
/// comptime; never instantiated.
const StreamLayout = extern struct {
    major_version: u8,
    minor_version: u8,
    header_size: u8,
    off_size: u8,
};

pub const CffError = error{
    MissingCharStrings,
    MissingFdArray,
    MissingFdSelect,
    UnsupportedCharstringType,
};

/// Resolves a DICT's Private operator ([size, offset]) and, from within it,
/// its Subrs operator (relative to the Private DICT's own start - a spec
/// quirk, unlike every other CFF offset). Null means no Private entry or no
/// local Subrs, both common.
fn resolvePrivateSubrs(cff_bytes: []const u8, dict: Dict, check: bool) bytes.BytesError!?Index {
    var priv: [2]i32 = undefined;
    if (dict.getInts(.private, &priv) != 2) return null;
    if (priv[0] <= 0 or priv[1] < 0) return null;

    const priv_size: usize = @intCast(priv[0]);
    const priv_offset: usize = @intCast(priv[1]);
    if (check) try bytes.checkSlice(cff_bytes, priv_offset, priv_size);
    const private_dict: Dict = .{ .table_data = cff_bytes[priv_offset..][0..priv_size] };

    const subrs_rel = private_dict.getInt(.subrs) orelse return null;
    if (subrs_rel < 0) return null;
    const subrs_offset = priv_offset + @as(usize, @intCast(subrs_rel));
    if (check) try bytes.checkOffset(cff_bytes, subrs_offset);
    return try Index.init(cff_bytes[subrs_offset..], check);
}

/// Compact Font Format: PostScript-flavored (Type 2 charstring) glyph
/// outline data, used in place of glyf+loca on CFF-flavored .otf fonts.
const CffTable = @This();

table_data: []const u8,
charstrings_index: Index,
global_subrs_index: Index,
is_cid: bool,
/// Set when !is_cid; the single local Subrs INDEX shared by every
/// glyph. Null if this font has no local subroutines at all.
local_subrs_index: ?Index,
/// Set when is_cid; the Font DICT INDEX, one entry per FD, each
/// glyph's own Private/Subrs selected via fdselect_offset.
fd_array_index: ?Index,
/// Set when is_cid; absolute offset of the FDSelect table.
fdselect_offset: ?usize,

pub fn init(table_data: []const u8, check: bool) (bytes.BytesError || CffError)!CffTable {
    if (check) try bytes.checkSize(table_data, 0, StreamLayout);

    var pos: usize = table_data[@offsetOf(StreamLayout, "header_size")];
    const name_index = try Index.init(table_data[pos..], check);
    pos += name_index.byteSize();
    const top_dict_index = try Index.init(table_data[pos..], check);
    pos += top_dict_index.byteSize();
    const string_index = try Index.init(table_data[pos..], check);
    pos += string_index.byteSize();
    const global_subrs_index = try Index.init(table_data[pos..], check);

    // Per spec, a CFF table could list multiple fonts in its Name/Top DICT
    // INDEXes; only the first is used - OpenType fonts never do this.
    if (top_dict_index.count() == 0) return CffError.MissingCharStrings;
    const top_dict: Dict = .{ .table_data = top_dict_index.get(0) };

    if (check) {
        const cstype = top_dict.getInt(.charstring_type) orelse 2;
        if (cstype != 2) return CffError.UnsupportedCharstringType;
    }

    const cs_off_raw = top_dict.getInt(.char_strings) orelse return CffError.MissingCharStrings;
    if (cs_off_raw < 0) return CffError.MissingCharStrings;
    const cs_off: usize = @intCast(cs_off_raw);
    if (check) try bytes.checkOffset(table_data, cs_off);
    const charstrings_index = try Index.init(table_data[cs_off..], check);

    const is_cid = top_dict.has(.ros);

    var local_subrs_index: ?Index = null;
    var fd_array_index: ?Index = null;
    var fdselect_offset: ?usize = null;

    if (is_cid) {
        const fdarray_off_raw = top_dict.getInt(.fd_array) orelse return CffError.MissingFdArray;
        const fdselect_off_raw = top_dict.getInt(.fd_select) orelse return CffError.MissingFdSelect;
        if (fdarray_off_raw < 0) return CffError.MissingFdArray;
        if (fdselect_off_raw < 0) return CffError.MissingFdSelect;

        const fdarray_off: usize = @intCast(fdarray_off_raw);
        const fdselect_off: usize = @intCast(fdselect_off_raw);
        if (check) try bytes.checkOffset(table_data, fdarray_off);
        fd_array_index = try Index.init(table_data[fdarray_off..], check);
        if (check) try bytes.checkOffset(table_data, fdselect_off);
        fdselect_offset = fdselect_off;
    } else {
        local_subrs_index = try resolvePrivateSubrs(table_data, top_dict, check);
    }

    return .{
        .table_data = table_data,
        .charstrings_index = charstrings_index,
        .global_subrs_index = global_subrs_index,
        .is_cid = is_cid,
        .local_subrs_index = local_subrs_index,
        .fd_array_index = fd_array_index,
        .fdselect_offset = fdselect_offset,
    };
}

pub inline fn majorVersion(self: CffTable) u8 {
    return self.table_data[@offsetOf(StreamLayout, "major_version")];
}

pub inline fn minorVersion(self: CffTable) u8 {
    return self.table_data[@offsetOf(StreamLayout, "minor_version")];
}

pub inline fn numGlyphs(self: CffTable) u16 {
    return self.charstrings_index.count();
}

/// One glyph's raw, unparsed Type 2 charstring bytes. Bounds are the
/// caller's responsibility (glyph < numGlyphs()), matching this
/// module's other raw index accessors.
pub fn getCharString(self: CffTable, glyph: GlyphIndex) []const u8 {
    return self.charstrings_index.get(@intCast(@backingInt(glyph)));
}

pub fn getGlobalSubr(self: CffTable, index: usize) []const u8 {
    return self.global_subrs_index.get(index);
}

/// Decodes FDSelect's mapping of `glyph` to a Font DICT index. Formats
/// 0 (a direct per-glyph array) and 3 (sorted ranges) are supported;
/// format 4 is CFF2-only. Only meaningful when is_cid is true.
pub fn getFdIndexForGlyph(self: CffTable, glyph: GlyphIndex) u8 {
    const off = self.fdselect_offset orelse return 0;
    const data = self.table_data[off..];
    const gid: usize = @intCast(@backingInt(glyph));

    return switch (data[0]) {
        0 => data[1 + gid],
        3 => blk: {
            const n_ranges = bytes.getU16BE(data[1..]);
            var lo: usize = 0;
            var hi: usize = n_ranges;
            while (lo < hi) {
                const mid = lo + (hi - lo) / 2;
                const first = bytes.getU16BE(data[3 + mid * 3 ..]);
                if (first <= gid) lo = mid + 1 else hi = mid;
            }
            break :blk data[3 + (lo - 1) * 3 + 2];
        },
        else => 0,
    };
}

/// The local Subrs INDEX applicable to `glyph`: the shared one for a
/// non-CID font, or the CID font's own FD's one, resolved fresh per call
/// since it genuinely varies per glyph. Trusts the font (no re-validation)
/// since init() already validated once.
pub fn getLocalSubrIndexForGlyph(self: CffTable, glyph: GlyphIndex) ?Index {
    const fd_array = self.fd_array_index orelse return self.local_subrs_index;
    const fd = self.getFdIndexForGlyph(glyph);
    const font_dict: Dict = .{ .table_data = fd_array.get(fd) };
    return resolvePrivateSubrs(self.table_data, font_dict, false) catch null;
}

const zigstd = @import("../../../../../zigstd.zig");

test "CFF StreamLayout matches the actual spec header (4 bytes, no padding)" {
    try zigstd.testExpectEqual(@as(usize, 4), @sizeOf(StreamLayout));
    try zigstd.testExpectEqual(@as(usize, 1), @offsetOf(StreamLayout, "minor_version"));
    try zigstd.testExpectEqual(@as(usize, 2), @offsetOf(StreamLayout, "header_size"));
    try zigstd.testExpectEqual(@as(usize, 3), @offsetOf(StreamLayout, "off_size"));
}

test "CffTable.getFdIndexForGlyph decodes FDSelect formats 0 and 3" {
    // Format 0: direct per-glyph array.
    const fmt0 = [_]u8{ 0, 0, 2, 1, 0 }; // glyph0->FD0, glyph1->FD2, glyph2->FD1, glyph3->FD0
    var t0: CffTable = undefined;
    t0.table_data = &fmt0;
    t0.fdselect_offset = 0;
    try zigstd.testExpectEqual(@as(u8, 0), t0.getFdIndexForGlyph(@fromBackingInt(@intCast(0))));
    try zigstd.testExpectEqual(@as(u8, 2), t0.getFdIndexForGlyph(@fromBackingInt(@intCast(1))));
    try zigstd.testExpectEqual(@as(u8, 1), t0.getFdIndexForGlyph(@fromBackingInt(@intCast(2))));

    // Format 3: ranges [0,5)->FD0, [5,10)->FD1, sentinel at 10.
    var fmt3: [1 + 2 + 2 * 3 + 2]u8 = undefined;
    fmt3[0] = 3;
    bytes.setU16BE(fmt3[1..], 2); // nRanges
    bytes.setU16BE(fmt3[3..], 0); // range0.first
    fmt3[5] = 0; // range0.fd
    bytes.setU16BE(fmt3[6..], 5); // range1.first
    fmt3[8] = 1; // range1.fd
    bytes.setU16BE(fmt3[9..], 10); // sentinel

    var t3: CffTable = undefined;
    t3.table_data = &fmt3;
    t3.fdselect_offset = 0;
    try zigstd.testExpectEqual(@as(u8, 0), t3.getFdIndexForGlyph(@fromBackingInt(@intCast(0))));
    try zigstd.testExpectEqual(@as(u8, 0), t3.getFdIndexForGlyph(@fromBackingInt(@intCast(4))));
    try zigstd.testExpectEqual(@as(u8, 1), t3.getFdIndexForGlyph(@fromBackingInt(@intCast(5))));
    try zigstd.testExpectEqual(@as(u8, 1), t3.getFdIndexForGlyph(@fromBackingInt(@intCast(9))));
}

// Builds a minimal, non-CID CFF table: a Top DICT with CharStrings+Private,
// a Private DICT with Subrs, a 1-entry local Subrs INDEX, and a 2-glyph
// CharStrings INDEX. DICT integers always use the fixed 5-byte form so
// sizes don't depend on the (forward-referencing) offset values.
fn buildNonCidCff(alloc: zigstd.Allocator) ![]u8 {
    const name_index = try Index.buildIndex(alloc, &.{"Test"});
    defer alloc.free(name_index);
    const string_index = try Index.buildIndex(alloc, &.{});
    defer alloc.free(string_index);
    const global_subrs_index = try Index.buildIndex(alloc, &.{});
    defer alloc.free(global_subrs_index);
    const charstrings_index = try Index.buildIndex(alloc, &.{ &.{0x0E}, &.{ 0x01, 0x02, 0x0E } });
    defer alloc.free(charstrings_index);
    const local_subrs_index = try Index.buildIndex(alloc, &.{&.{ 0xAA, 0xBB, 0xCC }});
    defer alloc.free(local_subrs_index);

    const header_size = 4;
    const private_size = 6; // appendDictInt (5) + Subrs operator byte (1)

    // Top DICT byte length is fixed regardless of operand values (always
    // uses the 5-byte int form), so its own size can be computed up front.
    const top_dict_len = 5 + 5 + 1 + 5 + 1; // Private(size,offset) + CharStrings

    const top_dict_index_len = blk: {
        const dummy = try alloc.alloc(u8, top_dict_len);
        defer alloc.free(dummy);
        @memset(dummy, 0);
        const idx = try Index.buildIndex(alloc, &.{dummy});
        defer alloc.free(idx);
        break :blk idx.len;
    };

    const name_off = header_size;
    const top_dict_off = name_off + name_index.len;
    const string_off = top_dict_off + top_dict_index_len;
    const gsubrs_off = string_off + string_index.len;
    const priv_off = gsubrs_off + global_subrs_index.len;
    const local_subrs_off = priv_off + private_size;
    const charstrings_off = local_subrs_off + local_subrs_index.len;
    const total = charstrings_off + charstrings_index.len;

    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);
    buf[0] = 1; // majorVersion
    buf[1] = 0; // minorVersion
    buf[2] = header_size;
    buf[3] = 4; // offSize (unused by our reader)

    @memcpy(buf[name_off..][0..name_index.len], name_index);

    var top_dict_bytes = try alloc.alloc(u8, top_dict_len);
    defer alloc.free(top_dict_bytes);
    var p: usize = 0;
    p = Dict.appendDictInt(top_dict_bytes, p, @intCast(private_size));
    p = Dict.appendDictInt(top_dict_bytes, p, @intCast(priv_off));
    top_dict_bytes[p] = 18; // Private
    p += 1;
    p = Dict.appendDictInt(top_dict_bytes, p, @intCast(charstrings_off));
    top_dict_bytes[p] = 17; // CharStrings
    p += 1;

    const top_dict_index = try Index.buildIndex(alloc, &.{top_dict_bytes});
    defer alloc.free(top_dict_index);
    @memcpy(buf[top_dict_off..][0..top_dict_index.len], top_dict_index);

    @memcpy(buf[string_off..][0..string_index.len], string_index);
    @memcpy(buf[gsubrs_off..][0..global_subrs_index.len], global_subrs_index);

    var priv_bytes: [private_size]u8 = undefined;
    const subrs_rel = local_subrs_off - priv_off;
    _ = Dict.appendDictInt(&priv_bytes, 0, @intCast(subrs_rel));
    priv_bytes[5] = 19; // Subrs
    @memcpy(buf[priv_off..][0..private_size], &priv_bytes);

    @memcpy(buf[local_subrs_off..][0..local_subrs_index.len], local_subrs_index);
    @memcpy(buf[charstrings_off..][0..charstrings_index.len], charstrings_index);

    return buf;
}

test "CffTable resolves CharStrings and (via Private -> Subrs) local Subrs for a non-CID font" {
    const alloc = zigstd.testAllocator;
    const buf = try buildNonCidCff(alloc);
    defer alloc.free(buf);

    const cff = try CffTable.init(buf, true);
    try zigstd.testExpectEqual(@as(u8, 1), cff.majorVersion());
    try zigstd.testExpectEqual(false, cff.is_cid);
    try zigstd.testExpectEqual(@as(u16, 2), cff.numGlyphs());

    try zigstd.testExpectEqual(@as(usize, 1), cff.getCharString(@fromBackingInt(@intCast(0))).len);
    try zigstd.testExpectEqual(@as(usize, 3), cff.getCharString(@fromBackingInt(@intCast(1))).len);

    const subrs = cff.getLocalSubrIndexForGlyph(@fromBackingInt(@intCast(0))).?;
    try zigstd.testExpectEqual(@as(u16, 1), subrs.count());
    try zigstd.testExpectEqual(true, @import("std").mem.eql(u8, &.{ 0xAA, 0xBB, 0xCC }, subrs.get(0)));
}

// Builds a minimal CID-keyed CFF table: ROS+FDArray+FDSelect in the Top
// DICT, FDSelect format 0 routing glyph 0/1 to FD 0/1, and a 2-entry Font
// DICT INDEX each with its own Private DICT and local Subrs INDEX.
fn buildCidCff(alloc: zigstd.Allocator) ![]u8 {
    const name_index = try Index.buildIndex(alloc, &.{"Test"});
    defer alloc.free(name_index);
    const string_index = try Index.buildIndex(alloc, &.{});
    defer alloc.free(string_index);
    const global_subrs_index = try Index.buildIndex(alloc, &.{});
    defer alloc.free(global_subrs_index);
    const charstrings_index = try Index.buildIndex(alloc, &.{ &.{0x0E}, &.{ 0x01, 0x02, 0x0E } });
    defer alloc.free(charstrings_index);
    const fd0_subrs_index = try Index.buildIndex(alloc, &.{&.{0x11}});
    defer alloc.free(fd0_subrs_index);
    const fd1_subrs_index = try Index.buildIndex(alloc, &.{&.{ 0x22, 0x33 }});
    defer alloc.free(fd1_subrs_index);

    const header_size = 4;
    const private_size = 6; // appendDictInt (5) + Subrs operator byte (1)
    const font_dict_len = 5 + 5 + 1; // Private(size,offset)
    const top_dict_len = (5 + 2) + (5 + 2) + (5 + 2) + (5 + 1); // ROS + FDArray + FDSelect + CharStrings

    const top_dict_index_len = blk: {
        const dummy = try alloc.alloc(u8, top_dict_len);
        defer alloc.free(dummy);
        @memset(dummy, 0);
        const idx = try Index.buildIndex(alloc, &.{dummy});
        defer alloc.free(idx);
        break :blk idx.len;
    };
    const font_dict_index_len = blk: {
        const dummy = try alloc.alloc(u8, font_dict_len);
        defer alloc.free(dummy);
        @memset(dummy, 0);
        const idx = try Index.buildIndex(alloc, &.{ dummy, dummy });
        defer alloc.free(idx);
        break :blk idx.len;
    };

    const name_off = header_size;
    const top_dict_off = name_off + name_index.len;
    const string_off = top_dict_off + top_dict_index_len;
    const gsubrs_off = string_off + string_index.len;
    const fdselect_off = gsubrs_off + global_subrs_index.len;
    const fdselect_len = 3; // format(1) + 2 glyphs
    const fd0_priv_off = fdselect_off + fdselect_len;
    const fd0_subrs_off = fd0_priv_off + private_size;
    const fd1_priv_off = fd0_subrs_off + fd0_subrs_index.len;
    const fd1_subrs_off = fd1_priv_off + private_size;
    const fdarray_off = fd1_subrs_off + fd1_subrs_index.len;
    const charstrings_off = fdarray_off + font_dict_index_len;
    const total = charstrings_off + charstrings_index.len;

    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);
    buf[0] = 1; // majorVersion
    buf[1] = 0; // minorVersion
    buf[2] = header_size;
    buf[3] = 4; // offSize (unused by our reader)

    @memcpy(buf[name_off..][0..name_index.len], name_index);

    var top_dict_bytes = try alloc.alloc(u8, top_dict_len);
    defer alloc.free(top_dict_bytes);
    var p: usize = 0;
    p = Dict.appendDictInt(top_dict_bytes, p, 0); // ROS operand (value unused - only presence matters)
    top_dict_bytes[p] = 12;
    top_dict_bytes[p + 1] = 30; // ROS
    p += 2;
    p = Dict.appendDictInt(top_dict_bytes, p, @intCast(fdarray_off));
    top_dict_bytes[p] = 12;
    top_dict_bytes[p + 1] = 36; // FDArray
    p += 2;
    p = Dict.appendDictInt(top_dict_bytes, p, @intCast(fdselect_off));
    top_dict_bytes[p] = 12;
    top_dict_bytes[p + 1] = 37; // FDSelect
    p += 2;
    p = Dict.appendDictInt(top_dict_bytes, p, @intCast(charstrings_off));
    top_dict_bytes[p] = 17; // CharStrings
    p += 1;

    const top_dict_index = try Index.buildIndex(alloc, &.{top_dict_bytes});
    defer alloc.free(top_dict_index);
    @memcpy(buf[top_dict_off..][0..top_dict_index.len], top_dict_index);

    @memcpy(buf[string_off..][0..string_index.len], string_index);
    @memcpy(buf[gsubrs_off..][0..global_subrs_index.len], global_subrs_index);

    buf[fdselect_off] = 0; // format 0
    buf[fdselect_off + 1] = 0; // glyph 0 -> FD 0
    buf[fdselect_off + 2] = 1; // glyph 1 -> FD 1

    var fd0_priv: [private_size]u8 = undefined;
    _ = Dict.appendDictInt(&fd0_priv, 0, @intCast(fd0_subrs_off - fd0_priv_off));
    fd0_priv[5] = 19; // Subrs
    @memcpy(buf[fd0_priv_off..][0..private_size], &fd0_priv);
    @memcpy(buf[fd0_subrs_off..][0..fd0_subrs_index.len], fd0_subrs_index);

    var fd1_priv: [private_size]u8 = undefined;
    _ = Dict.appendDictInt(&fd1_priv, 0, @intCast(fd1_subrs_off - fd1_priv_off));
    fd1_priv[5] = 19; // Subrs
    @memcpy(buf[fd1_priv_off..][0..private_size], &fd1_priv);
    @memcpy(buf[fd1_subrs_off..][0..fd1_subrs_index.len], fd1_subrs_index);

    var fd0_dict_bytes: [font_dict_len]u8 = undefined;
    var fp: usize = 0;
    fp = Dict.appendDictInt(&fd0_dict_bytes, fp, @intCast(private_size));
    fp = Dict.appendDictInt(&fd0_dict_bytes, fp, @intCast(fd0_priv_off));
    fd0_dict_bytes[fp] = 18; // Private

    var fd1_dict_bytes: [font_dict_len]u8 = undefined;
    fp = 0;
    fp = Dict.appendDictInt(&fd1_dict_bytes, fp, @intCast(private_size));
    fp = Dict.appendDictInt(&fd1_dict_bytes, fp, @intCast(fd1_priv_off));
    fd1_dict_bytes[fp] = 18; // Private

    const fd_array_index = try Index.buildIndex(alloc, &.{ &fd0_dict_bytes, &fd1_dict_bytes });
    defer alloc.free(fd_array_index);
    @memcpy(buf[fdarray_off..][0..fd_array_index.len], fd_array_index);

    @memcpy(buf[charstrings_off..][0..charstrings_index.len], charstrings_index);

    return buf;
}

test "CffTable resolves per-glyph Font DICT/Private/Subrs for a CID-keyed font" {
    const alloc = zigstd.testAllocator;
    const buf = try buildCidCff(alloc);
    defer alloc.free(buf);

    const cff = try CffTable.init(buf, true);
    try zigstd.testExpectEqual(true, cff.is_cid);
    try zigstd.testExpectEqual(@as(u16, 2), cff.numGlyphs());

    const glyph0: GlyphIndex = @fromBackingInt(@intCast(0));
    const glyph1: GlyphIndex = @fromBackingInt(@intCast(1));

    try zigstd.testExpectEqual(@as(u8, 0), cff.getFdIndexForGlyph(glyph0));
    try zigstd.testExpectEqual(@as(u8, 1), cff.getFdIndexForGlyph(glyph1));

    const subrs0 = cff.getLocalSubrIndexForGlyph(glyph0).?;
    try zigstd.testExpectEqual(true, @import("std").mem.eql(u8, &.{0x11}, subrs0.get(0)));

    const subrs1 = cff.getLocalSubrIndexForGlyph(glyph1).?;
    try zigstd.testExpectEqual(true, @import("std").mem.eql(u8, &.{ 0x22, 0x33 }, subrs1.get(0)));
}
