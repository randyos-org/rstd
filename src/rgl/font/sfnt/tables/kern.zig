const types = @import("../types.zig");
const bytes = @import("../../../../bytes.zig");

const GlyphIndex = types.GlyphIndex;
const zigstd = @import("../../../../zigstd.zig");

/// The 'kern' table header. This is the OpenType/Windows header (16-bit
/// version); Apple's own extended kern format (a 32-bit version field, plus
/// subtable formats 1 and 3) isn't part of the OpenType spec and isn't
/// supported here. Per spec, CFF-outline fonts don't use 'kern' at all; they
/// carry kerning in GPOS instead.
const StreamLayout = extern struct {
    version: u16,
    n_tables: u16,
};

/// One kerning subtable's header, shared by every subtable format.
const SubtableHeader = extern struct {
    version: u16,
    length: u16,
    coverage: u16,
};

/// The `coverage` field's format sub-field (bits 8-15) selects which
/// OpenType-defined subtable layout follows. Only format 0 (a sorted pair
/// list, "the only subtable format supported by Windows" per spec) is
/// implemented; format 2 (a class-pair 2D array) is rare in practice and
/// left unimplemented; KernSubtable.getKernValue() returns null for it.
pub const Coverage = packed struct(u16) {
    horizontal: bool,
    minimum: bool,
    cross_stream: bool,
    override: bool,
    _reserved: u4 = 0,
    format: u8,
};

/// Legacy pair-kerning table, superseded by GPOS's pair-adjustment positioning
/// (LookupType 2) but still checked as a fallback by real-world text renderers.
pub const KernTable = struct {
    const Self = @This();

    table_data: []const u8,

    pub fn init(table_data: []const u8, check: bool) !Self {
        if (check) try bytes.checkSize(table_data, 0, StreamLayout);
        return .{ .table_data = table_data };
    }

    pub inline fn version(self: Self) u16 {
        return types.field(StreamLayout, "version", self.table_data);
    }

    pub inline fn numTables(self: Self) u16 {
        return types.field(StreamLayout, "n_tables", self.table_data);
    }

    /// Byte offset (from the start of table_data) of subtable `index`.
    /// Subtables are variable-length (each carries its own `length`), so
    /// finding one past index 0 means walking every earlier subtable.
    pub fn getSubtableOffset(self: Self, index: usize) usize {
        var off: usize = @sizeOf(StreamLayout);
        var i: usize = 0;
        while (i < index) : (i += 1) {
            const length = bytes.getU16BE(self.table_data[off + @offsetOf(SubtableHeader, "length") ..]);
            off += length;
        }
        return off;
    }

    pub fn getSubtable(self: Self, index: usize) KernSubtable {
        return .{ .table_data = self.table_data[self.getSubtableOffset(index)..] };
    }

    /// Scans every subtable for an explicit kerning value between `left`
    /// and `right`, returning the first one found. Per spec, kerning-value
    /// subtables are additive and their order doesn't matter for a single
    /// lookup like this; minimum-value subtables (Coverage.minimum) are
    /// skipped since they clamp rather than contribute a kerning amount.
    pub fn getKernValue(self: Self, left: GlyphIndex, right: GlyphIndex) ?types.FWord {
        var i: usize = 0;
        while (i < self.numTables()) : (i += 1) {
            const subtable = self.getSubtable(i);
            if (subtable.coverage().minimum) continue;
            if (subtable.getKernValue(left, right)) |v| return v;
        }
        return null;
    }
};

/// One kerning subtable. table_data starts at this subtable's own header
/// (its `length` field bounds its extent, but table_data itself may run
/// further, into whatever follows in the parent 'kern' table).
pub const KernSubtable = struct {
    const Self = @This();

    table_data: []const u8,

    pub inline fn version(self: Self) u16 {
        return types.field(SubtableHeader, "version", self.table_data);
    }

    pub inline fn length(self: Self) u16 {
        return types.field(SubtableHeader, "length", self.table_data);
    }

    pub inline fn coverage(self: Self) Coverage {
        return @bitCast(types.field(SubtableHeader, "coverage", self.table_data));
    }

    const header_end = @sizeOf(SubtableHeader);

    /// Schema for format 0's own header (follows SubtableHeader), used only
    /// to compute field offsets/sizes at comptime.
    const Format0Header = extern struct {
        n_pairs: u16,
        search_range: u16,
        entry_selector: u16,
        range_shift: u16,
    };

    const pairs_start = header_end + @sizeOf(Format0Header);

    /// One entry in a format 0 subtable's sorted pair list. Its extern
    /// layout doubles as the wire schema for @sizeOf.
    pub const KernPair = extern struct {
        left: GlyphIndex,
        right: GlyphIndex,
        value: types.FWord,
    };

    pub inline fn numPairs(self: Self) u16 {
        return bytes.getU16BE(self.table_data[header_end + @offsetOf(Format0Header, "n_pairs") ..]);
    }

    pub fn getPair(self: Self, index: usize) KernPair {
        return bytes.typeBEToNative(KernPair, self.table_data[pairs_start + index * @sizeOf(KernPair) ..]);
    }

    /// Binary search for an exact (left, right) glyph pair, per format 0's
    /// combined-32-bit-key sort order (left in the high word). Returns null
    /// if this isn't a format 0 subtable, or the pair has no explicit entry
    /// (i.e. no kerning between them).
    pub fn getKernValue(self: Self, left: GlyphIndex, right: GlyphIndex) ?types.FWord {
        if (self.coverage().format != 0) return null;

        const target = (@as(u32, @backingInt(left)) << 16) | @backingInt(right);
        var lo: usize = 0;
        var hi: usize = self.numPairs();
        while (lo < hi) {
            const mid = lo + (hi - lo) / 2;
            const pair = self.getPair(mid);
            const key = (@as(u32, @backingInt(pair.left)) << 16) | @backingInt(pair.right);
            if (target < key) {
                hi = mid;
            } else if (target > key) {
                lo = mid + 1;
            } else {
                return pair.value;
            }
        }
        return null;
    }
};

test "kern StreamLayout/SubtableHeader/Format0Header match the OpenType stream layout" {
    try zigstd.testExpectEqual(@as(usize, 4), @sizeOf(StreamLayout));
    try zigstd.testExpectEqual(@as(usize, 2), @offsetOf(StreamLayout, "n_tables"));
    try zigstd.testExpectEqual(@as(usize, 6), @sizeOf(SubtableHeader));
    try zigstd.testExpectEqual(@as(usize, 4), @offsetOf(SubtableHeader, "coverage"));
    try zigstd.testExpectEqual(@as(usize, 6), @sizeOf(KernSubtable.KernPair));
}

fn buildKern(alloc: zigstd.Allocator, pairs: []const KernSubtable.KernPair) ![]u8 {
    const subtable_header_size = @sizeOf(SubtableHeader);
    const format0_header_size = @sizeOf(KernSubtable.Format0Header);
    const subtable_size = subtable_header_size + format0_header_size + pairs.len * @sizeOf(KernSubtable.KernPair);
    const total = @sizeOf(StreamLayout) + subtable_size;

    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);

    bytes.setU16BE(buf, 0); // version
    bytes.setU16BE(buf[@offsetOf(StreamLayout, "n_tables")..], 1);

    const sub = buf[@sizeOf(StreamLayout)..];
    bytes.setU16BE(sub[@offsetOf(SubtableHeader, "version")..], 0);
    bytes.setU16BE(sub[@offsetOf(SubtableHeader, "length")..], @intCast(subtable_size));
    const coverage: Coverage = .{ .horizontal = true, .minimum = false, .cross_stream = false, .override = false, .format = 0 };
    bytes.setU16BE(sub[@offsetOf(SubtableHeader, "coverage")..], @bitCast(coverage));

    const f0 = sub[subtable_header_size..];
    bytes.setU16BE(f0[@offsetOf(KernSubtable.Format0Header, "n_pairs")..], @intCast(pairs.len));

    const pairs_off = subtable_header_size + format0_header_size;
    var off: usize = pairs_off;
    for (pairs) |p| {
        bytes.setU16BE(sub[off..], @backingInt(p.left));
        bytes.setU16BE(sub[off + @sizeOf(u16) ..], @backingInt(p.right));
        bytes.setI16BE(sub[off + 2 * @sizeOf(u16) ..], p.value);
        off += @sizeOf(KernSubtable.KernPair);
    }

    return buf;
}

test "KernTable.getKernValue resolves an exact pair via format 0" {
    const alloc = zigstd.testAllocator;
    const buf = try buildKern(alloc, &.{
        .{ .left = @fromBackingInt(@intCast(5)), .right = @fromBackingInt(@intCast(10)), .value = -20 },
        .{ .left = @fromBackingInt(@intCast(5)), .right = @fromBackingInt(@intCast(20)), .value = -30 },
        .{ .left = @fromBackingInt(@intCast(7)), .right = @fromBackingInt(@intCast(1)), .value = 15 },
    });
    defer alloc.free(buf);

    const kern = try KernTable.init(buf, true);

    try zigstd.testExpectEqual(@as(?types.FWord, -20), kern.getKernValue(@fromBackingInt(@intCast(5)), @fromBackingInt(@intCast(10))));
    try zigstd.testExpectEqual(@as(?types.FWord, -30), kern.getKernValue(@fromBackingInt(@intCast(5)), @fromBackingInt(@intCast(20))));
    try zigstd.testExpectEqual(@as(?types.FWord, 15), kern.getKernValue(@fromBackingInt(@intCast(7)), @fromBackingInt(@intCast(1))));
}

test "KernTable.getKernValue returns null for a pair with no entry" {
    const alloc = zigstd.testAllocator;
    const buf = try buildKern(alloc, &.{
        .{ .left = @fromBackingInt(@intCast(5)), .right = @fromBackingInt(@intCast(10)), .value = -20 },
    });
    defer alloc.free(buf);

    const kern = try KernTable.init(buf, true);

    try zigstd.testExpectEqual(@as(?types.FWord, null), kern.getKernValue(@fromBackingInt(@intCast(5)), @fromBackingInt(@intCast(11))));
    try zigstd.testExpectEqual(@as(?types.FWord, null), kern.getKernValue(@fromBackingInt(@intCast(6)), @fromBackingInt(@intCast(10))));
}
