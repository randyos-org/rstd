const std = @import("std");
const flate = std.compress.flate;

const color = @import("../color.zig");
const RGBA = color.RGBA;
const RgbaImage = @import("../image.zig").RgbaImage;
const bytes = @import("../../bytes.zig");
const zigstd = @import("../../zigstd.zig");
const Allocator = zigstd.Allocator;
const StdIo = zigstd.StdIo;
const equal = zigstd.equal;

const appendIntBE = bytes.appendIntBE;
const appendArrayBytes = bytes.appendArrayBytes;
const getU32BE = bytes.getU32BE;
const setU32BE = bytes.setU32BE;
const ByteArrayList = bytes.ByteArrayList;
const Crc32 = std.hash.Crc32;
// const structAsBytes = bytes.structAsBytes;

const png_signature = [_]u8{ 0x89, 'P', 'N', 'G', '\r', '\n', 0x1A, '\n' };

pub const PngChunkType = enum {
    const Self = @This();

    IHDR,
    IDAT,
    IEND,

    pub fn toInt(id: Self) u32 {
        var array4: [4]u8 = undefined;
        @memcpy(&array4, @tagName(id));
        return getU32BE(&array4);
    }

    pub fn fromInt(v: u32) !Self {
        var array4: [4]u8 = undefined;
        setU32BE(&array4, v);
        return std.meta.stringToEnum(Self, &array4) orelse error.UnsupportedPngChunk;
    }
};

pub const PngChunk = struct {
    const Self = @This();
    const type_start = 4; // 4 bytes length
    const data_start = 8; // 4 bytes length + 4 bytes type
    const meta_size = data_start + 4; // 4 bytes CRC32

    chunk_type: PngChunkType,
    data: []const u8,

    pub fn len(self: *const Self) usize {
        return meta_size + self.data.len;
    }

    pub fn appendToByteArrayList(self: *const Self, alloc: Allocator, list: *ByteArrayList) !void {
        try appendIntBE(alloc, list, u32, self.data.len);
        const start = list.items.len;
        try appendIntBE(alloc, list, u32, self.chunk_type.toInt());
        try list.appendSlice(alloc, self.data);
        try appendIntBE(alloc, list, u32, Crc32.hash(list.items[start..]));
    }

    pub fn fromBytes(in: []const u8) !Self {
        if (in.len < meta_size) return error.InvalidPngChunk;
        const data_len = getU32BE(in);
        const total_len = meta_size + data_len;
        const crc_start = total_len - 4;

        if (in.len < total_len) return error.InvalidPngChunk;
        const chunk_type = try PngChunkType.fromInt(getU32BE(in[type_start..]));
        const data = in[data_start..crc_start];
        const crc = getU32BE(in[crc_start..]);
        if (crc != Crc32.hash(in[type_start..crc_start])) return error.InvalidPngChunk;
        return .{
            .chunk_type = chunk_type,
            .data = data,
        };
    }
};

pub const ColorType = enum(u8) {
    RGB = 2,
    RGBA = 6,
};

pub const IHdr = struct {
    const Self = @This();
    const size = 13;
    width: u32,
    height: u32,
    bit_depth: u8,
    color_type: ColorType,
    compression_method: u8,
    filter_method: u8,
    interlace_method: u8,

    pub fn toBytes(self: *Self, out: *[size]u8) void {
        setU32BE(out, self.width);
        setU32BE(out[4..], self.height);
        out[8] = self.bit_depth;
        out[9] = @backingInt(self.color_type);
        out[10] = self.compression_method;
        out[11] = self.filter_method;
        out[12] = self.interlace_method;
    }

    pub fn fromBytes(in: *const [size]u8) Self {
        return .{
            .width = getU32BE(in),
            .height = getU32BE(in[4..]),
            .bit_depth = in[8],
            .color_type = @fromBackingInt(in[9]),
            .compression_method = in[10],
            .filter_method = in[11],
            .interlace_method = in[12],
        };
    }

    pub fn bytesPerPixel(self: *const Self) usize {
        return switch (self.color_type) {
            .RGB => 3,
            .RGBA => 4,
        };
    }

    pub fn appendToByteArrayList(self: *Self, alloc: Allocator, list: *ByteArrayList) !void {
        var buf: [size]u8 = undefined;
        self.toBytes(&buf);
        try (PngChunk{
            .chunk_type = PngChunkType.IHDR,
            .data = buf[0..],
        }).appendToByteArrayList(alloc, list);
    }
};

/// a minimal, non-interlaced, filter-type-0 (None) RGBA8 PNG
/// encoder. See https://www.w3.org/TR/png/ for the chunk layout.
pub fn encode(alloc: Allocator, img: *const RgbaImage) ![]u8 {
    // Raw (pre-compression) scanlines: filter-type byte (0/None) + RGBA row bytes.
    const w: usize = @intCast(img.width);
    const h: usize = @intCast(img.height);
    var uncompressed: ByteArrayList = .empty;
    defer uncompressed.deinit(alloc);
    var y: usize = 0;
    while (y < h) : (y += 1) {
        try uncompressed.append(alloc, 0); // filter type: None
        // a packed RGBA row is PNG's R,G,B,A bytes as-is, unless `swapped`
        const row = img.data[y * w .. (y + 1) * w];
        if (color.swapped) {
            for (row) |p| try uncompressed.appendSlice(alloc, &p.toArray());
        } else try appendArrayBytes(alloc, &uncompressed, row);
    }
    return encodeScanlines(alloc, img.width, img.height, uncompressed.items);
}

/// An RGBA8 PNG of `scanlines`, each a filter-type byte and the row's bytes.
fn encodeScanlines(alloc: Allocator, width: u32, height: u32, scanlines: []const u8) ![]u8 {
    var out: ByteArrayList = .empty;
    defer out.deinit(alloc);
    try out.appendSlice(alloc, &png_signature);

    var ihdr: IHdr = .{
        .width = width,
        .height = height,
        .bit_depth = 8, // bit depth
        .color_type = ColorType.RGBA,
        .compression_method = 0, // compression method: deflate
        .filter_method = 0, // filter method
        .interlace_method = 0, // interlace method: none
    };
    try ihdr.appendToByteArrayList(alloc, &out);

    var idat: StdIo.Writer.Allocating = try .initCapacity(alloc, 4096);
    defer idat.deinit();
    var window: [flate.max_window_len]u8 = undefined;
    var compressor = try flate.Compress.init(&idat.writer, &window, .zlib, .default);
    try compressor.writer.writeAll(scanlines);
    try compressor.finish();
    try (PngChunk{
        .chunk_type = PngChunkType.IDAT,
        .data = idat.written(),
    }).appendToByteArrayList(alloc, &out);

    try (PngChunk{
        .chunk_type = PngChunkType.IEND,
        .data = &.{},
    }).appendToByteArrayList(alloc, &out);
    return out.toOwnedSlice(alloc);
}

pub const DecodeOptions = struct {
    /// Print a diagnostic when an unrecognized ancillary chunk is skipped.
    /// Off by default: this goes to stderr, which build-time callers (e.g.
    /// genicon) treat as a failure even at exit code 0.
    debug: bool = false,
};

/// a minimal decoder for 8-bit depth, non-interlaced, RGB or RGBA color type,
/// any per-scanline filter (None/Sub/Up/Average/Paeth), any chunk ordering
/// with unrecognized ancillary chunks (pHYs, tEXt, gAMA, iCCP, ...) skipped.
/// Palette, grayscale, 16-bit, and interlaced PNGs are out of scope.
pub fn decode(alloc: Allocator, raw_bytes: []const u8, options: DecodeOptions) !RgbaImage {
    const sig_len = png_signature.len;
    if (raw_bytes.len < sig_len or !equal(u8, raw_bytes[0..sig_len], &png_signature)) {
        return error.InvalidPngSignature;
    }
    var pos: usize = sig_len;

    const ihdr_chunk = try PngChunk.fromBytes(raw_bytes[pos..]);
    if (ihdr_chunk.chunk_type != PngChunkType.IHDR) return error.MissingIhdr;
    if (ihdr_chunk.data.len != IHdr.size) return error.InvalidPngChunk;
    const ihdr: IHdr = .fromBytes(ihdr_chunk.data[0..IHdr.size]);
    if (ihdr.bit_depth != 8 or ihdr.interlace_method != 0) return error.UnsupportedPngFormat;
    const width = ihdr.width;
    const height = ihdr.height;
    const bpp = ihdr.bytesPerPixel();
    pos += ihdr_chunk.len();

    var idat = try collectIdat(alloc, raw_bytes, pos, options);
    defer idat.deinit(alloc);

    if (width == 0 or height == 0 or bpp == 0) return error.MissingIhdr;
    if (idat.items.len == 0) return error.MissingIdat;

    const w: usize = @intCast(width);
    const h: usize = @intCast(height);
    const stride = w * bpp;

    var idat_reader: StdIo.Reader = .fixed(idat.items);
    var inflate_buf: [flate.max_window_len]u8 = undefined;
    var decompress: flate.Decompress = .init(&idat_reader, .zlib, &inflate_buf);

    var uncompressed: StdIo.Writer.Allocating = try .initCapacity(alloc, h * (stride + 1));
    defer uncompressed.deinit();
    _ = try decompress.reader.streamRemaining(&uncompressed.writer);
    const filtered = uncompressed.written();
    if (filtered.len != h * (stride + 1)) return error.TruncatedPngData;

    var img = try RgbaImage.init(alloc, width, height, .blank);
    errdefer img.deinit(alloc);
    try unfilterRows(alloc, img, filtered, w, h, bpp);
    return img;
}

/// Every IDAT chunk's data from `start` to IEND, concatenated (PNG may split the
/// compressed stream). PngChunkType.fromInt rejects any chunk type other than
/// IHDR/IDAT/IEND; ancillary chunks (pHYs, tEXt, gAMA, iCCP, ...) are skipped.
fn collectIdat(alloc: Allocator, raw_bytes: []const u8, start: usize, options: DecodeOptions) !ByteArrayList {
    var pos = start;
    var idat: ByteArrayList = .empty;
    errdefer idat.deinit(alloc);
    while (true) {
        const chunk = PngChunk.fromBytes(raw_bytes[pos..]) catch |err| switch (err) {
            error.UnsupportedPngChunk => {
                const data_len = getU32BE(raw_bytes[pos..]);
                if (options.debug) {
                    const tag = raw_bytes[pos..][PngChunk.type_start..PngChunk.data_start];
                    std.debug.print("png.decode: ignoring unsupported chunk \"{s}\"\n", .{tag});
                }
                pos += PngChunk.meta_size + data_len;
                continue;
            },
            else => return err,
        };
        pos += chunk.len();
        switch (chunk.chunk_type) {
            .IDAT => try idat.appendSlice(alloc, chunk.data),
            .IEND => break,
            .IHDR => return error.InvalidPngChunk,
        }
    }

    return idat;
}

/// Un-filters each scanline into `img`: PNG filters are relative to the pixel
/// to the left, above and above-left, all in the already-unfiltered image.
fn unfilterRows(alloc: Allocator, img: RgbaImage, filtered: []const u8, w: usize, h: usize, bpp: usize) !void {
    const stride = w * bpp;
    const line_a = try alloc.alloc(u8, stride);
    defer alloc.free(line_a);
    @memset(line_a, 0); // the "row above" row 0 is defined as all zero
    const line_b = try alloc.alloc(u8, stride);
    defer alloc.free(line_b);

    // prev_line/cur_line get swapped each iteration below; line_a/line_b keep
    // their original allocations so the defers above free them exactly once.
    var prev_line: []u8 = line_a;
    var cur_line: []u8 = line_b;

    var y: usize = 0;
    while (y < h) : (y += 1) {
        const row_start = y * (stride + 1);
        const filter_type = filtered[row_start];
        const src_row = filtered[row_start + 1 ..][0..stride];

        var x: usize = 0;
        while (x < stride) : (x += 1) {
            const a: u8 = if (x >= bpp) cur_line[x - bpp] else 0;
            const b: u8 = prev_line[x];
            const c: u8 = if (x >= bpp) prev_line[x - bpp] else 0;
            cur_line[x] = switch (filter_type) {
                0 => src_row[x],
                1 => src_row[x] +% a,
                2 => src_row[x] +% b,
                3 => src_row[x] +% @as(u8, @intCast((@as(u16, a) + @as(u16, b)) / 2)),
                4 => src_row[x] +% paeth(a, b, c),
                else => return error.UnsupportedPngFilter,
            };
        }

        var px: usize = 0;
        while (px < w) : (px += 1) {
            const c = cur_line[px * bpp ..][0..bpp];
            img.setPixel(@intCast(px), @intCast(y), .{
                .r = c[0],
                .g = c[1],
                .b = c[2],
                .a = if (bpp == 4) c[3] else 0xFF,
            });
        }

        std.mem.swap([]u8, &prev_line, &cur_line);
    }
}

fn paeth(a: u8, b: u8, c: u8) u8 {
    const ia: i32 = @intCast(a);
    const ib: i32 = @intCast(b);
    const ic: i32 = @intCast(c);
    const p = ia + ib - ic;
    const pa = @abs(p - ia);
    const pb = @abs(p - ib);
    const pc = @abs(p - ic);
    if (pa <= pb and pa <= pc) return a;
    if (pb <= pc) return b;
    return c;
}

test "encode/decode round-trips a small RGBA image" {
    const alloc = zigstd.testAllocator;
    var img = try RgbaImage.init(alloc, 3, 2, null);
    defer img.deinit(alloc);
    img.setPixel(0, 0, .{ .r = 255, .g = 0, .b = 0, .a = 255 });
    img.setPixel(1, 0, .{ .r = 0, .g = 255, .b = 0, .a = 128 });
    img.setPixel(2, 0, .{ .r = 0, .g = 0, .b = 255, .a = 0 });
    img.setPixel(0, 1, .{ .r = 10, .g = 20, .b = 30, .a = 40 });
    img.setPixel(1, 1, .{ .r = 255, .g = 255, .b = 255, .a = 255 });
    img.setPixel(2, 1, .{ .r = 0, .g = 0, .b = 0, .a = 0 });

    const raw_bytes = try encode(alloc, &img);
    defer alloc.free(raw_bytes);

    var decoded = try decode(alloc, raw_bytes, .{});
    defer decoded.deinit(alloc);

    try zigstd.testExpectEqual(img.width, decoded.width);
    try zigstd.testExpectEqual(img.height, decoded.height);
    for (img.data, decoded.data) |a, b| try zigstd.testExpectTrue(a.equals(b));
}

test "decode frees the image when a scanline's filter type is invalid" {
    const alloc = zigstd.testAllocator;
    const raw_bytes = try encodeScanlines(alloc, 1, 1, &.{ 5, 1, 2, 3, 4 });
    defer alloc.free(raw_bytes);
    try zigstd.testExpectError(error.UnsupportedPngFilter, decode(alloc, raw_bytes, .{}));
}

test "decode rejects a non-PNG signature" {
    const alloc = zigstd.testAllocator;
    try zigstd.testExpectError(error.InvalidPngSignature, decode(alloc, "not a png", .{}));
}
