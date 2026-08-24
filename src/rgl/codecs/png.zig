const std = @import("std");

const RGBA = @import("../color.zig").RGBA;
const Image = @import("../Image.zig");

fn appendU32BE(alloc: std.mem.Allocator, list: *std.ArrayList(u8), v: u32) !void {
    try list.appendSlice(alloc, &std.mem.toBytes(std.mem.nativeToBig(u32, v)));
}

fn crc32(data: []const u8) u32 {
    var h: std.hash.Crc32 = .init();
    h.update(data);
    return h.final();
}

fn writeChunk(alloc: std.mem.Allocator, list: *std.ArrayList(u8), chunk_type: *const [4]u8, data: []const u8) !void {
    try appendU32BE(alloc, list, @intCast(data.len));
    const start = list.items.len;
    try list.appendSlice(alloc, chunk_type);
    try list.appendSlice(alloc, data);
    try appendU32BE(alloc, list, crc32(list.items[start..]));
}

/// encodePNG is a minimal, non-interlaced, filter-type-0 (None) RGBA8 PNG
/// encoder - just enough to round-trip genicon's own output through
/// icon.png/icon.ico. See https://www.w3.org/TR/png/ for the chunk layout.
pub fn encode(alloc: std.mem.Allocator, img: *const Image) ![]u8 {
    var out: std.ArrayList(u8) = .empty;
    defer out.deinit(alloc);
    try out.appendSlice(alloc, &[_]u8{ 0x89, 'P', 'N', 'G', '\r', '\n', 0x1A, '\n' });

    var ihdr: [13]u8 = undefined;
    std.mem.writeInt(u32, ihdr[0..4], img.width, .big);
    std.mem.writeInt(u32, ihdr[4..8], img.height, .big);
    ihdr[8] = 8; // bit depth
    ihdr[9] = 6; // color type: truecolor + alpha (RGBA)
    ihdr[10] = 0; // compression method: deflate
    ihdr[11] = 0; // filter method
    ihdr[12] = 0; // interlace method: none
    try writeChunk(alloc, &out, "IHDR", &ihdr);

    // Raw (pre-compression) scanlines: filter-type byte (0/None) + RGBA row bytes.
    const w: usize = @intCast(img.width);
    const h: usize = @intCast(img.height);
    var raw: std.ArrayList(u8) = .empty;
    defer raw.deinit(alloc);
    var y: usize = 0;
    while (y < h) : (y += 1) {
        try raw.append(alloc, 0); // filter type: None
        // RGBA is a packed struct of four u8 fields, so each row's pixels are
        // already byte-for-byte identical to PNG's expected R,G,B,A order.
        try raw.appendSlice(alloc, std.mem.sliceAsBytes(img.data[y * w .. (y + 1) * w]));
    }

    var idat: std.Io.Writer.Allocating = try .initCapacity(alloc, 4096);
    defer idat.deinit();
    var window: [std.compress.flate.max_window_len]u8 = undefined;
    var compressor = try std.compress.flate.Compress.init(&idat.writer, &window, .zlib, .default);
    try compressor.writer.writeAll(raw.items);
    try compressor.finish();
    try writeChunk(alloc, &out, "IDAT", idat.written());

    try writeChunk(alloc, &out, "IEND", &.{});
    return out.toOwnedSlice(alloc);
}

/// decodePNG is a minimal decoder covering what real-world PNG exporters
/// actually produce: 8-bit depth, non-interlaced, RGB or RGBA color type,
/// any per-scanline filter (None/Sub/Up/Average/Paeth), any chunk ordering
/// with unknown ancillary chunks (pHYs, tEXt, gAMA, iCCP, ...) skipped.
/// Palette, grayscale, 16-bit, and interlaced PNGs are out of scope.
pub fn decode(alloc: std.mem.Allocator, bytes: []const u8) !Image {
    const sig = [_]u8{ 0x89, 'P', 'N', 'G', '\r', '\n', 0x1A, '\n' };
    if (bytes.len < sig.len or !std.mem.eql(u8, bytes[0..sig.len], &sig)) {
        return error.InvalidPngSignature;
    }

    var width: u32 = 0;
    var height: u32 = 0;
    var bpp: usize = 0; // bytes per pixel: 3 (RGB) or 4 (RGBA), at 8-bit depth
    var idat: std.ArrayList(u8) = .empty;
    defer idat.deinit(alloc);

    var pos: usize = sig.len;
    while (pos + 8 <= bytes.len) {
        const len: usize = std.mem.readInt(u32, bytes[pos..][0..4], .big);
        const chunk_type = bytes[pos + 4 ..][0..4];
        const data_start = pos + 8;
        if (data_start + len + 4 > bytes.len) return error.InvalidPngChunk;
        const data = bytes[data_start..][0..len];

        if (std.mem.eql(u8, chunk_type, "IHDR")) {
            if (len != 13) return error.InvalidPngChunk;
            width = std.mem.readInt(u32, data[0..4], .big);
            height = std.mem.readInt(u32, data[4..8], .big);
            const bit_depth = data[8];
            const color_type = data[9];
            const interlace = data[12];
            if (bit_depth != 8 or interlace != 0) return error.UnsupportedPngFormat;
            bpp = switch (color_type) {
                2 => 3, // RGB
                6 => 4, // RGBA
                else => return error.UnsupportedPngFormat,
            };
        } else if (std.mem.eql(u8, chunk_type, "IDAT")) {
            try idat.appendSlice(alloc, data);
        } else if (std.mem.eql(u8, chunk_type, "IEND")) {
            break;
        }
        // else: skip unknown/ancillary chunk

        pos = data_start + len + 4;
    }
    if (width == 0 or height == 0 or bpp == 0) return error.MissingIhdr;
    if (idat.items.len == 0) return error.MissingIdat;

    const w: usize = @intCast(width);
    const h: usize = @intCast(height);
    const stride = w * bpp;

    var idat_reader: std.Io.Reader = .fixed(idat.items);
    var inflate_buf: [std.compress.flate.max_window_len]u8 = undefined;
    var decompress: std.compress.flate.Decompress = .init(&idat_reader, .zlib, &inflate_buf);

    var raw: std.Io.Writer.Allocating = try .initCapacity(alloc, h * (stride + 1));
    defer raw.deinit();
    _ = try decompress.reader.streamRemaining(&raw.writer);
    const filtered = raw.written();
    if (filtered.len != h * (stride + 1)) return error.TruncatedPngData;

    var img = try Image.init(alloc, width, height, .{ .r = 0, .g = 0, .b = 0, .a = 0 });

    // Un-filter each scanline (PNG filters are defined relative to the
    // pixel to the left, the pixel above, and the pixel above-left - all
    // in the already-unfiltered image, not the raw compressed bytes).
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

    return img;
}

fn paeth(a: u8, b: u8, c: u8) u8 {
    const p = @as(i32, a) + @as(i32, b) - @as(i32, c);
    const pa = @abs(p - @as(i32, a));
    const pb = @abs(p - @as(i32, b));
    const pc = @abs(p - @as(i32, c));
    if (pa <= pb and pa <= pc) return a;
    if (pb <= pc) return b;
    return c;
}

test "encode/decode round-trips a small RGBA image" {
    const alloc = std.testing.allocator;
    var img = try Image.init(alloc, 3, 2, null);
    defer img.destroy(alloc);
    img.setPixel(0, 0, .{ .r = 255, .g = 0, .b = 0, .a = 255 });
    img.setPixel(1, 0, .{ .r = 0, .g = 255, .b = 0, .a = 128 });
    img.setPixel(2, 0, .{ .r = 0, .g = 0, .b = 255, .a = 0 });
    img.setPixel(0, 1, .{ .r = 10, .g = 20, .b = 30, .a = 40 });
    img.setPixel(1, 1, .{ .r = 255, .g = 255, .b = 255, .a = 255 });
    img.setPixel(2, 1, .{ .r = 0, .g = 0, .b = 0, .a = 0 });

    const bytes = try encode(alloc, &img);
    defer alloc.free(bytes);

    var decoded = try decode(alloc, bytes);
    defer decoded.destroy(alloc);

    try std.testing.expectEqual(img.width, decoded.width);
    try std.testing.expectEqual(img.height, decoded.height);
    for (img.data, decoded.data) |a, b| try std.testing.expect(a.equals(b));
}

test "decode rejects a non-PNG signature" {
    const alloc = std.testing.allocator;
    try std.testing.expectError(error.InvalidPngSignature, decode(alloc, "not a png"));
}
