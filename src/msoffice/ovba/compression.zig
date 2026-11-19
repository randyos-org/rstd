//! [MS-OVBA] 2.4.1 compression: the container of every VBA stream.

const std = @import("std");
const zigstd = @import("../../zigstd.zig");
const bytes = @import("../../bytes.zig");
const Allocator = zigstd.Allocator;
const ByteArrayList = bytes.ByteArrayList;

const chunk_size = 4096;

const CopyTokenHelp = struct {
    /// Bits holding the offset; 4..12 for positions within a 4096-byte chunk.
    bit_count: u4,
    max_length: u16,

    /// `difference` is the current position relative to the chunk start.
    fn init(difference: usize) CopyTokenHelp {
        // Smallest integer >= log2(difference), but at least 4.
        const bits: usize = if (difference <= 1) 0 else std.math.log2_int(usize, difference - 1) + 1;
        const bit_count: u4 = @intCast(@max(bits, 4));
        return .{ .bit_count = bit_count, .max_length = (@as(u16, 0xFFFF) >> bit_count) + 3 };
    }

    fn lengthMask(self: CopyTokenHelp) u16 {
        return @as(u16, 0xFFFF) >> self.bit_count;
    }

    /// Shift that positions the offset in the high bits of a copy token.
    fn offsetShift(self: CopyTokenHelp) u4 {
        return @intCast(16 - @as(u5, self.bit_count));
    }
};

/// Compresses `src` into a CompressedContainer. Caller owns the result.
pub fn compress(gpa: Allocator, src: []const u8) ![]u8 {
    var out: ByteArrayList = .empty;
    errdefer out.deinit(gpa);
    try out.append(gpa, 0x01); // SignatureByte
    var start: usize = 0;
    while (start < src.len) {
        const end = @min(start + chunk_size, src.len);
        try compressChunk(gpa, &out, src[start..end]);
        start = end;
    }
    return out.toOwnedSlice(gpa);
}

fn compressChunk(gpa: Allocator, out: *ByteArrayList, chunk: []const u8) !void {
    const header_pos = out.items.len;
    try out.appendSlice(gpa, &.{ 0, 0 });
    var pos: usize = 0;
    while (pos < chunk.len) {
        const flag_pos = out.items.len;
        try out.append(gpa, 0);
        var flags: u8 = 0;
        for (0..8) |bit| {
            if (pos >= chunk.len) break;
            const help: CopyTokenHelp = .init(pos);
            var best_len: usize = 0;
            var best_off: usize = 0;
            var cand = pos;
            while (cand > 0) {
                cand -= 1;
                var len: usize = 0;
                while (pos + len < chunk.len and len < help.max_length and chunk[cand + len] == chunk[pos + len]) len += 1;
                if (len > best_len) {
                    best_len = len;
                    best_off = pos - cand;
                }
            }
            if (best_len >= 3) {
                const token = (@as(u16, @intCast(best_off - 1)) << help.offsetShift()) | @as(u16, @intCast(best_len - 3));
                try bytes.appendIntLE(gpa, out, u16, token);
                flags |= @as(u8, 1) << @intCast(bit);
                pos += best_len;
            } else {
                try out.append(gpa, chunk[pos]);
                pos += 1;
            }
        }
        out.items[flag_pos] = flags;
    }

    const data_len = out.items.len - header_pos - 2;
    if (data_len > chunk_size) {
        // Incompressible: store raw, which must be exactly 4096 bytes.
        if (chunk.len != chunk_size) return error.VbaChunkIncompressible;
        out.shrinkRetainingCapacity(header_pos);
        try bytes.appendIntLE(gpa, out, u16, 0x3000 | (chunk_size - 1));
        try out.appendSlice(gpa, chunk);
        return;
    }
    // Header: chunk size minus 3 (12 bits), signature 0b011, compressed flag.
    const header: u16 = 0x8000 | 0x3000 | @as(u16, @intCast(data_len + 2 - 3));
    bytes.writeIntLE(u16, out.items[header_pos..][0..2], header);
}

/// Decompresses a CompressedContainer. Caller owns the result.
pub fn decompress(gpa: Allocator, src: []const u8) ![]u8 {
    if (src.len == 0 or src[0] != 0x01) return error.VbaBadContainer;
    var out: ByteArrayList = .empty;
    errdefer out.deinit(gpa);
    var pos: usize = 1;
    while (pos < src.len) {
        if (pos + 2 > src.len) return error.VbaTruncated;
        const header = bytes.readIntLE(u16, src[pos..][0..2]);
        if ((header >> 12) & 0x7 != 0x3) return error.VbaBadChunk;
        const end = @min(pos + (header & 0x0FFF) + 3, src.len);
        pos += 2;
        const chunk_start = out.items.len;
        if (header & 0x8000 == 0) {
            if (pos + chunk_size > src.len) return error.VbaTruncated;
            try out.appendSlice(gpa, src[pos..][0..chunk_size]);
            pos += chunk_size;
            continue;
        }
        while (pos < end) {
            const flags = src[pos];
            pos += 1;
            for (0..8) |bit| {
                if (pos >= end) break;
                if ((flags >> @intCast(bit)) & 1 == 0) {
                    try out.append(gpa, src[pos]);
                    pos += 1;
                } else {
                    if (pos + 2 > end) return error.VbaTruncated;
                    const token = bytes.readIntLE(u16, src[pos..][0..2]);
                    pos += 2;
                    const help: CopyTokenHelp = .init(out.items.len - chunk_start);
                    const len = (token & help.lengthMask()) + 3;
                    const offset = (token >> help.offsetShift()) + 1;
                    if (offset > out.items.len - chunk_start) return error.VbaBadCopyToken;
                    for (0..len) |_| try out.append(gpa, out.items[out.items.len - offset]);
                }
            }
        }
    }
    return out.toOwnedSlice(gpa);
}

fn expectRoundTrip(src: []const u8) !void {
    const gpa = zigstd.testAllocator;
    const c = try compress(gpa, src);
    defer gpa.free(c);
    const d = try decompress(gpa, c);
    defer gpa.free(d);
    try zigstd.testExpectEqualSlices(u8, src, d);
}

test "compression round trips" {
    try expectRoundTrip("");
    try expectRoundTrip("a");
    // cspell:disable-next-line
    try expectRoundTrip("#aaabcdefaaaaghijaaaaaklaaamnopqaaaaaaaaaaaarstuvwxyzaaa");
    var big: [10_000]u8 = undefined;
    for (&big, 0..) |*c, i| c.* = "Option Explicit\r\nPublic Sub X()\r\n"[i % 33];
    try expectRoundTrip(&big);
    var noise: [9000]u8 = undefined;
    var prng: std.Random.DefaultPrng = .init(1);
    prng.random().bytes(&noise);
    try expectRoundTrip(&noise);
}
