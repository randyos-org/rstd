//! A test `gvar`: one axis, one glyph, two tuples peaking at 1.0.

const bytes = @import("../../../../../bytes.zig");

/// A gvar for one axis and one glyph with two tuples, both peaking at 1:
/// points 0 and 2 moved by (10, 4) and (20, 8), then every point by x 1.
pub fn twoTuples(buf: []u8) []u8 {
    const tuple1 = [_]u8{
        2, 0x01, 0, 2, // private points 0 and 2
        0x01, 10, 20, // x deltas
        0x01, 4, 8, // y deltas
    };
    const tuple2 = [_]u8{
        0, // every point
        0x07, 1, 1, 1, 1, 1, 1, 1, 1, // x deltas
        0x87, // y deltas, zero
    };
    const array_off = 20 + 8;
    const gvd_header = 4 + 2 * 6;
    @memset(buf, 0);
    bytes.setU16BE(buf[0..], 1); // majorVersion
    bytes.setU16BE(buf[4..], 1); // axisCount
    bytes.setU32BE(buf[8..], array_off); // sharedTuplesOffset (none)
    bytes.setU16BE(buf[12..], 1); // glyphCount
    bytes.setU16BE(buf[14..], 1); // flags: long offsets
    bytes.setU32BE(buf[16..], array_off);
    const total = gvd_header + tuple1.len + tuple2.len;
    bytes.setU32BE(buf[24..], total);
    const g = buf[array_off..];
    bytes.setU16BE(g[0..], 2); // tupleVariationCount
    bytes.setU16BE(g[2..], gvd_header); // dataOffset
    for ([_]usize{ tuple1.len, tuple2.len }, 0..) |size, k| {
        bytes.setU16BE(g[4 + k * 6 ..], @intCast(size));
        bytes.setU16BE(g[6 + k * 6 ..], 0xa000); // embedded peak, private points
        bytes.setU16BE(g[8 + k * 6 ..], 0x4000); // peak 1.0
    }
    @memcpy(g[gvd_header..][0..tuple1.len], &tuple1);
    @memcpy(g[gvd_header + tuple1.len ..][0..tuple2.len], &tuple2);
    return buf[0 .. array_off + total];
}
