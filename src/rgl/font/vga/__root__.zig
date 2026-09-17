//! Embedded 8x16 VGA bitmap font (CP437), for cases that want a fixed
//! console-style font.
//!
//! IBM_VGA_8x16.bin is derived from IBM MCGA/VGA (PS/2) video BIOS,
//! 1 byte/row MSB-first (bit 7 = leftmost pixel).

pub const width: u8 = 8;
pub const height: u8 = 16;
pub const glyph_count: usize = 256;

const data = @embedFile("IBM_VGA_8x16.bin");

const zigstd = @import("../../../zigstd.zig");

comptime {
    zigstd.assert(data.len == glyph_count * height);
}

/// The 16 row-bytes for CP437 code point `cp`, MSB-first (bit 7 leftmost).
pub fn glyph(cp: u8) *const [height]u8 {
    return data[@as(usize, cp) * height ..][0..height];
}

test "glyph(0x41) is the letter A" {
    try zigstd.testExpectEqualSlices(u8, &[height]u8{
        0,    0,    0x10, 0x38, 0x6c, 0xc6, 0xc6, 0xfe,
        0xc6, 0xc6, 0xc6, 0xc6, 0,    0,    0,    0,
    }, glyph('A'));
}

test "glyph(0x20) (space) is blank" {
    const blank: [height]u8 = @splat(0);
    try zigstd.testExpectEqualSlices(u8, &blank, glyph(' '));
}
