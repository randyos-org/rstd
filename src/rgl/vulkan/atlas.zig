//! The glyph atlas `drawText` lays text out on: the embedded VGA font as one
//! row of 256 glyphs, one coverage byte per texel. No Vulkan here, so any
//! consumer can build or test it.

const vga = @import("../font/vga/__root__.zig");
const zigstd = @import("../../zigstd.zig");
const Allocator = zigstd.Allocator;

pub const width: u32 = vga.glyph_count * vga.width;
pub const height: u32 = vga.height;

/// The atlas's texels, `width` by `height`, row by row; the caller frees them.
pub fn vgaAtlas(gpa: Allocator) Allocator.Error![]u8 {
    const atlas = try gpa.alloc(u8, width * height);
    for (0..vga.glyph_count) |cp| for (vga.glyph(@intCast(cp)), 0..) |bits, row| {
        for (0..vga.width) |col| {
            const on = (bits >> @intCast(7 - col)) & 1 != 0;
            atlas[row * width + cp * vga.width + col] = if (on) 255 else 0;
        }
    };
    return atlas;
}

test "the VGA atlas puts glyph c at column c * 8, and a space is empty" {
    const gpa = zigstd.testAllocator;
    const atlas = try vgaAtlas(gpa);
    defer gpa.free(atlas);
    try zigstd.testExpectEqual(@as(usize, width * height), atlas.len);
    // 'A' has ink somewhere in its cell, ' ' none
    var ink_a: usize = 0;
    var ink_space: usize = 0;
    for (0..height) |row| for (0..vga.width) |col| {
        if (atlas[row * width + @as(usize, 'A') * vga.width + col] != 0) ink_a += 1;
        if (atlas[row * width + @as(usize, ' ') * vga.width + col] != 0) ink_space += 1;
    };
    try zigstd.testExpectTrue(ink_a > 0);
    try zigstd.testExpectEqual(@as(usize, 0), ink_space);
}
