// const sfnt = try Sfnt.init(alloc, @embedFile("GoNotoCurrent-Regular.ttf"), .{});
// const glyph = sfnt.getGlyphIndex('A') orelse ...;
// const scale = sfnt.scaleForPixelHeight(20.0);
//
// var pixels: std.ArrayListUnmanaged(u8) = .empty;
// defer pixels.deinit(alloc);
// const bmp = try sfnt.glyphBitmap(alloc, &pixels, glyph, scale, scale);
// // pixels.items[0 .. bmp.width * bmp.height] is the 8-bit coverage bitmap;
// // bmp.off_x/off_y is where it sits relative to the glyph origin.

// var atlas = try sfnt.buildTextAtlas(alloc, "Hello, world!", 32, 0, 0, .{ .type = .msdf });
// defer atlas.deinit(alloc);
// atlas.pixels is one shared texture (kerning-aware, left to right);
// atlas.glyphs[i].{pos,size,atlas_pos} place each glyph both on screen
// and within the texture.

const std = @import("std");

pub const Sfnt = @import("Sfnt.zig");
pub const render = @import("render/__root__.zig");
pub const types = @import("stream_types.zig");
pub const aux_types = @import("aux_types.zig");
pub const tables = @import("tables/__root__.zig");
pub const msdf = @import("msdf/__root__.zig");

pub const Vertex = aux_types.Vertex;
pub const BitmapBox = render.bounds.BitmapBox;
pub const CharstringCtx = tables.CharstringCtx;

test {
    std.testing.refAllDecls(@This());
}
