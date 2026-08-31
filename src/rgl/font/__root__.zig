// pub const TrueType = @import("TrueType.zig");
// const tt = try TrueType.load(@embedFile("GoNotoCurrent-Regular.ttf"));

// const glyph = tt.codepointGlyphIndex('A');
// const scale = tt.scaleForPixelHeight(20.0);

// var pixels: std.ArrayListUnmanaged(u8) = .empty;
// defer pixels.deinit(alloc);
// const bmp = try tt.glyphBitmap(alloc, &pixels, glyph, scale, scale);
// // pixels.items[0 .. bmp.width * bmp.height] is the 8-bit coverage bitmap;
// // bmp.off_x/off_y is where it sits relative to the glyph origin.

// var font = try Font.init(@embedFile("GoNotoCurrent-Regular.ttf"));
// defer font.deinit(alloc);

// const atlas = try font.generateAtlas(alloc, "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz", .{
//     .type = .msdf,
//     .pixel_size = 32,
// });
// defer atlas.deinit(alloc);
// atlas.pixels is the packed RGB texture; atlas.glyphs[i].tex_bounds gives the
// UV rectangle within it for each requested codepoint, alongside its metrics.
