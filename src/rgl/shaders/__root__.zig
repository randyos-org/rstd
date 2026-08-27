pub const lib = @import("lib/__root__.zig");

pub const rgl_vert = @embedFile("rgl.vert.spv");
pub const rgl_frag = @embedFile("rgl.frag.spv");

pub const glyph_vert = @embedFile("rgl.glyph.vert.spv");
pub const glyph_frag = @embedFile("rgl.glyph.frag.spv");

pub const image_vert = @embedFile("rgl.image.vert.spv");
pub const image_frag = @embedFile("rgl.image.frag.spv");
