pub const lib = @import("lib/__root__.zig");
pub const types = @import("types.zig");
pub const layout = @import("layout.zig");

pub const has_shader_opts = @hasDecl(@import("root"), "shader_opts");

pub const rgl_vert = if (has_shader_opts) null else @embedFile("rgl.shapes.vert.spv");
pub const rgl_frag = if (has_shader_opts) null else @embedFile("rgl.shapes.frag.spv");

pub const glyph_vert = if (has_shader_opts) null else @embedFile("rgl.glyph.vert.spv");
pub const glyph_frag = if (has_shader_opts) null else @embedFile("rgl.glyph.frag.spv");

pub const image_vert = if (has_shader_opts) null else @embedFile("rgl.image.vert.spv");
pub const image_frag = if (has_shader_opts) null else @embedFile("rgl.image.frag.spv");
