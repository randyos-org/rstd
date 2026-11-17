pub const lib = @import("lib/__root__.zig");
pub const types = @import("types.zig");
pub const layout = @import("layout.zig");

pub const has_shader_opts = @hasDecl(@import("root"), "shader_opts");

pub const rgl_vert = if (has_shader_opts) null else spirv("rgl.shapes.vert.spv");
pub const rgl_frag = if (has_shader_opts) null else spirv("rgl.shapes.frag.spv");

pub const glyph_vert = if (has_shader_opts) null else spirv("rgl.glyph.vert.spv");
pub const glyph_frag = if (has_shader_opts) null else spirv("rgl.glyph.frag.spv");

pub const image_vert = if (has_shader_opts) null else spirv("rgl.image.vert.spv");
pub const image_frag = if (has_shader_opts) null else spirv("rgl.image.frag.spv");

/// SPIR-V `name` on a word boundary: `vkCreateShaderModule` reads words, and
/// `@embedFile` promises bytes only.
fn spirv(comptime name: []const u8) []align(4) const u8 {
    return &struct {
        const words: [@embedFile(name).len]u8 align(4) = @embedFile(name).*;
    }.words;
}
