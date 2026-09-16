//! Glyph rendering

const std = @import("std");
const spirv = std.spirv;
pub const shader_opts = @import("shader_opts");

const shader_mod = @import("rstd").rgl.shaders;
const lib = shader_mod.lib;
const types = shader_mod.types;
const Vec2 = types.Vec2;
const Vec4 = types.Vec4;
const Image = types.Image;
const SampledImage = types.SampledImage;
const param = lib.compile.param;
const input = lib.compile.input;

const param_tex_coord = param("tex_coord", Vec2, 0);
const param_color = param("color", Vec4, 1);

const main = lib.compile.exportShaderMain(&vert, &frag);
comptime {
    _ = main;
}

const atlas = lib.compile.constImage("atlas", 0, 0);

const in_glyph_pos = input("glyph_pos", Vec2, 0);
const in_glyph_size = input("glyph_size", Vec2, 1);
const in_atlas_pos = input("atlas_pos", Vec2, 2);
const in_color = input("instance_color", Vec4, 3);

fn vert() callconv(.spirv_vertex) void {
    const corner = lib.coords.quadCorner(spirv.vertex_index);
    const pixel_pos = in_glyph_pos.* + in_glyph_size.* * corner;
    spirv.position_out.* = lib.coords.fromScreenToNdc(pixel_pos, lib.globals.screen_size);

    // Normalized [0,1]: the sampler uses normalized coordinates, since
    // implicit-Lod sampling (the only variant this Zig version exposes)
    // is illegal with unnormalized ones.
    param_tex_coord.* = (in_atlas_pos.* + in_glyph_size.* * corner) / lib.globals.atlas_size;
    param_color.* = in_color.*;
}

fn frag() callconv(.{ .spirv_fragment = .{} }) void {
    const sampled = spirv.imageSampleImplicitLod(atlas, param_tex_coord.*);
    const color = param_color.*;
    const px_range = lib.globals.glyph_msdf_px_range;

    // 0 (Globals' default) means the bound atlas is plain single-channel
    // coverage - sample .r directly as alpha, same as before. Otherwise
    // it's an msdf.Generator atlas: median(r,g,b) recovers the signed
    // distance, thresholded against px_range for an antialiased edge. This
    // px_range*(distance-0.5) form is only correct at 1:1 texel:screen-pixel
    // scale (true here, since a glyph's on-screen quad is exactly the size
    // buildTextAtlas rasterized it at) - a renderer that scales glyphs up
    // or down would need the fwidth-based screenPxRange() correction
    // msdfgen's own example shader uses instead.
    const coverage = if (px_range == 0) sampled[0] else blk: {
        const distance = lib.math.median3(sampled[0], sampled[1], sampled[2]);
        const screen_px_distance = px_range * (distance - 0.5);
        break :blk lib.math.clamp(screen_px_distance + 0.5, 0.0, 1.0);
    };

    const alpha = color[3] * coverage;
    lib.compile.frag_color_out.* = lib.coords.rgbPreMulAlpha(lib.coords.rgb(color), alpha);
}
