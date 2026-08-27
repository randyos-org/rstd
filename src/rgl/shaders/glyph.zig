//! Glyph rendering

const std = @import("std");
const spirv = std.spirv;
pub const shader_opts = @import("shader_opts");

const lib = @import("lib/__root__.zig");
const Vec2 = lib.types.Vec2;
const Vec4 = lib.types.Vec4;
const Image = lib.types.Image;
const SampledImage = lib.types.SampledImage;
const param = lib.compile.param;
const input = lib.compile.input;

const param_tex_coord = param("tex_coord", Vec2, 0);
const param_color = param("color", Vec4, 1);

const main = lib.compile.exportShaderMain(&vert, &frag);
comptime { _ = main; }

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
    // Atlas is single-channel coverage in `.r`; tint by the instance
    // color and premultiply.
    const alpha = color[3] * sampled[0];
    lib.compile.frag_color_out.* = lib.coords.rgbPreMulAlpha(lib.coords.rgb(color), alpha);
}
