//! One textured quad per instance, normalized [0,1] UVs.
//! Each draw binds its own texture.

const std = @import("std");
const spirv = std.spirv;
pub const shader_opts = @import("shader_opts");

const shader_mod = @import("rstd").rgl.shaders;
const lib = shader_mod.lib;
const types = shader_mod.types;
const Vec2 = types.Vec2;
const Vec4 = types.Vec4;
const param = lib.compile.param;
const input = lib.compile.input;

const param_uv = param("uv", Vec2, 0);
const param_color = param("color", Vec4, 1);

const main = lib.compile.exportShaderMain(struct {
    pub const vert = vertMain;
    pub const frag = fragMain;
});
comptime {
    _ = main;
}

const image = lib.compile.constImage("image", 0, 0);

const in_rect_pos = input("rect_pos", Vec2, 0);
const in_rect_size = input("rect_size", Vec2, 1);
const in_color = input("instance_color", Vec4, 2);

fn vertMain() callconv(.spirv_vertex) void {
    const corner = lib.coords.quadCorner(spirv.vertex_index);
    const pixel_pos = in_rect_pos.* + in_rect_size.* * corner;
    spirv.position_out.* = lib.coords.fromScreenToNdc(pixel_pos, lib.globals.screen_size);

    param_uv.* = corner;
    param_color.* = in_color.*;
}

fn fragMain() callconv(.{ .spirv_fragment = .{} }) void {
    // Sampled image holds straight alpha; premultiply to match rgl.zig.
    const sampled = spirv.imageSampleImplicitLod(image, param_uv.*);
    const color = param_color.*;
    const straight = sampled * color;
    const alpha = straight[3];
    lib.compile.frag_color_out.* = lib.coords.rgbPreMulAlpha(lib.coords.rgb(straight), alpha);
}
