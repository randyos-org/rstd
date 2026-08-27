//! General-purpose 2D SDF-based primitives

const std = @import("std");
const spirv = std.spirv;
pub const shader_opts = @import("shader_opts");
const lib = @import("lib/__root__.zig");

const Vec2 = lib.types.Vec2;
const Vec4 = lib.types.Vec4;
const splatVec = lib.types.splatVec;
const param = lib.compile.param;
const paramFlat = lib.compile.paramFlat;
const input = lib.compile.input;

pub const ShapeKind = enum(u32) {
    /// `corner_radius` 0 is `Rectangle`, nonzero is `RoundRect`.
    rect = 0,
    /// `Ellipse`, restricted to circles; see common.zig.
    circle = 1,
};

/// Half-width, in pixels, of the antialiasing band around an edge.
const aa_width: f32 = 0.75;

// local_pos varies per pixel; the rest is constant per instance.
// shape_kind is an integer and must be flat.
const param_local_pos = param("local_pos", Vec2, 0);
const param_half_size = paramFlat("half_size", Vec2, 1);
const param_corner_radius = paramFlat("corner_radius", f32, 2);
const param_border_width = paramFlat("border_width", f32, 3);
const param_fill_color = paramFlat("fill_color", Vec4, 4);
const param_border_color = paramFlat("border_color", Vec4, 5);
const param_shape_kind = paramFlat("shape_kind", u32, 6);

const in_rect_pos = input("rect_pos", Vec2, 0);
const in_rect_size = input("rect_size", Vec2, 1);
const in_corner_radius = input("instance_corner_radius", f32, 2);
const in_border_width = input("instance_border_width", f32, 3);
const in_fill_color = input("instance_fill_color", Vec4, 4);
const in_border_color = input("instance_border_color", Vec4, 5);
const in_shape_kind = input("instance_shape_kind", u32, 6);

fn vert() callconv(.spirv_vertex) void {
    const corner = lib.coords.quadCorner(spirv.vertex_index);
    const rect_size = in_rect_size.*;
    const pixel_pos = in_rect_pos.* + rect_size * corner;
    spirv.position_out.* = lib.coords.fromScreenToNdc(pixel_pos, lib.globals.screen_size);

    param_local_pos.* = (corner - splatVec(2, 0.5)) * rect_size;
    param_half_size.* = rect_size * splatVec(2, 0.5);
    param_corner_radius.* = in_corner_radius.*;
    param_border_width.* = in_border_width.*;
    param_fill_color.* = in_fill_color.*;
    param_border_color.* = in_border_color.*;
    param_shape_kind.* = in_shape_kind.*;
}

fn frag() callconv(.{ .spirv_fragment = .{} }) void {
    const p = param_local_pos.*;
    const half_size = param_half_size.*;
    const border_width = param_border_width.*;
    const kind: ShapeKind = @fromBackingInt(@intCast(param_shape_kind.*));

    const dist = switch (kind) {
        .rect => lib.sdf.sdRoundedBox(p, half_size, param_corner_radius.*),
        .circle => lib.sdf.sdCircle(p, half_size[0]),
    };

    // 1 inside, 0 outside, ramped across dist == 0.
    const outer_coverage = 1.0 - lib.math.smoothstep_scalar(-aa_width, aa_width, dist);

    // Border/fill split is a color mix, not an alpha change.
    const base_color = if (border_width > 0.0) blk: {
        const inner_dist = dist + border_width; // 0 at the seam
        const fill_mix = lib.math.smoothstep_scalar(-aa_width, aa_width, inner_dist);
        break :blk lib.math.mix(param_border_color.*, param_fill_color.*, fill_mix);
    } else param_fill_color.*;

    const alpha = base_color[3] * outer_coverage;
    lib.compile.frag_color_out.* = lib.coords.rgbPreMulAlpha(lib.coords.rgb(base_color), alpha);
}

const main = lib.compile.exportShaderMain(&vert, &frag);
comptime { _ = main; }
