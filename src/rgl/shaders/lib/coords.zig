const std = @import("std");
const types = @import("../types.zig");
const Vec2 = types.Vec2;
const Vec3 = types.Vec3;
const Vec4 = types.Vec4;
const splatVec = types.splatVec;

/// Corner of an instanced quad from `vertex_index` alone; no vertex
/// buffer needed. 0..3 form a triangle strip: (0,0) (1,0) (0,1) (1,1).
pub fn quadCorner(vertex_index: u32) Vec2 {
    const x: f32 = @floatFromInt(vertex_index & 1);
    const y: f32 = @floatFromInt((vertex_index >> 1) & 1);
    return .{ x, y };
}

/// Pixel-space position to Vulkan clip space.
/// Vulkan's Normalized Device Coordinates (NDC) convention is already
/// x-right/y-down, matching screen coordinates; no axis flip needed.
pub fn fromScreenToNdc(pixel_pos: Vec2, screen_size: Vec2) Vec4 {
    const ndc = pixel_pos / screen_size * splatVec(2, 2) - splatVec(2, 1);
    return .{ ndc[0], ndc[1], 0, 1 };
}

pub inline fn xy(v: anytype) Vec2 {
    return swizzle2(v, 0, 1);
}

pub inline fn zw(v: Vec4) Vec2 {
    return swizzle2(v, 2, 3);
}

pub inline fn rgb(v: anytype) Vec3 {
    comptime std.debug.assert(@sizeOf(@TypeOf(v)) >= 3 * @sizeOf(@TypeOf(v[0])));
    return swizzle3(v, 0, 1, 2);
}

pub inline fn swizzle2(v: anytype, a: usize, b: usize) Vec2 {
    return .{ v[a], v[b] };
}

pub inline fn swizzle3(v: anytype, a: usize, b: usize, c: usize) Vec3 {
    comptime std.debug.assert(@sizeOf(@TypeOf(v)) >= 3 * @sizeOf(@TypeOf(v[0])));
    return .{ v[a], v[b], v[c] };
}

pub inline fn rgbPreMulAlpha(rgb_vec: Vec3, alpha: f32) Vec4 {
    const v: Vec4 = .{ rgb_vec[0], rgb_vec[1], rgb_vec[2], 1 };
    return v * splatVec(4, alpha);
}
