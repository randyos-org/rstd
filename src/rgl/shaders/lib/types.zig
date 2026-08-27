const std = @import("std");
const spirv = std.spirv;

pub const ShaderStage = enum { frag, vert };

pub const V2f32 = @Vector(2, f32);
pub const V3f32 = @Vector(3, f32);
pub const V4f32 = @Vector(4, f32);

// aliases to resemble glsl
pub const Float = f32;
pub const Vec2 = V2f32;
pub const Vec3 = V3f32;
pub const Vec4 = V4f32;
// pub const Mat2 = @SpirvType(.{ .runtime_array =  })
pub const Mat2 = [2]Vec2; // @Matrix(2, 2, f32);
pub const Mat3 = [3]Vec3; // @Matrix(3, 3, f32);
pub const Mat4 = [4]Vec4; // @Matrix(4, 4, f32);
pub const Sampler2D = Sampler;
pub const SamplerCube = Sampler;

// improve readability
pub const Rgba = Vec4;

pub inline fn splatVec(N: usize, v: f32) @Vector(N, f32) {
    return @as(@Vector(N, f32), @splat(v));
}

pub const Sampler = @SpirvType(.sampler);
pub const AnyImage = @SpirvType(.{ .image = .{
    .usage = .{ .sampled = u32 },
    .format = .unknown,
    .dim = .@"2d",
    .depth = .unknown,
    .arrayed = false,
    .multisampled = false,
    .access = .unknown,
} });
pub const DepthImage = @SpirvType(.{ .image = .{
    .usage = .{ .sampled = u32 },
    .format = .unknown,
    .dim = .@"2d",
    .depth = .depth,
    .arrayed = false,
    .multisampled = false,
    .access = .unknown,
} });
pub const Image = @SpirvType(.{ .image = .{
    .usage = .{ .sampled = f32 },
    .format = .unknown,
    .dim = .@"2d",
    .depth = .not_depth,
    .arrayed = false,
    .multisampled = false,
    .access = .unknown,
} });
pub const SampledImage = @SpirvType(.{ .sampled_image = Image });

// const RuntimeArray = @SpirvType(.{ .runtime_array = u32 });
// const sampled_image = @extern(*addrspace(.constant) const SampledImage, .{
//     .name = "sampled_image",
//     .decoration = .{ .descriptor = .{ .set = 0, .binding = 1 } },
// });
