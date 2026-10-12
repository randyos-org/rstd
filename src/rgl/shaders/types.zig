const std = @import("std");
const spirv = std.spirv;

const rsimd = @import("../../simd/__root__.zig");
const zigvec = rsimd.zigvec;

pub const V2f32 = zigvec.V2f32;
pub const V3f32 = zigvec.V3f32;
pub const V4f32 = zigvec.V4f32;
pub const splatVec = zigvec.splatVecF32;

pub const ShaderStage = @import("../../buildutils/shadertypes.zig").Stage;

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
