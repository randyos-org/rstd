//! A `VkSurface`'s built-in draws on the GPU (#409): rects, glyphs on the
//! VGA atlas (`atlas.zig`), and images, each with its own texture. One per
//! frame in flight: `prepare` before the render pass (it uploads), `record`
//! inside it, `release` once the frame's fence has signaled. Custom
//! pipelines record their own batches (`gpu.recordBatch`). Analyzed only
//! with `-Dvulkan`, as `gpu.zig`.

const std = @import("std");
const vk = @import("vulkan");
const gpu = @import("gpu.zig");
const atlas = @import("atlas.zig");
const rglvk = @import("__root__.zig");
const Texture = @import("Texture.zig");
const RgbaImage = @import("../image.zig").RgbaImage;
const zigstd = @import("../../zigstd.zig");

const Allocator = zigstd.Allocator;
const SurfaceDraws = @This();

sampler: vk.Sampler = .null,
/// The VGA atlas, uploaded once.
glyphs: ?Texture = null,
atlas_pool: vk.DescriptorPool = .null,
atlas_set: vk.DescriptorSet = .null,
/// This frame's textures, one per distinct image, and their sets' pool.
images: std.ArrayList(Image) = .empty,
image_pool: vk.DescriptorPool = .null,
/// Runs of consecutive draws of one image, in draw order.
runs: std.ArrayList(Run) = .empty,
/// Every image draw's instance, in draw order.
instances: std.ArrayList(u8) = .empty,

const Image = struct { source: *const RgbaImage, tex: Texture, set: vk.DescriptorSet };
const Run = struct { image: u32, first: u32, count: u32 };

/// Uploads the atlas (the first time) and `surface`'s images, outside a
/// render pass. `release` the last frame's first.
pub fn prepare(d: *SurfaceDraws, gpa: Allocator, ctx: Texture.Context, pipelines: gpu.Pipelines, surface: *const rglvk.VkSurface) !void {
    const dev = ctx.dev;
    zigstd.assert(d.images.items.len == 0 and d.runs.items.len == 0);
    if (d.sampler == .null) d.sampler = try Texture.nearestSampler(dev);
    if (d.glyphs == null) {
        const bytes = try atlas.vgaAtlas(gpa);
        defer gpa.free(bytes);
        const tex = try Texture.create(ctx, .r8_unorm, atlas.width, atlas.height, bytes);
        errdefer tex.deinit(dev);
        d.atlas_pool = try onePool(dev, 1);
        d.atlas_set = try tex.descriptorSet(dev, d.atlas_pool, pipelines.glyph, d.sampler);
        d.glyphs = tex;
    }
    const draws = surface.images.items;
    if (draws.len == 0) return;
    errdefer d.release(dev);
    for (draws, 0..) |draw, i| {
        const at = for (d.images.items, 0..) |im, k| {
            if (im.source == draw.image) break k;
        } else blk: {
            const img = draw.image;
            const tex = try Texture.create(ctx, .r8g8b8a8_unorm, img.width, img.height, std.mem.sliceAsBytes(img.data));
            errdefer tex.deinit(dev);
            try d.images.append(gpa, .{ .source = img, .tex = tex, .set = .null });
            break :blk d.images.items.len - 1;
        };
        const runs = d.runs.items;
        if (runs.len > 0 and runs[runs.len - 1].image == at) {
            runs[runs.len - 1].count += 1;
        } else try d.runs.append(gpa, .{ .image = @intCast(at), .first = @intCast(i), .count = 1 });
        try d.instances.appendSlice(gpa, std.mem.asBytes(&draw.instance));
    }
    d.image_pool = try onePool(dev, @intCast(d.images.items.len));
    for (d.images.items) |*im| im.set = try im.tex.descriptorSet(dev, d.image_pool, pipelines.image, d.sampler);
}

/// Records `surface`'s rects, glyphs and images, in that order, inside an
/// active render pass whose viewport and scissor the caller set.
pub fn record(d: *const SurfaceDraws, dev: vk.DeviceProxy, mem_props: vk.PhysicalDeviceMemoryProperties, cmd: vk.CommandBuffer, pipelines: gpu.Pipelines, res: *gpu.FrameResources, surface: *const rglvk.VkSurface, screen_size: [2]f32) !void {
    const globals: gpu.Globals = .{ .screen_size = screen_size, .atlas_size = .{ atlas.width, atlas.height } };
    try gpu.recordBatch(dev, mem_props, cmd, pipelines.rect, &res.rect, surface.batch(rglvk.rect_pipeline_id), globals, null);
    try gpu.recordBatch(dev, mem_props, cmd, pipelines.glyph, &res.glyph, surface.batch(rglvk.glyph_pipeline_id), globals, d.atlas_set);
    if (d.runs.items.len == 0) return;
    try res.image.upload(dev, mem_props, d.instances.items);
    const image = pipelines.image;
    dev.cmdBindPipeline(cmd, .graphics, image.handle);
    dev.cmdPushConstants(cmd, image.layout, .{ .vertex = true, .fragment = true }, 0, @sizeOf(gpu.Globals), &globals);
    const offsets = [_]vk.DeviceSize{0};
    dev.cmdBindVertexBuffers(cmd, 0, (&res.image.buf)[0..1], &offsets);
    for (d.runs.items) |run| {
        dev.cmdBindDescriptorSets(cmd, .graphics, image.layout, 0, (&d.images.items[run.image].set)[0..1], null);
        dev.cmdDraw(cmd, 4, run.count, 0, run.first);
    }
}

/// Frees this frame's image textures, once the GPU is done with them.
pub fn release(d: *SurfaceDraws, dev: vk.DeviceProxy) void {
    for (d.images.items) |im| im.tex.deinit(dev);
    d.images.clearRetainingCapacity();
    d.runs.clearRetainingCapacity();
    d.instances.clearRetainingCapacity();
    if (d.image_pool != .null) dev.destroyDescriptorPool(d.image_pool, null);
    d.image_pool = .null;
}

pub fn deinit(d: *SurfaceDraws, gpa: Allocator, dev: vk.DeviceProxy) void {
    d.release(dev);
    d.images.deinit(gpa);
    d.runs.deinit(gpa);
    d.instances.deinit(gpa);
    if (d.glyphs) |t| t.deinit(dev);
    if (d.atlas_pool != .null) dev.destroyDescriptorPool(d.atlas_pool, null);
    if (d.sampler != .null) dev.destroySampler(d.sampler, null);
    d.* = .{};
}

/// A pool of `sets` combined-image-sampler sets.
fn onePool(dev: vk.DeviceProxy, sets: u32) !vk.DescriptorPool {
    const sizes = [_]vk.DescriptorPoolSize{.{ .type = .combined_image_sampler, .descriptor_count = sets }};
    return dev.createDescriptorPool(&.{ .max_sets = sets, .pool_size_count = sizes.len, .p_pool_sizes = &sizes }, null);
}
