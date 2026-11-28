//! A sampled 2D texture in device memory, uploaded through a staging buffer
//! and left in `shader_read_only_optimal`; its sampler and descriptor set.
//! Analyzed only with `-Dvulkan`, as `gpu.zig`.

const vk = @import("vulkan");
const gpu = @import("gpu.zig");

const Texture = @This();

image: vk.Image,
mem: vk.DeviceMemory,
view: vk.ImageView,
width: u32,
height: u32,

/// What an upload needs: the device, its memory types, and a queue with a
/// command pool on its family.
pub const Context = struct {
    dev: vk.DeviceProxy,
    mem_props: vk.PhysicalDeviceMemoryProperties,
    queue: vk.Queue,
    pool: vk.CommandPool,
};

const color_range: vk.ImageSubresourceRange = .{ .aspect_mask = .{ .color = true }, .base_mip_level = 0, .level_count = 1, .base_array_layer = 0, .layer_count = 1 };
const color_layers: vk.ImageSubresourceLayers = .{ .aspect_mask = .{ .color = true }, .mip_level = 0, .base_array_layer = 0, .layer_count = 1 };

/// A `width` by `height` texture of `format` from `bytes` (tightly packed);
/// waits for the upload, so call it outside a render pass.
pub fn create(ctx: Context, format: vk.Format, width: u32, height: u32, bytes: []const u8) !Texture {
    const dev = ctx.dev;
    const image = try dev.createImage(&.{
        .image_type = .@"2d",
        .format = format,
        .extent = .{ .width = width, .height = height, .depth = 1 },
        .mip_levels = 1,
        .array_layers = 1,
        .samples = .{ .@"1" = true },
        .tiling = .optimal,
        .usage = .{ .sampled = true, .transfer_dst = true },
        .sharing_mode = .exclusive,
        .initial_layout = .undefined,
    }, null);
    errdefer dev.destroyImage(image, null);
    const reqs = dev.getImageMemoryRequirements(image);
    const kind = gpu.findMemoryType(ctx.mem_props, reqs, .{ .device_local = true }) orelse return error.NoSuitableMemoryType;
    const mem = try dev.allocateMemory(&.{ .allocation_size = reqs.size, .memory_type_index = kind }, null);
    errdefer dev.freeMemory(mem, null);
    try dev.bindImageMemory(image, mem, 0);
    const view = try dev.createImageView(&.{
        .image = image,
        .view_type = .@"2d",
        .format = format,
        .components = .{ .r = .identity, .g = .identity, .b = .identity, .a = .identity },
        .subresource_range = color_range,
    }, null);
    errdefer dev.destroyImageView(view, null);
    const t: Texture = .{ .image = image, .mem = mem, .view = view, .width = width, .height = height };
    try t.upload(ctx, bytes);
    return t;
}

pub fn deinit(t: Texture, dev: vk.DeviceProxy) void {
    dev.destroyImageView(t.view, null);
    dev.destroyImage(t.image, null);
    dev.freeMemory(t.mem, null);
}

/// Copies `bytes` in through a staging buffer and one command buffer.
fn upload(t: Texture, ctx: Context, bytes: []const u8) !void {
    const dev = ctx.dev;
    var staging: gpu.InstanceBuffer = .{};
    defer staging.deinit(dev);
    try staging.uploadUsage(dev, ctx.mem_props, bytes, .{ .transfer_src = true });
    var cmd: [1]vk.CommandBuffer = undefined;
    try dev.allocateCommandBuffers(&.{ .command_pool = ctx.pool, .level = .primary, .command_buffer_count = 1 }, &cmd);
    defer dev.freeCommandBuffers(ctx.pool, &cmd);
    try dev.beginCommandBuffer(cmd[0], &.{ .flags = .{ .one_time_submit = true } });
    barrier(dev, cmd[0], t.image, .undefined, .transfer_dst_optimal);
    const region = [_]vk.BufferImageCopy{.{
        .buffer_offset = 0,
        .buffer_row_length = 0,
        .buffer_image_height = 0,
        .image_subresource = color_layers,
        .image_offset = .{ .x = 0, .y = 0, .z = 0 },
        .image_extent = .{ .width = t.width, .height = t.height, .depth = 1 },
    }};
    dev.cmdCopyBufferToImage(cmd[0], staging.buf, t.image, .transfer_dst_optimal, &region);
    barrier(dev, cmd[0], t.image, .transfer_dst_optimal, .shader_read_only_optimal);
    try dev.endCommandBuffer(cmd[0]);
    const info = [_]vk.SubmitInfo{.{ .command_buffer_count = 1, .p_command_buffers = &cmd }};
    try dev.queueSubmit(ctx.queue, &info, .null);
    try dev.queueWaitIdle(ctx.queue);
}

/// A layout transition of `img`, ordered after every earlier command.
pub fn barrier(dev: vk.DeviceProxy, cmd: vk.CommandBuffer, img: vk.Image, from: vk.ImageLayout, to: vk.ImageLayout) void {
    const all: vk.AccessFlags = .{ .memory_read = true, .memory_write = true };
    const b = [_]vk.ImageMemoryBarrier{.{
        .src_access_mask = all,
        .dst_access_mask = all,
        .old_layout = from,
        .new_layout = to,
        .src_queue_family_index = vk.QUEUE_FAMILY_IGNORED,
        .dst_queue_family_index = vk.QUEUE_FAMILY_IGNORED,
        .image = img,
        .subresource_range = color_range,
    }};
    dev.cmdPipelineBarrier(cmd, .{ .all_commands = true }, .{ .all_commands = true }, .{}, null, null, &b);
}

/// The nearest-texel, clamped sampler the glyph and image pipelines read with.
pub fn nearestSampler(dev: vk.DeviceProxy) !vk.Sampler {
    return dev.createSampler(&.{
        .mag_filter = .nearest,
        .min_filter = .nearest,
        .mipmap_mode = .nearest,
        .address_mode_u = .clamp_to_edge,
        .address_mode_v = .clamp_to_edge,
        .address_mode_w = .clamp_to_edge,
        .mip_lod_bias = 0,
        .anisotropy_enable = .false,
        .max_anisotropy = 1,
        .compare_enable = .false,
        .compare_op = .always,
        .min_lod = 0,
        .max_lod = 0,
        .border_color = .float_transparent_black,
        .unnormalized_coordinates = .false,
    }, null);
}

/// A descriptor set from `pool` for `pipeline`'s layout, binding this
/// texture with `sampler`.
pub fn descriptorSet(t: Texture, dev: vk.DeviceProxy, pool: vk.DescriptorPool, pipeline: gpu.Pipeline, sampler: vk.Sampler) !vk.DescriptorSet {
    var set: [1]vk.DescriptorSet = undefined;
    try dev.allocateDescriptorSets(&.{ .descriptor_pool = pool, .descriptor_set_count = 1, .p_set_layouts = &.{pipeline.descriptor_set_layout} }, &set);
    const info = [_]vk.DescriptorImageInfo{.{ .sampler = sampler, .image_view = t.view, .image_layout = .shader_read_only_optimal }};
    const write = [_]vk.WriteDescriptorSet{.{
        .dst_set = set[0],
        .dst_binding = 0,
        .dst_array_element = 0,
        .descriptor_count = 1,
        .descriptor_type = .combined_image_sampler,
        .p_image_info = &info,
        .p_buffer_info = &.{},
        .p_texel_buffer_view = &.{},
    }};
    dev.updateDescriptorSets(&write, null);
    return set[0];
}
