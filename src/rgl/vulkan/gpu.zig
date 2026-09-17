//! Real Vulkan-API glue for the batches `VkSurface` (__root__.zig)
//! accumulates: pipeline objects and command-buffer recording. Compiled
//! only when something actually imports it -- see this module's `@import`
//! of the vulkan-zig bindings below, only available when rstd/rgui is
//! built with `-Dvulkan=true` (src/buildutils/vulkan.zig).
//!
//! Deliberately does NOT create or own a VkInstance, physical/logical
//! device, VkSurfaceKHR, swapchain, render pass, framebuffers, command
//! pool, or any sync objects (semaphores/fences) -- all of that stays the
//! caller's responsibility. Every function here takes an already-created
//! `vk.DeviceProxy` and, where relevant, an already-active render pass
//! instance (a `vk.CommandBuffer` mid-`vkCmdBeginRenderPass`) to record
//! into.
//!
//! CAVEAT: this file's `@import("vulkan")` resolves to whichever vk-zig
//! module rstd's own build.zig generated (see addVulkan there). If a
//! consumer (rgui) pulls rstd in as a `b.dependency(...)` -- as it does
//! today -- that's a *different* generated module than the one rgui's own
//! top-level code imports, even though both come from the same generator
//! commit: vk-zig's types are nominal, so e.g. `vk.Device` from here and
//! `vk.Device` from rgui's src/render/vulkan/vulkan.zig are NOT the same
//! type, and handles can't cross that boundary as-is. Fixing that means
//! either folding rstd into rgui's build in-graph (like wireRglShaders
//! already does for the SPIR-V sub-target, reusing one `*Module`) instead
//! of fetching it as an isolated dependency, or moving this file's
//! contents into rgui's own tree instead of rstd's. Left for whoever
//! reconciles the real swapchain/device code to decide.

const std = @import("std");
const vk = @import("vulkan");

const batch = @import("__root__.zig");
const RectInstance = batch.RectInstance;
const GlyphInstance = batch.GlyphInstance;
const ImageInstance = batch.ImageInstance;
const VkSurface = batch.VkSurface;

const shaders = @import("../shaders/__root__.zig");

/// Push constant shared by all three pipelines, matching every shader's
/// `lib.globals.*` reads (shaders/shaders/lib/__root__.zig's `Globals`).
/// Caller fills `screen_size` from the current swapchain extent each frame.
pub const Globals = extern struct {
    screen_size: [2]f32,
    atlas_size: [2]f32 = .{ 0, 0 },
    glyph_msdf_px_range: f32 = 0,
};

pub const Pipeline = struct {
    handle: vk.Pipeline,
    layout: vk.PipelineLayout,
    descriptor_set_layout: vk.DescriptorSetLayout = .null_handle,

    pub fn destroy(self: Pipeline, dev: vk.DeviceProxy) void {
        dev.destroyPipeline(self.handle, null);
        dev.destroyPipelineLayout(self.layout, null);
        if (self.descriptor_set_layout != .null_handle) dev.destroyDescriptorSetLayout(self.descriptor_set_layout, null);
    }
};

fn vertexAttr(location: u32, format: vk.Format, offset: u32) vk.VertexInputAttributeDescription {
    return .{ .binding = 0, .location = location, .format = format, .offset = offset };
}

pub const rect_attrs = [_]vk.VertexInputAttributeDescription{
    vertexAttr(0, .r32g32_sfloat, @offsetOf(RectInstance, "rect_pos")),
    vertexAttr(1, .r32g32_sfloat, @offsetOf(RectInstance, "rect_size")),
    vertexAttr(2, .r32_sfloat, @offsetOf(RectInstance, "rotation")),
    vertexAttr(3, .r32_sfloat, @offsetOf(RectInstance, "corner_radius")),
    vertexAttr(4, .r32_sfloat, @offsetOf(RectInstance, "border_width")),
    vertexAttr(5, .r32g32b32a32_sfloat, @offsetOf(RectInstance, "fill_color")),
    vertexAttr(6, .r32g32b32a32_sfloat, @offsetOf(RectInstance, "border_color")),
    vertexAttr(7, .r32_uint, @offsetOf(RectInstance, "shape_kind")),
};

pub const glyph_attrs = [_]vk.VertexInputAttributeDescription{
    vertexAttr(0, .r32g32_sfloat, @offsetOf(GlyphInstance, "glyph_pos")),
    vertexAttr(1, .r32g32_sfloat, @offsetOf(GlyphInstance, "glyph_size")),
    vertexAttr(2, .r32g32_sfloat, @offsetOf(GlyphInstance, "atlas_pos")),
    vertexAttr(3, .r32g32b32a32_sfloat, @offsetOf(GlyphInstance, "color")),
};

pub const image_attrs = [_]vk.VertexInputAttributeDescription{
    vertexAttr(0, .r32g32_sfloat, @offsetOf(ImageInstance, "rect_pos")),
    vertexAttr(1, .r32g32_sfloat, @offsetOf(ImageInstance, "rect_size")),
    vertexAttr(2, .r32g32b32a32_sfloat, @offsetOf(ImageInstance, "color")),
};

/// Premultiplied-alpha "over" blend, shared by all three pipelines -- see
/// shaders/rgl.zig's header for why they output premultiplied alpha.
const premultiplied_blend: vk.PipelineColorBlendAttachmentState = .{
    .blend_enable = .true,
    .src_color_blend_factor = .one,
    .dst_color_blend_factor = .one_minus_src_alpha,
    .color_blend_op = .add,
    .src_alpha_blend_factor = .one,
    .dst_alpha_blend_factor = .one_minus_src_alpha,
    .alpha_blend_op = .add,
    .color_write_mask = .{ .r = true, .g = true, .b = true, .a = true },
};

/// Builds one instanced-quad pipeline (triangle-strip, 4 vertices, no
/// vertex buffer of its own -- `quadCorner` derives the corner from
/// `spirv.vertex_index`). `stride`/`attrs` describe the per-instance data
/// bound at binding 0 (`input_rate = .instance`); `needs_sampler` adds a
/// single combined-image-sampler descriptor set (glyph/image; rect needs
/// none). Caller owns `render_pass` and destroys the returned `Pipeline`
/// via `.destroy()` before destroying the device.
pub fn createPipeline(
    dev: vk.DeviceProxy,
    render_pass: vk.RenderPass,
    vert_spv: []const u8,
    frag_spv: []const u8,
    attrs: []const vk.VertexInputAttributeDescription,
    stride: u32,
    needs_sampler: bool,
) !Pipeline {
    var descriptor_set_layout: vk.DescriptorSetLayout = .null_handle;
    if (needs_sampler) {
        const binding: [1]vk.DescriptorSetLayoutBinding = .{.{
            .binding = 0,
            .descriptor_type = .combined_image_sampler,
            .descriptor_count = 1,
            .stage_flags = .{ .fragment = true },
        }};
        descriptor_set_layout = try dev.createDescriptorSetLayout(&.{
            .binding_count = binding.len,
            .p_bindings = &binding,
        }, null);
    }
    errdefer if (descriptor_set_layout != .null_handle) dev.destroyDescriptorSetLayout(descriptor_set_layout, null);

    const push_constant_range: [1]vk.PushConstantRange = .{.{
        .stage_flags = .{ .vertex = true, .fragment = true },
        .offset = 0,
        .size = @sizeOf(Globals),
    }};

    const layout = try dev.createPipelineLayout(&.{
        .set_layout_count = if (needs_sampler) 1 else 0,
        .p_set_layouts = if (needs_sampler) (&descriptor_set_layout)[0..1] else null,
        .push_constant_range_count = push_constant_range.len,
        .p_push_constant_ranges = &push_constant_range,
    }, null);
    errdefer dev.destroyPipelineLayout(layout, null);

    const vert = try dev.createShaderModule(&.{ .code_size = vert_spv.len, .p_code = @ptrCast(@alignCast(vert_spv.ptr)) }, null);
    defer dev.destroyShaderModule(vert, null);
    const frag = try dev.createShaderModule(&.{ .code_size = frag_spv.len, .p_code = @ptrCast(@alignCast(frag_spv.ptr)) }, null);
    defer dev.destroyShaderModule(frag, null);

    const stages: [2]vk.PipelineShaderStageCreateInfo = .{
        .{ .stage = .{ .vertex = true }, .module = vert, .p_name = "main" },
        .{ .stage = .{ .fragment = true }, .module = frag, .p_name = "main" },
    };

    const binding_desc: [1]vk.VertexInputBindingDescription = .{.{
        .binding = 0,
        .stride = stride,
        .input_rate = .instance,
    }};

    const vertex_input: vk.PipelineVertexInputStateCreateInfo = .{
        .vertex_binding_description_count = binding_desc.len,
        .p_vertex_binding_descriptions = &binding_desc,
        .vertex_attribute_description_count = @intCast(attrs.len),
        .p_vertex_attribute_descriptions = attrs.ptr,
    };

    const input_assembly: vk.PipelineInputAssemblyStateCreateInfo = .{
        .topology = .triangle_strip,
        .primitive_restart_enable = .false,
    };

    const viewport_state: vk.PipelineViewportStateCreateInfo = .{
        .viewport_count = 1,
        .p_viewports = undefined,
        .scissor_count = 1,
        .p_scissors = undefined,
    };

    const rasterization: vk.PipelineRasterizationStateCreateInfo = .{
        .depth_clamp_enable = .false,
        .rasterizer_discard_enable = .false,
        .polygon_mode = .fill,
        .cull_mode = .{},
        .front_face = .clockwise,
        .depth_bias_enable = .false,
        .depth_bias_constant_factor = 0,
        .depth_bias_clamp = 0,
        .depth_bias_slope_factor = 0,
        .line_width = 1,
    };

    const multisample: vk.PipelineMultisampleStateCreateInfo = .{
        .rasterization_samples = .{ .@"1" = true },
        .sample_shading_enable = .false,
        .min_sample_shading = 1,
        .alpha_to_coverage_enable = .false,
        .alpha_to_one_enable = .false,
    };

    const blend_attachments: [1]vk.PipelineColorBlendAttachmentState = .{premultiplied_blend};
    const color_blend: vk.PipelineColorBlendStateCreateInfo = .{
        .logic_op_enable = .false,
        .logic_op = .copy,
        .attachment_count = blend_attachments.len,
        .p_attachments = &blend_attachments,
        .blend_constants = .{ 0, 0, 0, 0 },
    };

    const dyn_states: [2]vk.DynamicState = .{ .viewport, .scissor };
    const dynamic_state: vk.PipelineDynamicStateCreateInfo = .{
        .dynamic_state_count = dyn_states.len,
        .p_dynamic_states = &dyn_states,
    };

    const create_info: [1]vk.GraphicsPipelineCreateInfo = .{.{
        .stage_count = stages.len,
        .p_stages = &stages,
        .p_vertex_input_state = &vertex_input,
        .p_input_assembly_state = &input_assembly,
        .p_tessellation_state = null,
        .p_viewport_state = &viewport_state,
        .p_rasterization_state = &rasterization,
        .p_multisample_state = &multisample,
        .p_depth_stencil_state = null,
        .p_color_blend_state = &color_blend,
        .p_dynamic_state = &dynamic_state,
        .layout = layout,
        .render_pass = render_pass,
        .subpass = 0,
        .base_pipeline_handle = .null_handle,
        .base_pipeline_index = -1,
    }};

    var pipeline: [1]vk.Pipeline = undefined;
    _ = try dev.createGraphicsPipelines(.null_handle, &create_info, null, &pipeline);

    return .{ .handle = pipeline[0], .layout = layout, .descriptor_set_layout = descriptor_set_layout };
}

pub const Pipelines = struct {
    rect: Pipeline,
    glyph: Pipeline,
    image: Pipeline,

    pub fn create(dev: vk.DeviceProxy, render_pass: vk.RenderPass) !Pipelines {
        const rect = try createPipeline(dev, render_pass, shaders.rgl_vert, shaders.rgl_frag, &rect_attrs, @sizeOf(RectInstance), false);
        errdefer rect.destroy(dev);
        const glyph = try createPipeline(dev, render_pass, shaders.glyph_vert, shaders.glyph_frag, &glyph_attrs, @sizeOf(GlyphInstance), true);
        errdefer glyph.destroy(dev);
        const image = try createPipeline(dev, render_pass, shaders.image_vert, shaders.image_frag, &image_attrs, @sizeOf(ImageInstance), true);
        return .{ .rect = rect, .glyph = glyph, .image = image };
    }

    pub fn destroy(self: Pipelines, dev: vk.DeviceProxy) void {
        self.rect.destroy(dev);
        self.glyph.destroy(dev);
        self.image.destroy(dev);
    }
};

fn findMemoryType(mem_props: vk.PhysicalDeviceMemoryProperties, reqs: vk.MemoryRequirements, flags: vk.MemoryPropertyFlags) ?u32 {
    for (mem_props.memory_types[0..mem_props.memory_type_count], 0..) |mem_type, i| {
        const bit: u32 = @as(u32, 1) << @truncate(i);
        if (reqs.memory_type_bits & bit != 0 and mem_type.property_flags.contains(flags)) return @truncate(i);
    }
    return null;
}

/// A host-visible vertex buffer that grows (recreating itself) to fit
/// whatever's uploaded into it. One of these per instance-struct type, per
/// frame-in-flight, is the intended usage.
///
/// TODO: this recreates the whole buffer on every capacity increase rather
/// than a smarter growth policy, and offers no double-buffering of its
/// own -- fine for now since nothing calls this against a live swapchain
/// yet, but revisit once a real per-frame-in-flight caller exists.
pub fn InstanceBuffer(comptime T: type) type {
    return struct {
        const Self = @This();

        buf: vk.Buffer = .null_handle,
        mem: vk.DeviceMemory = .null_handle,
        capacity: usize = 0,

        pub fn destroy(self: *Self, dev: vk.DeviceProxy) void {
            if (self.buf != .null_handle) dev.destroyBuffer(self.buf, null);
            if (self.mem != .null_handle) dev.freeMemory(self.mem, null);
            self.* = .{};
        }

        /// Recreates the buffer if `items` doesn't fit the current one,
        /// then uploads via a map/copy/unmap (no persistent mapping).
        pub fn upload(self: *Self, dev: vk.DeviceProxy, mem_props: vk.PhysicalDeviceMemoryProperties, items: []const T) !void {
            if (items.len == 0) return;
            if (items.len > self.capacity) {
                self.destroy(dev);
                const size: u64 = @sizeOf(T) * items.len;
                self.buf = try dev.createBuffer(&.{
                    .size = size,
                    .usage = .{ .vertex_buffer = true },
                    .sharing_mode = .exclusive,
                }, null);
                errdefer dev.destroyBuffer(self.buf, null);

                const reqs = dev.getBufferMemoryRequirements(self.buf);
                const mem_type = findMemoryType(mem_props, reqs, .{ .host_visible = true, .host_coherent = true }) orelse return error.NoSuitableMemoryType;
                self.mem = try dev.allocateMemory(&.{ .allocation_size = reqs.size, .memory_type_index = mem_type }, null);
                try dev.bindBufferMemory(self.buf, self.mem, 0);
                self.capacity = items.len;
            }

            const mapped = try dev.mapMemory(self.mem, 0, vk.WHOLE_SIZE, .{});
            defer dev.unmapMemory(self.mem);
            const dst: [*]T = @ptrCast(@alignCast(mapped.?));
            @memcpy(dst[0..items.len], items);
        }
    };
}

/// Per-frame-in-flight GPU resources for the rect/roundrect/circle/line
/// pipeline -- the one primitive family that needs no texture/descriptor
/// set, so it's the only one fully wired up end to end (see the TODOs
/// below for glyph/image).
pub const RectFrameResources = struct {
    vbuf: InstanceBuffer(RectInstance) = .{},

    pub fn destroy(self: *RectFrameResources, dev: vk.DeviceProxy) void {
        self.vbuf.destroy(dev);
    }
};

/// Uploads `surface.rects` and records the draw for it into `cmd`, which
/// must already be recording inside an active render pass instance whose
/// viewport/scissor were set by the caller (both are dynamic state here,
/// same as shader_smoke_test.zig). No-op if there's nothing to draw.
///
/// Caller is responsible for calling `batch.clearBatches(surface)` once
/// this frame's submission is known to have been consumed (e.g. after the
/// corresponding fence signals) -- recording alone doesn't clear it, since
/// the instance data has to stay valid until the GPU has actually read it.
pub fn recordRectBatch(
    dev: vk.DeviceProxy,
    mem_props: vk.PhysicalDeviceMemoryProperties,
    cmd: vk.CommandBuffer,
    pipeline: Pipeline,
    resources: *RectFrameResources,
    surface: *const VkSurface,
    globals: Globals,
) !void {
    const instances = surface.rects.items;
    if (instances.len == 0) return;

    try resources.vbuf.upload(dev, mem_props, instances);

    dev.cmdBindPipeline(cmd, .graphics, pipeline.handle);
    dev.cmdPushConstants(cmd, pipeline.layout, .{ .vertex = true, .fragment = true }, 0, @sizeOf(Globals), &globals);
    const offsets = [_]vk.DeviceSize{0};
    dev.cmdBindVertexBuffers(cmd, 0, (&resources.vbuf.buf)[0..1], &offsets);
    dev.cmdDraw(cmd, 4, @intCast(instances.len), 0, 0);
}

/// TODO: glyph batches need a live glyph atlas texture (upload + image
/// view + sampler + descriptor set, rebuilt whenever the set of on-screen
/// glyphs changes) before this can bind one and draw -- see
/// shader_smoke_test.zig's buildStringAtlasViaSfnt/glyph_tex_* setup for
/// the reference version of that plumbing. Left unrecorded for now.
pub fn recordGlyphBatch(
    dev: vk.DeviceProxy,
    cmd: vk.CommandBuffer,
    pipeline: Pipeline,
    surface: *const VkSurface,
    globals: Globals,
) void {
    _ = dev;
    _ = cmd;
    _ = pipeline;
    _ = globals;
    if (surface.glyphs.items.len == 0) return;
    // TODO: bind an atlas descriptor set + instance buffer, then cmdDraw.
}

/// TODO: each unique source `Image` needs its own descriptor set (one
/// combined-image-sampler per drawImage call site, per ImageDraw's own doc
/// comment) -- allocating/caching/uploading those per-image textures isn't
/// implemented yet, so instances are grouped but never actually drawn.
pub fn recordImageBatch(
    dev: vk.DeviceProxy,
    cmd: vk.CommandBuffer,
    pipeline: Pipeline,
    surface: *const VkSurface,
    globals: Globals,
) void {
    _ = dev;
    _ = cmd;
    _ = pipeline;
    _ = globals;
    if (surface.images.items.len == 0) return;
    // TODO: for each unique `image` pointer, upload/bind its texture +
    // descriptor set, then cmdDraw its instances as one batch.
}
