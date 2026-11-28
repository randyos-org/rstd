//! A Vulkan device with no window, for the shader harness (issue #380): the
//! loader, an instance, the first device with a graphics queue, and helpers
//! for buffers, images, render targets, one-shot submissions and readback.
//! Lives in place: `init` fills a `Headless` the caller owns, as the
//! proxies point into it.

const std = @import("std");
const sysinfo = @import("builtin");
const vk = @import("vulkan");
const rstd = @import("rstd");

const gpu = rstd.rgl.vulkan.gpu;
const Allocator = std.mem.Allocator;

const Headless = @This();

const loader_name = switch (sysinfo.os.tag) {
    .windows => "vulkan-1.dll",
    .macos => "libvulkan.1.dylib",
    else => "libvulkan.so.1",
};

var get_instance_proc_addr: vk.PfnGetInstanceProcAddr = undefined;

fn getInstanceProcAddr(instance: ?vk.Instance, proc: [*:0]const u8) vk.PfnVoidFunction {
    return get_instance_proc_addr(instance, proc);
}

lib: rstd.dll.Dll,
base: vk.BaseWrapper,
inst_wrapper: vk.InstanceWrapper,
inst: vk.InstanceProxy,
phys: vk.PhysicalDevice,
props: vk.PhysicalDeviceProperties,
mem_props: vk.PhysicalDeviceMemoryProperties,
dev_wrapper: vk.DeviceWrapper,
dev: vk.DeviceProxy,
queue: vk.Queue,
pool: vk.CommandPool,

pub fn init(h: *Headless, gpa: Allocator) !void {
    h.lib = rstd.dll.Dll.openZ(loader_name) catch return error.NoVulkanLoader;
    errdefer h.lib.close();
    get_instance_proc_addr = h.lib.lookup(vk.PfnGetInstanceProcAddr, "vkGetInstanceProcAddr") orelse return error.NoVulkanLoader;
    h.base = vk.BaseWrapper.load(getInstanceProcAddr);

    const app: vk.ApplicationInfo = .{
        .p_application_name = "rstd-shader-run",
        .application_version = 0,
        .p_engine_name = "rgl",
        .engine_version = 0,
        .api_version = vk.API_VERSION_1_2.toU32(),
    };
    const inst = try h.base.createInstance(&.{ .p_application_info = &app }, null);
    h.inst_wrapper = vk.InstanceWrapper.load(inst, h.base.dispatch.vkGetInstanceProcAddr.?);
    h.inst = .init(inst, &h.inst_wrapper);
    errdefer h.inst.destroyInstance(null);

    // the first device with a graphics queue: lavapipe in CI
    const devices = try h.inst.enumeratePhysicalDevicesAlloc(gpa);
    defer gpa.free(devices);
    const family = for (devices) |d| {
        if (try graphicsFamily(h.inst, gpa, d)) |f| {
            h.phys = d;
            break f;
        }
    } else return error.NoVulkanDevice;
    h.props = h.inst.getPhysicalDeviceProperties(h.phys);
    h.mem_props = h.inst.getPhysicalDeviceMemoryProperties(h.phys);

    const priority = [_]f32{1};
    const queues = [_]vk.DeviceQueueCreateInfo{.{ .queue_family_index = family, .queue_count = 1, .p_queue_priorities = &priority }};
    const dev = try h.inst.createDevice(h.phys, &.{ .queue_create_info_count = queues.len, .p_queue_create_infos = &queues }, null);
    h.dev_wrapper = vk.DeviceWrapper.load(dev, h.inst_wrapper.dispatch.vkGetDeviceProcAddr.?);
    h.dev = .init(dev, &h.dev_wrapper);
    errdefer h.dev.destroyDevice(null);
    h.queue = h.dev.getDeviceQueue(family, 0);
    h.pool = try h.dev.createCommandPool(&.{ .queue_family_index = family }, null);
}

pub fn deinit(h: *Headless) void {
    h.dev.destroyCommandPool(h.pool, null);
    h.dev.destroyDevice(null);
    h.inst.destroyInstance(null);
    h.lib.close();
}

/// The device's name, as the driver reports it.
pub fn deviceName(h: *const Headless) []const u8 {
    return std.mem.sliceTo(&h.props.device_name, 0);
}

fn graphicsFamily(inst: vk.InstanceProxy, gpa: Allocator, d: vk.PhysicalDevice) !?u32 {
    const families = try inst.getPhysicalDeviceQueueFamilyPropertiesAlloc(d, gpa);
    defer gpa.free(families);
    for (families, 0..) |f, i| if (f.queue_flags.graphics) return @intCast(i);
    return null;
}

/// Records with `record(ctx, cmd)` into a one-time command buffer, submits
/// it and waits for the queue to drain.
pub fn submit(h: *const Headless, ctx: anytype, comptime record: fn (@TypeOf(ctx), vk.CommandBuffer) anyerror!void) !void {
    var cmd: [1]vk.CommandBuffer = undefined;
    try h.dev.allocateCommandBuffers(&.{ .command_pool = h.pool, .level = .primary, .command_buffer_count = 1 }, &cmd);
    defer h.dev.freeCommandBuffers(h.pool, &cmd);
    try h.dev.beginCommandBuffer(cmd[0], &.{ .flags = .{ .one_time_submit = true } });
    try record(ctx, cmd[0]);
    try h.dev.endCommandBuffer(cmd[0]);
    const info = [_]vk.SubmitInfo{.{ .command_buffer_count = 1, .p_command_buffers = &cmd }};
    try h.dev.queueSubmit(h.queue, &info, .null);
    try h.dev.queueWaitIdle(h.queue);
}

/// A buffer bound to its own memory.
pub const Buffer = struct {
    buf: vk.Buffer,
    mem: vk.DeviceMemory,
    size: usize,

    pub fn deinit(b: Buffer, h: *const Headless) void {
        h.dev.destroyBuffer(b.buf, null);
        h.dev.freeMemory(b.mem, null);
    }

    /// Its memory, mapped; `unmap` after.
    pub fn map(b: Buffer, h: *const Headless) ![]u8 {
        const p = try h.dev.mapMemory(b.mem, 0, vk.WHOLE_SIZE, .{});
        return @as([*]u8, @ptrCast(p.?))[0..b.size];
    }

    pub fn unmap(b: Buffer, h: *const Headless) void {
        h.dev.unmapMemory(b.mem);
    }
};

/// A host-visible, coherent buffer of `size` bytes.
pub fn buffer(h: *const Headless, size: usize, usage: vk.BufferUsageFlags) !Buffer {
    const buf = try h.dev.createBuffer(&.{ .size = size, .usage = usage, .sharing_mode = .exclusive }, null);
    errdefer h.dev.destroyBuffer(buf, null);
    const reqs = h.dev.getBufferMemoryRequirements(buf);
    const mem = try h.alloc(reqs, .{ .host_visible = true, .host_coherent = true });
    errdefer h.dev.freeMemory(mem, null);
    try h.dev.bindBufferMemory(buf, mem, 0);
    return .{ .buf = buf, .mem = mem, .size = size };
}

/// A host-visible buffer holding `bytes`.
pub fn bufferWith(h: *const Headless, bytes: []const u8, usage: vk.BufferUsageFlags) !Buffer {
    const b = try h.buffer(bytes.len, usage);
    errdefer b.deinit(h);
    @memcpy(try b.map(h), bytes);
    b.unmap(h);
    return b;
}

fn alloc(h: *const Headless, reqs: vk.MemoryRequirements, flags: vk.MemoryPropertyFlags) !vk.DeviceMemory {
    const kind = gpu.findMemoryType(h.mem_props, reqs, flags) orelse return error.NoSuitableMemoryType;
    return h.dev.allocateMemory(&.{ .allocation_size = reqs.size, .memory_type_index = kind }, null);
}

/// A 2D image in device memory and its view.
pub const Image = struct {
    image: vk.Image,
    mem: vk.DeviceMemory,
    view: vk.ImageView,
    format: vk.Format,
    width: u32,
    height: u32,

    pub fn deinit(i: Image, h: *const Headless) void {
        h.dev.destroyImageView(i.view, null);
        h.dev.destroyImage(i.image, null);
        h.dev.freeMemory(i.mem, null);
    }
};

const color_range: vk.ImageSubresourceRange = .{ .aspect_mask = .{ .color = true }, .base_mip_level = 0, .level_count = 1, .base_array_layer = 0, .layer_count = 1 };
const color_layers: vk.ImageSubresourceLayers = .{ .aspect_mask = .{ .color = true }, .mip_level = 0, .base_array_layer = 0, .layer_count = 1 };

pub fn image(h: *const Headless, format: vk.Format, width: u32, height: u32, usage: vk.ImageUsageFlags) !Image {
    const img = try h.dev.createImage(&.{
        .image_type = .@"2d",
        .format = format,
        .extent = .{ .width = width, .height = height, .depth = 1 },
        .mip_levels = 1,
        .array_layers = 1,
        .samples = .{ .@"1" = true },
        .tiling = .optimal,
        .usage = usage,
        .sharing_mode = .exclusive,
        .initial_layout = .undefined,
    }, null);
    errdefer h.dev.destroyImage(img, null);
    const mem = try h.alloc(h.dev.getImageMemoryRequirements(img), .{ .device_local = true });
    errdefer h.dev.freeMemory(mem, null);
    try h.dev.bindImageMemory(img, mem, 0);
    const view = try h.dev.createImageView(&.{
        .image = img,
        .view_type = .@"2d",
        .format = format,
        .components = .{ .r = .identity, .g = .identity, .b = .identity, .a = .identity },
        .subresource_range = color_range,
    }, null);
    return .{ .image = img, .mem = mem, .view = view, .format = format, .width = width, .height = height };
}

/// A color image to draw into: one subpass that clears it and leaves it
/// ready to copy out.
pub const Target = struct {
    color: Image,
    render_pass: vk.RenderPass,
    framebuffer: vk.Framebuffer,

    pub fn init(h: *const Headless, format: vk.Format, width: u32, height: u32) !Target {
        const color = try h.image(format, width, height, .{ .color_attachment = true, .transfer_src = true });
        errdefer color.deinit(h);
        const attachment = [_]vk.AttachmentDescription{.{
            .format = format,
            .samples = .{ .@"1" = true },
            .load_op = .clear,
            .store_op = .store,
            .stencil_load_op = .dont_care,
            .stencil_store_op = .dont_care,
            .initial_layout = .undefined,
            .final_layout = .transfer_src_optimal,
        }};
        const ref = [_]vk.AttachmentReference{.{ .attachment = 0, .layout = .color_attachment_optimal }};
        const subpass = [_]vk.SubpassDescription{.{ .pipeline_bind_point = .graphics, .color_attachment_count = 1, .p_color_attachments = &ref }};
        const render_pass = try h.dev.createRenderPass(&.{
            .attachment_count = 1,
            .p_attachments = &attachment,
            .subpass_count = 1,
            .p_subpasses = &subpass,
        }, null);
        errdefer h.dev.destroyRenderPass(render_pass, null);
        const views = [_]vk.ImageView{color.view};
        const framebuffer = try h.dev.createFramebuffer(&.{
            .render_pass = render_pass,
            .attachment_count = 1,
            .p_attachments = &views,
            .width = width,
            .height = height,
            .layers = 1,
        }, null);
        return .{ .color = color, .render_pass = render_pass, .framebuffer = framebuffer };
    }

    pub fn deinit(t: Target, h: *const Headless) void {
        h.dev.destroyFramebuffer(t.framebuffer, null);
        h.dev.destroyRenderPass(t.render_pass, null);
        t.color.deinit(h);
    }

    /// Begins the render pass, cleared to `clear`, with the whole target as
    /// the viewport and scissor.
    pub fn begin(t: Target, h: *const Headless, cmd: vk.CommandBuffer, clear: [4]f32) void {
        const area: vk.Rect2D = .{ .offset = .{ .x = 0, .y = 0 }, .extent = .{ .width = t.color.width, .height = t.color.height } };
        const values = [_]vk.ClearValue{.{ .color = .{ .float_32 = clear } }};
        h.dev.cmdBeginRenderPass(cmd, &.{
            .render_pass = t.render_pass,
            .framebuffer = t.framebuffer,
            .render_area = area,
            .clear_value_count = 1,
            .p_clear_values = &values,
        }, .@"inline");
        const viewport = [_]vk.Viewport{.{ .x = 0, .y = 0, .width = @floatFromInt(t.color.width), .height = @floatFromInt(t.color.height), .min_depth = 0, .max_depth = 1 }};
        h.dev.cmdSetViewport(cmd, 0, &viewport);
        h.dev.cmdSetScissor(cmd, 0, &.{area});
    }

    /// Ends the render pass and copies the image into `out`, rows packed.
    pub fn end(t: Target, h: *const Headless, cmd: vk.CommandBuffer, out: Buffer) void {
        h.dev.cmdEndRenderPass(cmd);
        const region = [_]vk.BufferImageCopy{.{
            .buffer_offset = 0,
            .buffer_row_length = 0,
            .buffer_image_height = 0,
            .image_subresource = color_layers,
            .image_offset = .{ .x = 0, .y = 0, .z = 0 },
            .image_extent = .{ .width = t.color.width, .height = t.color.height, .depth = 1 },
        }};
        h.dev.cmdCopyImageToBuffer(cmd, t.color.image, .transfer_src_optimal, out.buf, &region);
    }
};
