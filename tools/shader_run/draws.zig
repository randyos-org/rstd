//! rgl's three pipelines (`shapes`, `glyph`, `image`) on the device, checked
//! against the software renderer (issue #380): one scene drawn through
//! `rgl.vulkan` and recorded with `gpu.zig`, and through `rgl.sw`, then
//! compared pixel by pixel. Text and images sample texels 1:1 and must match
//! everywhere; shapes, whose edges the shaders antialias and `sw` does not,
//! must match wherever `sw`'s image is uniform `margin` pixels around.

const std = @import("std");
const vk = @import("vulkan");
const rstd = @import("rstd");
const Headless = @import("Headless.zig");

const rgl = rstd.rgl;
const gpu = rgl.vulkan.gpu;
const vga = rgl.font.vga;
const RGBA = rgl.color.RGBA;
const RgbaImage = rgl.RgbaImage;
const Rect = rgl.geometry.Rect;
const Allocator = std.mem.Allocator;

const width = 96;
const height = 64;
/// Pixels a shape's edge may reach past `sw`'s: the shaders' antialiasing.
const margin = 2;
const background: RGBA = .{ .r = 16, .g = 32, .b = 48, .a = 255 };

/// The scene's draws, none overlapping, so the GPU's per-pipeline order and
/// `sw`'s call order give one image. Colors are opaque: blended over the
/// background, they come out as they went in.
fn scene(comptime S: type, s: *S, image: *const RgbaImage) !void {
    const m = if (S == rgl.sw.SwSurface) rgl.sw else rgl.vulkan;
    const call = struct {
        fn f(r: anytype) !void {
            if (@typeInfo(@TypeOf(r)) == .error_union) _ = try r;
        }
    }.f;
    try call(m.fillRect(s, .{ .x = 4, .y = 4, .w = 20, .h = 12 }, .{ .r = 255, .g = 0, .b = 0, .a = 255 }));
    try call(m.strokeRect(s, .{ .x = 28, .y = 4, .w = 20, .h = 14 }, 3, .{ .r = 255, .g = 255, .b = 255, .a = 255 }));
    try call(m.fillRoundRect(s, .{ .x = 52, .y = 4, .w = 24, .h = 16 }, 5, .{ .r = 0, .g = 200, .b = 0, .a = 255 }));
    try call(m.fillCircle(s, .{ .x = 86, .y = 12 }, 8, .{ .r = 40, .g = 80, .b = 255, .a = 255 }));
    try call(m.drawLine(s, .{ .x = 4, .y = 28 }, .{ .x = 40, .y = 40 }, 3, .{ .r = 255, .g = 220, .b = 0, .a = 255 }));
    try call(m.drawText(s, .{ .x = 48, .y = 24 }, "Hi rgl", .{ .r = 0, .g = 230, .b = 230, .a = 255 }));
    try call(m.drawImage(s, .{ .x = 8, .y = 44 }, image));
}

/// Where the text and image land: compared everywhere.
const exact = [_]Rect{
    .{ .x = 48, .y = 24, .w = 6 * vga.width, .h = vga.height },
    .{ .x = 8, .y = 44, .w = image_size, .h = image_size },
};
const image_size = 12;

/// A test image of opaque, distinct texels.
fn testImage(gpa: Allocator) !RgbaImage {
    const img = try RgbaImage.init(gpa, image_size, image_size, null);
    for (0..image_size) |y| for (0..image_size) |x| {
        img.data[y * image_size + x] = .{ .r = @intCast(x * 20), .g = @intCast(y * 20), .b = @intCast(255 - (x + y) * 9), .a = 255 };
    };
    return img;
}

/// The VGA font as one row of 256 glyphs, one coverage byte per texel:
/// the atlas `rgl.vulkan.drawText` lays text out on.
fn vgaAtlas(gpa: Allocator) ![]u8 {
    const w = vga.glyph_count * vga.width;
    const atlas = try gpa.alloc(u8, w * vga.height);
    for (0..vga.glyph_count) |cp| for (vga.glyph(@intCast(cp)), 0..) |bits, row| {
        for (0..vga.width) |col| {
            const on = (bits >> @intCast(7 - col)) & 1 != 0;
            atlas[row * w + cp * vga.width + col] = if (on) 255 else 0;
        }
    };
    return atlas;
}

/// Draws the scene both ways; the count of pixels that differ.
pub fn run(h: *Headless, gpa: Allocator) !u32 {
    var image = try testImage(gpa);
    defer image.deinit(gpa);

    var cpu = try rgl.sw.SwSurface.init(gpa, width, height, .{ .fill = background });
    defer cpu.deinit(gpa);
    try scene(rgl.sw.SwSurface, &cpu, &image);

    var surface = try rgl.vulkan.VkSurface.init(gpa, width, height, .{});
    defer surface.deinit(gpa);
    try scene(rgl.vulkan.VkSurface, &surface, &image);

    const target = try Headless.Target.init(h, .r8g8b8a8_unorm, width, height);
    defer target.deinit(h);
    const pipelines = try gpu.Pipelines.create(h.dev, target.render_pass);
    defer pipelines.deinit(h.dev);
    var res: gpu.FrameResources = .{};
    defer res.deinit(h.dev);

    const sampler = try h.dev.createSampler(&.{
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
    defer h.dev.destroySampler(sampler, null);

    const atlas_bytes = try vgaAtlas(gpa);
    defer gpa.free(atlas_bytes);
    const atlas = try h.texture(.r8_unorm, vga.glyph_count * vga.width, vga.height, atlas_bytes);
    defer atlas.deinit(h);
    const img_tex = try h.texture(.r8g8b8a8_unorm, image_size, image_size, std.mem.sliceAsBytes(image.data));
    defer img_tex.deinit(h);

    const sizes = [_]vk.DescriptorPoolSize{.{ .type = .combined_image_sampler, .descriptor_count = 2 }};
    const pool = try h.dev.createDescriptorPool(&.{ .max_sets = 2, .pool_size_count = sizes.len, .p_pool_sizes = &sizes }, null);
    defer h.dev.destroyDescriptorPool(pool, null);
    const atlas_set = try textureSet(h, pool, pipelines.glyph, sampler, atlas);
    const image_set = try textureSet(h, pool, pipelines.image, sampler, img_tex);

    // the image pipeline's instances, one draw per image (`rgl.vulkan`
    // groups them by source image)
    if (surface.images.items.len != 1) return error.UnexpectedImageDraws;
    const image_instance = try h.bufferWith(std.mem.asBytes(&surface.images.items[0].instance), .{ .vertex_buffer = true });
    defer image_instance.deinit(h);

    const out = try h.buffer(width * height * 4, .{ .transfer_dst = true });
    defer out.deinit(h);
    try h.submit(Frame{
        .h = h,
        .target = target,
        .pipelines = pipelines,
        .res = &res,
        .surface = &surface,
        .atlas_set = atlas_set,
        .image_set = image_set,
        .image_instance = image_instance.buf,
        .out = out,
    }, Frame.record);
    const pixels = try out.map(h);
    defer out.unmap(h);
    return compare(pixels, cpu.target);
}

fn textureSet(h: *const Headless, pool: vk.DescriptorPool, pipeline: gpu.Pipeline, sampler: vk.Sampler, tex: Headless.Image) !vk.DescriptorSet {
    var set: [1]vk.DescriptorSet = undefined;
    try h.dev.allocateDescriptorSets(&.{ .descriptor_pool = pool, .descriptor_set_count = 1, .p_set_layouts = &.{pipeline.descriptor_set_layout} }, &set);
    const info = [_]vk.DescriptorImageInfo{.{ .sampler = sampler, .image_view = tex.view, .image_layout = .shader_read_only_optimal }};
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
    h.dev.updateDescriptorSets(&write, null);
    return set[0];
}

const Frame = struct {
    h: *const Headless,
    target: Headless.Target,
    pipelines: gpu.Pipelines,
    res: *gpu.FrameResources,
    surface: *const rgl.vulkan.VkSurface,
    atlas_set: vk.DescriptorSet,
    image_set: vk.DescriptorSet,
    image_instance: vk.Buffer,
    out: Headless.Buffer,

    fn record(f: Frame, cmd: vk.CommandBuffer) !void {
        const h = f.h;
        const bg = background.toRGBAf32();
        f.target.begin(h, cmd, .{ bg.r, bg.g, bg.b, bg.a });
        const globals: gpu.Globals = .{
            .screen_size = .{ width, height },
            .atlas_size = .{ vga.glyph_count * vga.width, vga.height },
        };
        try gpu.recordBatch(h.dev, h.mem_props, cmd, f.pipelines.rect, &f.res.rect, f.surface.batch(rgl.vulkan.rect_pipeline_id), globals, null);
        try gpu.recordBatch(h.dev, h.mem_props, cmd, f.pipelines.glyph, &f.res.glyph, f.surface.batch(rgl.vulkan.glyph_pipeline_id), globals, f.atlas_set);
        const image = f.pipelines.image;
        h.dev.cmdBindPipeline(cmd, .graphics, image.handle);
        h.dev.cmdBindDescriptorSets(cmd, .graphics, image.layout, 0, &.{f.image_set}, null);
        h.dev.cmdPushConstants(cmd, image.layout, .{ .vertex = true, .fragment = true }, 0, @sizeOf(gpu.Globals), &globals);
        h.dev.cmdBindVertexBuffers(cmd, 0, &.{f.image_instance}, &.{0});
        h.dev.cmdDraw(cmd, 4, 1, 0, 0);
        f.target.end(h, cmd, f.out);
    }
};

fn inside(r: Rect, x: usize, y: usize) bool {
    const ix: i32 = @intCast(x);
    const iy: i32 = @intCast(y);
    return ix >= r.x and iy >= r.y and ix < r.x + @as(i32, @intCast(r.w)) and iy < r.y + @as(i32, @intCast(r.h));
}

/// Whether `sw`'s pixels within `margin` of `(x, y)` are all one color.
fn uniform(img: RgbaImage, x: usize, y: usize) bool {
    const c = img.data[y * width + x];
    if (x < margin or y < margin or x + margin >= width or y + margin >= height) return false;
    for (y - margin..y + margin + 1) |yy| for (x - margin..x + margin + 1) |xx| {
        if (!img.data[yy * width + xx].equals(c)) return false;
    };
    return true;
}

/// The pixels where `got` (the device's, RGBA8 rows) differs from `want` by
/// more than one step in a channel, of those compared.
fn compare(got: []const u8, want: RgbaImage) u32 {
    var failures: u32 = 0;
    var compared: u32 = 0;
    for (0..height) |y| for (0..width) |x| {
        const is_exact = for (exact) |r| {
            if (inside(r, x, y)) break true;
        } else false;
        if (!is_exact and !uniform(want, x, y)) continue;
        compared += 1;
        const w = want.data[y * width + x];
        const g = got[(y * width + x) * 4 ..][0..4];
        const diff = @max(@max(absDiff(g[0], w.r), absDiff(g[1], w.g)), @max(absDiff(g[2], w.b), absDiff(g[3], w.a)));
        if (diff <= 1) continue;
        failures += 1;
        if (failures <= 20) std.debug.print("draw ({d},{d}): got {d} {d} {d} {d}, want {d} {d} {d} {d}\n", .{ x, y, g[0], g[1], g[2], g[3], w.r, w.g, w.b, w.a });
    };
    std.debug.print("draws: {d} of {d} pixels compared\n", .{ compared, width * height });
    return failures;
}

fn absDiff(a: u8, b: u8) u8 {
    return if (a > b) a - b else b - a;
}
