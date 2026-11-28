//! rgl's three pipelines (`shapes`, `glyph`, `image`) on the device, checked
//! against the software renderer (issue #380): one scene drawn through
//! `rgl.vulkan` and recorded with `gpu.SurfaceDraws`, and through `rgl.sw`, then
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
    var draws: gpu.SurfaceDraws = .{};
    defer draws.deinit(gpa, h.dev);
    const ctx: gpu.Texture.Context = .{ .dev = h.dev, .mem_props = h.mem_props, .queue = h.queue, .pool = h.pool };
    try draws.prepare(gpa, ctx, pipelines, &surface);

    const out = try h.buffer(width * height * 4, .{ .transfer_dst = true });
    defer out.deinit(h);
    try h.submit(Frame{
        .h = h,
        .target = target,
        .pipelines = pipelines,
        .res = &res,
        .draws = &draws,
        .surface = &surface,
        .out = out,
    }, Frame.record);
    const pixels = try out.map(h);
    defer out.unmap(h);
    return compare(pixels, cpu.target);
}

const Frame = struct {
    h: *const Headless,
    target: Headless.Target,
    pipelines: gpu.Pipelines,
    res: *gpu.FrameResources,
    draws: *const gpu.SurfaceDraws,
    surface: *const rgl.vulkan.VkSurface,
    out: Headless.Buffer,

    fn record(f: Frame, cmd: vk.CommandBuffer) !void {
        const h = f.h;
        const bg = background.toRGBAf32();
        f.target.begin(h, cmd, .{ bg.r, bg.g, bg.b, bg.a });
        try f.draws.record(h.dev, h.mem_props, cmd, f.pipelines, f.res, f.surface, .{ width, height });
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
