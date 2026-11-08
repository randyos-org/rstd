const image = @import("image.zig");
pub const Image = image.Image;
pub const RgbaImage = image.RgbaImage;
pub const Surface = @import("Surface.zig");
pub const codecs = @import("codecs/__root__.zig");
pub const color = @import("color.zig");
pub const draw = @import("draw.zig");
pub const geometry = @import("geometry.zig");
pub const pack = @import("pack.zig");
pub const shaders = @import("shaders/__root__.zig");
pub const render = @import("render.zig");
pub const renderdoc = @import("renderdoc.zig");
pub const font = @import("font/__root__.zig");
pub const sw = @import("sw.zig");
pub const vulkan = @import("vulkan/__root__.zig");

const zigstd = @import("../zigstd.zig");

test {
    // Not here: shaders (its @embedFiles exist only with -Dvulkan=true),
    // vulkan.gpu (needs the vulkan-zig bindings, below) and renderdoc
    // (conditionally included).
    zigstd.testRefAllDecls(image);
    zigstd.testRefAllDecls(Surface);
    zigstd.testRefAllDecls(codecs);
    zigstd.testRefAllDecls(color);
    zigstd.testRefAllDecls(draw);
    zigstd.testRefAllDecls(geometry);
    zigstd.testRefAllDecls(pack);
    zigstd.testRefAllDecls(render);
    zigstd.testRefAllDecls(font);
    zigstd.testRefAllDecls(sw);
    zigstd.testRefAllDecls(vulkan);

    // with -Dvulkan=true, type-check gpu.zig too (no GPU needed)
    const rstd_options = @import("rstd_options");
    if (rstd_options.has_vulkan) {
        zigstd.testRefAllDecls(vulkan.gpu);
    }
}
