const zigstd = @import("../zigstd.zig");

pub const Image = @import("Image.zig");
pub const Surface = @import("Surface.zig");
pub const codecs = @import("codecs/__root__.zig");
pub const color = @import("color.zig");
pub const draw = @import("draw.zig");
pub const geometry = @import("geometry.zig");
pub const shaders = @import("shaders/__root__.zig");
pub const render = @import("render.zig");
pub const renderdoc = @import("renderdoc.zig");
pub const font = @import("font/__root__.zig");
pub const sw = @import("sw/__root__.zig");
pub const vulkan = @import("vulkan/__root__.zig");

test {
    // shaders excluded: its @embedFiles only resolve when build.zig
    // compiled rgl's shaders (-Dvulkan=true), which a plain `zig build
    // rstd-test` doesn't do by default. vulkan.gpu (real vk-zig-backed
    // pipeline/command-buffer code) is excluded for the same reason --
    // it's only ever analyzed by something that actually imports it.
    //
    // renderdoc is also excluded since it is only conditionally included
    zigstd.testRefAllDecls(Image);
    zigstd.testRefAllDecls(Surface);
    zigstd.testRefAllDecls(codecs);
    zigstd.testRefAllDecls(color);
    zigstd.testRefAllDecls(draw);
    zigstd.testRefAllDecls(geometry);
    zigstd.testRefAllDecls(render);
    zigstd.testRefAllDecls(font);
    zigstd.testRefAllDecls(sw);
    zigstd.testRefAllDecls(vulkan);

    // gpu.zig needs the real vulkan-zig bindings, only wired up with
    // -Dvulkan=true (see src/buildutils/vulkan.zig); refAllDecls it too
    // when that's on, so its pipeline/command-buffer code at least
    // type-checks even without a physical GPU to run it against.
    const rstd_options = @import("rstd_options");
    if (rstd_options.has_vulkan) {
        zigstd.testRefAllDecls(vulkan.gpu);
    }
}
