const std = @import("std");

pub const Image = @import("Image.zig");
pub const codecs = @import("codecs/__root__.zig");
pub const color = @import("color.zig");
pub const draw = @import("draw.zig");
pub const shaders = @import("shaders/__root__.zig");
pub const render = @import("render.zig");
pub const renderdoc = @import("renderdoc.zig");
pub const font = @import("font/__root__.zig");

test {
    // shaders excluded: its @embedFiles only resolve when build.zig
    // compiled rgl's shaders (-Dvulkan=true), which a plain `zig build
    // rstd-test` doesn't do by default.
    //
    // renderdoc is also excluded since it is only conditionally included
    std.testing.refAllDecls(Image);
    std.testing.refAllDecls(codecs);
    std.testing.refAllDecls(color);
    std.testing.refAllDecls(draw);
    std.testing.refAllDecls(render);
    std.testing.refAllDecls(font);
}
