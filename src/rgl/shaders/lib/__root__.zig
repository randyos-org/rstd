pub const compile = @import("compile.zig");
pub const coords = @import("coords.zig");
pub const math = @import("math.zig");
pub const sdf = @import("sdf.zig");
pub const textures = @import("textures.zig");
pub const types = @import("types.zig");

const V2f32 = types.V2f32;

/// Per-draw globals, bound as a push constant by every pipeline here.
pub const Globals = extern struct {
    screen_size: V2f32,
    /// only meaningful to `glyph.zig` (normalizes atlas texture coordinates)
    atlas_size: V2f32,
};
pub const globals = @extern(*addrspace(.push_constant) const Globals, .{
    .name = "globals",
});
