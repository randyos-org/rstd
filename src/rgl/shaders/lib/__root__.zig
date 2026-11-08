pub const compile = @import("compile.zig");
pub const coords = @import("coords.zig");
pub const math = @import("math.zig");
pub const sdf = @import("sdf.zig");

pub const types = @import("../types.zig");
const layout = @import("../layout.zig");
const V2f32 = types.V2f32;

pub const Globals = layout.Globals(V2f32);
pub const globals = @extern(*addrspace(.push_constant) const Globals, .{
    .name = "globals",
});
