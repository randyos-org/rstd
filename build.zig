//! Build script entry point
pub const build = @import("build/__root__.zig").build;

pub const buildutils = @import("src/buildutils.zig");
pub const machine = @import("src/machine.zig");
pub const io = @import("src/io.zig");
