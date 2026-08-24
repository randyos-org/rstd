//! Build script entry point
pub const buildroot = @import("build/__root__.zig");
pub const build = buildroot.build;
pub const zon = buildroot.zon;

pub const buildutils = @import("src/buildutils.zig");
pub const machine = @import("src/machine.zig");
pub const io = @import("src/io.zig");
