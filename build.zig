//! Build script entry point
pub const buildroot = @import("build/__root__.zig");
pub const build = buildroot.build;

pub const buildutils = @import("src/buildutils/__root__.zig");
