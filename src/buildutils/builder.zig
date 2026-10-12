const rstdbuild = @import("__root__.zig");
const Build = rstdbuild.Build;
const Optimize = rstdbuild.Optimize;
const ResolvedTarget = rstdbuild.ResolvedTarget;

pub fn getOptimize(b: *Build) Optimize {
    return b.standardOptimizeOption(.{}) orelse .debug;
}

pub fn getTarget(b: *Build) ResolvedTarget {
    return b.standardTargetOptions(.{});
}
