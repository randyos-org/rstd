const rstdbuild = @import("__root__.zig");
const Build = rstdbuild.Build;
const OptimizeMode = rstdbuild.OptimizeMode;
const ResolvedTarget = rstdbuild.ResolvedTarget;

pub fn getOptimize(b: *Build) OptimizeMode {
    return b.standardOptimizeOption(.{});
}

pub fn getTarget(b: *Build) ResolvedTarget {
    return b.standardTargetOptions(.{});
}
