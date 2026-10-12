const rstdbuild = @import("__root__.zig");
const Build = rstdbuild.Build;
const Optimize = rstdbuild.Optimize;
const ResolvedTarget = rstdbuild.ResolvedTarget;

pub fn getOptimize(b: *Build) Optimize {
    // zig 0.18-dev returns `?Optimize`; 0.17.0 returns `Optimize`. Gate on the
    // actual return type rather than `builtin.zig_version` so any dev build
    // works regardless of exactly when the change landed. The condition is
    // comptime-known, so the un-taken branch is never analyzed.
    const opt = b.standardOptimizeOption(.{});
    return if (@typeInfo(@TypeOf(opt)) == .optional) opt orelse .debug else opt;
}

pub fn getTarget(b: *Build) ResolvedTarget {
    return b.standardTargetOptions(.{});
}
