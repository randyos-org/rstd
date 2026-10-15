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

/// A command run from the build root on the user's console, for one that may
/// ask questions, open a window, or keep running (a VM, a server).
pub fn console(b: *Build, argv: []const []const u8) *Build.Step.Run {
    const run = b.addSystemCommand(argv);
    run.setCwd(b.path("."));
    run.stdio = .inherit;
    return run;
}
