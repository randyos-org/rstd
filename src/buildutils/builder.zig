const rstdbuild = @import("__root__.zig");
const Build = rstdbuild.Build;
const Optimize = rstdbuild.Optimize;
const ResolvedTarget = rstdbuild.ResolvedTarget;

pub fn getOptimize(b: *Build) Optimize {
    return b.standardOptimizeOption(.{});
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
