const rstdbuild = @import("__root__.zig");
const Build = rstdbuild.Build;
const ResolvedTarget = rstdbuild.ResolvedTarget;
const Module = rstdbuild.Module;
const OptimizeMode = rstdbuild.OptimizeMode;
const DepModule = @import("deps.zig").DepModule;

pub fn addRenderDoc(b: *Build, target: ResolvedTarget, optimize: OptimizeMode) DepModule {
    const tc = b.addTranslateC(.{
        .root_source_file = b.dependency("renderdoc", .{}).path("renderdoc/api/app/renderdoc_app.h"),
        .target = target,
        .optimize = optimize,
        .link_libc = true,
    });
    // renderdoc_app.h assumes `bool` is already defined
    // (true under C++, but not plain C); supply what stdbool.h would.
    tc.defineCMacro("bool", "_Bool");
    tc.defineCMacro("true", "1");
    tc.defineCMacro("false", "0");
    return .{
        .mod = tc.createModule(),
        .name = "renderdoc",
    };
}
