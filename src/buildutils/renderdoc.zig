const rstdbuild = @import("__root__.zig");
const Build = rstdbuild.Build;
const ResolvedTarget = rstdbuild.ResolvedTarget;
const Module = rstdbuild.Module;
const Optimize = rstdbuild.Optimize;
const DepModule = @import("deps.zig").DepModule;

/// RenderDoc's in-app API header translated to a module, or null while its
/// lazy dependencies (renderdoc, translate_c) are still being fetched; the
/// build reruns once they are. `b` and `build_zig` must be rstd's own
/// builder and build.zig, since those dependencies are rstd's.
pub fn addRenderDoc(b: *Build, comptime build_zig: type, target: ResolvedTarget, optimize: Optimize) ?DepModule {
    // Ask for everything before bailing, so one fetch round gets it all.
    const renderdoc_dep = b.dependencyLazy("renderdoc", .{}) catch return null;
    // const translate_c_dep = b.lazyDependency("translate_c", .{}) orelse return null;
    // // std.Build.Step.TranslateC is deprecated (0.17.0) for this package.
    // const translate_c = b.lazyImport(build_zig, "translate_c") orelse return null;
    _ = build_zig; // only the translate_c package needs it

    // const tc: translate_c.Translator = .init(translate_c_dep orelse return null, .{
    const tc = b.addTranslateC(.{
        // .c_source_file
        .root_source_file = renderdoc_dep.path("renderdoc/api/app/renderdoc_app.h"),
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
        .mod = tc.createModule(), // tc.mod
        .name = "renderdoc",
    };
}
