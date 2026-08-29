const rstdbuild = @import("__root__.zig");
const docs_mod = @import("docs.zig");
const module = @import("module.zig");

const Build = rstdbuild.Build;
const ResolvedTarget = rstdbuild.ResolvedTarget;
const Docs = docs_mod.Docs;
const DepModule = @import("deps.zig").DepModule;

pub fn addRstd(
    b: *Build,
    target: ?ResolvedTarget,
    docs: ?Docs,
    vulkan: bool,
) !DepModule {
    const rstd_dep = b.dependency("rstd", .{
        .target = target,
        .vulkan = vulkan,
    });
    // rstd_mod already has its own "build_options" import, wired up by
    // rstd's own build.zig when it configured itself; don't re-add one
    // here. `Module.addOptions`/`Options.createModule()` isn't memoized, so
    // doing that would stamp out a second Module rooted at the exact same
    // generated file the caller's own module wraps too, which the compiler
    // now rejects ("file exists in modules 'build_options' and
    // 'build_options0'").
    const rstd_mod = rstd_dep.module("rstd");

    if (docs) |d| {
        // The registered "rstd" module has no target (consumers set their
        // own via addImport), so clone it with one bound rather than
        // redeclaring root_source_file/imports/options.
        const rstd_docs_mod = try module.cloneModuleWithTarget(b, rstd_mod, target);
        const rstd_docs_obj = b.addObject(.{
            .name = "rstd",
            .root_module = rstd_docs_mod,
        });
        d.addCompileStepDocs(b, rstd_docs_obj, "rstd");
    }

    return .{
        .mod = rstd_mod,
        .name = "rstd",
        .dep = rstd_dep,
    };
}
