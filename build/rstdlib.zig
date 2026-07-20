const std = @import("std");
const buildroot = @import("__root__.zig");
const rstdbuild = buildroot.rstd.buildutils;

const Build = rstdbuild.Build;
const RunStep = rstdbuild.RunStep;
const ResolvedTarget = rstdbuild.ResolvedTarget;
const Module = rstdbuild.Module;
const OptimizeMode = rstdbuild.OptimizeMode;
const BuildOptions = rstdbuild.BuildOptions;
const Docs = rstdbuild.Docs;

pub fn addRstd(
    b: *Build,
    docs: Docs,
    options: *BuildOptions,
    target: ?ResolvedTarget,
) *Module {
    const rstd_root = b.path("src/rstd/__root__.zig");
    const rstd_mod = b.createModule(.{
        .root_source_file = rstd_root,
    });
    rstdbuild.addBuildOptsModToModule(options, rstd_mod);

    // We want the rstd docs to be built separately from the other packages,
    // so we will make a separate compile step for them and attach the docs.
    const rstd_docs_mod = b.createModule(.{
        .root_source_file = rstd_root,
        .target = target,
    });
    rstdbuild.addBuildOptsModToModule(options, rstd_docs_mod);
    // add the standalone rstd build step, just generates an unlinked object file
    const rstd_docs_obj = b.addObject(.{
        .name = "rstd",
        .root_module = rstd_docs_mod,
    });
    docs.addCompileStepDocs(b, rstd_docs_obj, "rstd");

    return rstd_mod;
}
