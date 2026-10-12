pub const zon = @import("../build.zig.zon");

const rstdbuild = @import("../src/buildutils/__root__.zig");
const blaskern = @import("blaskern.zig");
const qemutest = @import("qemutest.zig");
const wasmtest = @import("wasmtest.zig");
const Build = rstdbuild.Build;

pub fn build(b: *Build) !void {
    // standard build setup
    const optimize = rstdbuild.getOptimize(b);
    const target = rstdbuild.getTarget(b);
    const docs = rstdbuild.addDocs(b);

    const rstd_root = b.path("src/__root__.zig");
    const rstd_mod = b.addModule("rstd", .{
        .root_source_file = rstd_root,
    });
    const rstd_options = try rstdbuild.addBuildOptionsModule(b, @tagName(zon.name), zon.version, null);
    rstd_mod.addImport("rstd_options", try rstdbuild.optionsModule(rstd_options));
    // Debug builds link only the target's BLAS kernels, which is far faster
    // to compile; `-Dblas-cores=all` builds every core to test dispatch.
    const blas_cores = b.option(rstdbuild.BlasCores, "blas-cores", "BLAS kernel objects to build: the target's core only (debug default), all, or none") orelse
        @as(rstdbuild.BlasCores, if (optimize == .debug) .target else .all);
    blaskern.addBlasKernels(b, rstd_mod, target, blas_cores);
    rstdbuild.addWin32ToModIfWinTgt(b, target, rstd_mod);
    // only add RenderDoc to non-Vulkan, desktop-OS targets
    switch (target.result.os.tag) {
        .windows,
        .linux,
        .macos,
        => rstdbuild.addRenderDoc(b, target, optimize).addImportToModule(rstd_mod),
        else => {},
    }

    // rgl's Vulkan/SPIR-V shaders pull in no extra dependency (the
    // spirv32-vulkan target ships in the compiler), but rstd is a
    // general std lib supplement used by headless/CLI code too, so
    // compiling them only happens when asked for.
    const vulkan_opt = b.option(bool, "vulkan", "Compile rgl's Vulkan/SPIR-V shaders") orelse false;
    rstd_options.addOption(bool, "has_vulkan", vulkan_opt);
    if (vulkan_opt) {
        try rstdbuild.wireRglShaders(b, optimize, rstd_mod);
        rstdbuild.addVulkan(b).addImportToModule(rstd_mod);
    }

    // docs
    // ----
    // We want the rstd docs to be built separately from the other packages,
    // so we will make a separate compile step for them and attach the docs.
    // A Compile step needs its root module bound to a target, but rstd_mod is
    // deliberately target-less (consumers set their own via addImport) - so
    // clone it rather than redeclaring root_source_file/imports/options.
    const rstd_docs_mod = try rstdbuild.cloneModuleWithTarget(b, rstd_mod, target);
    // Docs only analyze the source; the kernel objects would make this object
    // link others, which the self-hosted COFF linker rejects.
    rstd_docs_mod.link_objects = .empty;
    // add the standalone rstd build step, just generates an unlinked object file
    const rstd_docs_obj = b.addObject(.{
        .name = "rstd",
        .root_module = rstd_docs_mod,
    });
    docs.addCompileStepDocs(b, rstd_docs_obj, "rstd");

    // tests
    // -----
    const rstd_test_mod = try rstdbuild.cloneModuleWithTarget(b, rstd_mod, target);
    const rstd_test_exe = b.addTest(.{
        .root_module = rstd_test_mod,
    });
    const run_rstd_test = b.addRunArtifact(rstd_test_exe);
    const test_step = b.step("rstd-test", "Run rstd's unit tests");
    test_step.dependOn(&run_rstd_test.step);

    // the same tests that need no OS, on emulated machines
    qemutest.addQemuTests(b, optimize);
    wasmtest.addWasmTests(b, optimize);
}
