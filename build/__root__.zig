pub const zon = @import("../build.zig.zon");

const std = @import("std");
const rstdbuild = @import("../src/buildutils/__root__.zig");
const blaskern = @import("blaskern.zig");
const qemutest = @import("qemutest.zig");
const wasmtest = @import("wasmtest.zig");
const Build = rstdbuild.Build;
const Module = rstdbuild.Module;
const Compile = rstdbuild.CompileStep;
const ResolvedTarget = rstdbuild.ResolvedTarget;
const Optimize = rstdbuild.Optimize;

/// rstd's module for `target`: its options, BLAS kernels and platform imports.
/// `name` registers it for consumers; null keeps it private to this build (a
/// second target, as `qemutest.zig`'s macOS VM uses).
pub fn addRstdModule(
    b: *Build,
    name: ?[]const u8,
    target: ResolvedTarget,
    optimize: Optimize,
    blas_cores: rstdbuild.BlasCores,
    vulkan: bool,
    /// RenderDoc's in-app API for `rgl.renderdoc`; desktop targets only.
    renderdoc: bool,
) !*Module {
    const root = b.path("src/__root__.zig");
    const mod = if (name) |n|
        b.addModule(n, .{ .root_source_file = root })
    else
        b.createModule(.{ .root_source_file = root, .target = target, .optimize = optimize });
    const has_renderdoc = renderdoc and switch (target.result.os.tag) {
        .windows, .linux, .macos => true,
        else => false,
    };
    const options = try rstdbuild.addBuildOptionsModule(b, @tagName(zon.name), zon.version, null);
    options.addOption(bool, "has_vulkan", vulkan);
    options.addOption(bool, "has_renderdoc", has_renderdoc);
    mod.addImport("rstd_options", try rstdbuild.optionsModule(options));
    blaskern.addBlasKernels(b, mod, target, blas_cores);
    rstdbuild.addWin32ToModIfWinTgt(b, target, mod);
    if (has_renderdoc) {
        // null until its lazy dependencies are fetched; the build then reruns
        if (rstdbuild.addRenderDoc(b, @import("../build.zig"), target, optimize)) |rd| rd.addImportToModule(mod);
    }
    return mod;
}

/// rstd's unit tests for `target`, from a module of `addRstdModule`.
pub fn addRstdTest(b: *Build, mod: *Module, target: ResolvedTarget, filters: []const []const u8) !*Compile {
    return b.addTest(.{ .root_module = try rstdbuild.cloneModuleWithTarget(b, mod, target), .filters = filters });
}

pub fn build(b: *Build) !void {
    // standard build setup
    const optimize = rstdbuild.getOptimize(b);
    const target = rstdbuild.getTarget(b);
    const docs = rstdbuild.addDocs(b);

    // Debug builds link only the target's BLAS kernels, which is far faster
    // to compile; `-Dblas-cores=all` builds every core to test dispatch.
    const blas_cores = b.option(rstdbuild.BlasCores, "blas-cores", "BLAS kernel objects to build: the target's core only (debug default), all, or none") orelse
        @as(rstdbuild.BlasCores, if (optimize == .debug) .target else .all);
    // rgl's Vulkan/SPIR-V shaders pull in no extra dependency (the
    // spirv32-vulkan target ships in the compiler), but rstd is a
    // general std lib supplement used by headless/CLI code too, so
    // compiling them only happens when asked for.
    const vulkan_opt = b.option(bool, "vulkan", "Compile rgl's Vulkan/SPIR-V shaders") orelse false;
    const renderdoc_opt = b.option(bool, "renderdoc", "Link RenderDoc's in-app capture API (fetches renderdoc and translate-c)") orelse false;
    const rstd_mod = try addRstdModule(b, "rstd", target, optimize, blas_cores, vulkan_opt, renderdoc_opt);
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
    const test_filter = b.option([]const []const u8, "rstd-test-filter", "Only run rstd tests whose name contains this (repeatable)") orelse &.{};
    const rstd_test_exe = try addRstdTest(b, rstd_mod, target, test_filter);
    const run_rstd_test = b.addRunArtifact(rstd_test_exe);
    const test_step = b.step("rstd-test", "Run rstd's unit tests");
    test_step.dependOn(&run_rstd_test.step);
    // the same tests on real machines over SSH, for targets no emulator covers
    try @import("sshtest.zig").addSshTests(b, optimize, test_filter);

    // Mach-O, which no machine here runs: building is what catches assembly
    // that other object formats accept (`.align`), so build every core's
    // kernels for both architectures. `qemu-test-x86_64-macos` runs them.
    const macho_step = b.step("rstd-check-macho", "Build rstd's unit tests for macOS (x86_64, aarch64) without running them");
    for ([_]std.Target.Cpu.Arch{ .x86_64, .aarch64 }) |arch| {
        const macho = b.resolveTargetQuery(.{ .cpu_arch = arch, .os_tag = .macos });
        const macho_mod = try addRstdModule(b, null, macho, optimize, .all, false, false);
        macho_step.dependOn(&(try addRstdTest(b, macho_mod, macho, &.{})).step);
    }

    // UEFI, which no test runner supports: the OS-free tests (the qemu rigs'
    // root: simd and BLAS, every interface at every element type) built, not
    // run, through a runner that only links them
    const uefi_step = b.step("rstd-check-uefi", "Build rstd's OS-free tests (simd, BLAS) for UEFI (x86_64, aarch64) without running them");
    for ([_]std.Target.Cpu.Arch{ .x86_64, .aarch64 }) |arch| {
        const uefi = b.resolveTargetQuery(.{ .cpu_arch = arch, .os_tag = .uefi });
        const uefi_mod = b.createModule(.{ .root_source_file = b.path("src/_test_qemu.zig"), .target = uefi, .optimize = optimize });
        blaskern.addBlasKernels(b, uefi_mod, uefi, .all);
        uefi_step.dependOn(&b.addTest(.{
            .name = b.fmt("rstd-check-uefi-{t}", .{arch}),
            .root_module = uefi_mod,
            .test_runner = .{ .path = b.path("test/uefi/runner.zig"), .mode = .simple },
        }).step);
    }

    // the same tests that need no OS, on emulated machines
    try qemutest.addQemuTests(b, optimize);
    wasmtest.addWasmTests(b, optimize);
}
