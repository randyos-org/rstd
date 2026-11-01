pub const zon = @import("../build.zig.zon");

const std = @import("std");
const sysinfo = @import("builtin");
const rstdbuild = @import("../src/buildutils/__root__.zig");
const build_root = @import("../build.zig");
const blaskern = @import("blaskern.zig");
const qemutest = @import("qemutest.zig");
const wasmtest = @import("wasmtest.zig");
const sshtest = @import("sshtest.zig");
const shader_probe = @import("../src/rgl/shaders/probe.zig");
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
    // Stock Zig rejects `@Matrix` even in a file it never analyzes, so the
    // file spelling it is imported only by compilers that have it.
    const matrix = if (@hasField(std.lang.Type, "matrix")) "native" else "fallback";
    mod.addImport("rstd_matrix", b.createModule(.{ .root_source_file = b.path("src/linalg/_matrix/" ++ matrix ++ ".zig") }));
    try blaskern.addBlasKernels(b, mod, target, blas_cores, null);
    rstdbuild.addWin32ToModIfWinTgt(b, target, mod);
    if (has_renderdoc) {
        // null until its lazy dependencies are fetched; the build then reruns
        if (rstdbuild.addRenderDoc(b, build_root, target, optimize)) |rd| rd.addImportToModule(mod);
    }
    return mod;
}

/// `-Drstd-tsan`, which `build` reads before any test is added.
pub var tsan = false;

/// rstd's unit tests for `target` and `optimize`, from a module of
/// `addRstdModule` (the registered one has no mode; its consumers set theirs).
pub fn addRstdTest(b: *Build, mod: *Module, target: ResolvedTarget, optimize: Optimize, filters: []const []const u8) !*Compile {
    const test_mod = try rstdbuild.cloneModuleWithTarget(b, mod, target);
    test_mod.optimize = optimize;
    if (tsan and target.result.os.tag == .linux and target.result.abi.isGnu()) test_mod.sanitize_thread = true;
    return b.addTest(.{ .root_module = test_mod, .filters = filters });
}

/// zfls's release targets: the `Cross-compile and package` step of zfls's
/// `.github/workflows/ci.yml` (macOS 11 and later, as there).
const zfls_targets = blk: {
    const macos11: std.Target.Query.OsVersion = .{ .semver = .{ .major = 11, .minor = 0, .patch = 0 } };
    break :blk [_]std.Target.Query{
        .{ .cpu_arch = .x86_64, .os_tag = .linux },
        .{ .cpu_arch = .aarch64, .os_tag = .linux },
        .{ .cpu_arch = .x86_64, .os_tag = .macos, .os_version_min = macos11 },
        .{ .cpu_arch = .aarch64, .os_tag = .macos, .os_version_min = macos11 },
        .{ .cpu_arch = .x86_64, .os_tag = .windows },
        .{ .cpu_arch = .aarch64, .os_tag = .windows },
    };
};

pub fn build(b: *Build) !void {
    // standard build setup
    const optimize = rstdbuild.getOptimize(b);
    const target = rstdbuild.getTarget(b);
    const docs = rstdbuild.addDocs(b);

    // Debug builds link only the target's BLAS kernels, which is far faster
    // to compile; `-Dblas-cores=all` builds every core to test dispatch.
    const blas_cores = b.option(rstdbuild.BlasCores, "blas-cores", "BLAS kernel objects to build: the target's core only (debug default), all, or none") orelse
        blaskern.defaultCores(target, optimize);
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
        // one shader per shaderlib operation, and rgl's own, to disassemble
        // (src/rgl/shaders/lib/readme.md)
        const probes = b.step("shader-probes", "Build each shaderlib probe and rgl's shaders to SPIR-V in zig-out/spv");
        const only = b.option([]const u8, "shader-probe", "shader-probes builds only this probe");
        inline for (@typeInfo(shader_probe.Op).@"enum".field_names) |op| {
            if (only == null or std.mem.eql(u8, only.?, op)) {
                const opts = b.addOptions();
                opts.addOption([]const u8, "op", op);
                const s = try rstdbuild.addShaderWith(b, optimize, b.path("src/rgl/shaders/probe.zig"), op, .{ ._rstd_mod = rstd_mod }, &.{.{ .name = "probe_op", .module = opts.createModule() }});
                probes.dependOn(&b.addInstallFile(s.spv_frag_path, "spv/" ++ op ++ ".frag.spv").step);
            }
        }
        for ([_][]const u8{ "rgl", "glyph", "image" }) |n| {
            const s = try rstdbuild.addShader(b, optimize, b.path(b.fmt("src/rgl/shaders/{s}.zig", .{n})), n, .{ ._rstd_mod = rstd_mod });
            probes.dependOn(&b.addInstallFile(s.spv_frag_path, b.fmt("spv/{s}.frag.spv", .{n})).step);
            probes.dependOn(&b.addInstallFile(s.spv_vert_path, b.fmt("spv/{s}.vert.spv", .{n})).step);
        }
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
    // glibc only: static on musl it faults at start; macOS needs SDK headers
    // the toolchain lacks (issues.md)
    tsan = b.option(bool, "rstd-tsan", "Build rstd's linux-gnu tests with ThreadSanitizer; the x86_64 and aarch64 Linux rigs then take glibc from Nix on -Dtest-host") orelse false;
    const test_filter = b.option([]const []const u8, "rstd-test-filter", "Only run rstd tests whose name contains this (repeatable)") orelse &.{};
    const rstd_test_exe = try addRstdTest(b, rstd_mod, target, optimize, test_filter);
    const run_rstd_test = b.addRunArtifact(rstd_test_exe);
    const test_step = b.step("rstd-test", "Run rstd's unit tests");
    test_step.dependOn(&run_rstd_test.step);
    // compiles exactly what `rstd-test` runs, without running it, so a build
    // queue can hold its slot for compiling only and run the tests after
    b.step("rstd-test-build", "Build rstd's unit tests as `rstd-test` would, without running them").dependOn(&rstd_test_exe.step);
    // the same executable and its debug info in <prefix>/rstd-test, so one
    // cross-built elsewhere can be fetched and run here
    const install_test = b.addInstallArtifact(rstd_test_exe, .{ .dest_dir = .{ .override = .{ .custom = "rstd-test" } } });
    b.step("rstd-test-install", "Install the executable `rstd-test` runs, with its debug info, to <prefix>/rstd-test").dependOn(&install_test.step);
    // as 32-bit x87 code, in Debug and ReleaseSafe: Debug's stores hide x87
    // width bugs (T87); Windows and Linux on x86 hosts
    const x87_step = b.step("rstd-test-x87", "Run rstd's unit tests as pentium3 code on this x86 host, in Debug and ReleaseSafe");
    const x87_build = b.step("rstd-test-x87-build", "Build what rstd-test-x87 runs, without running it");
    const host = b.graph.host.result;
    if (host.cpu.arch.isX86() and (host.os.tag == .windows or host.os.tag == .linux)) {
        const x87 = b.resolveTargetQuery(.{ .cpu_arch = .x86, .os_tag = host.os.tag, .abi = if (host.os.tag == .windows) .gnu else .musl, .cpu_model = .{ .explicit = &std.Target.x86.cpu.pentium3 } });
        for ([_]Optimize{ .debug, .safe }) |mode| {
            const x87_mod = try addRstdModule(b, null, x87, mode, blas_cores, false, false);
            const x87_exe = try addRstdTest(b, x87_mod, x87, mode, test_filter);
            x87_step.dependOn(&b.addRunArtifact(x87_exe).step);
            x87_build.dependOn(&x87_exe.step);
        }
    } else x87_step.dependOn(&b.addFail("rstd-test-x87 runs on x86 Windows and Linux hosts").step);
    // the same tests on real machines over SSH, for targets no emulator covers
    try sshtest.addSshTests(b, optimize, test_filter);

    // Mach-O, which no machine here runs: building is what catches assembly
    // that other object formats accept (`.align`), so build every core's
    // kernels for both architectures. `qemu-test-x86_64-macos` runs them.
    const macho_step = b.step("rstd-check-macho", "Build rstd's unit tests for macOS (x86_64, aarch64) without running them");
    for ([_]std.Target.Cpu.Arch{ .x86_64, .aarch64 }) |arch| {
        const macho = b.resolveTargetQuery(.{ .cpu_arch = arch, .os_tag = .macos });
        const macho_mod = try addRstdModule(b, null, macho, optimize, .all, false, false);
        macho_step.dependOn(&(try addRstdTest(b, macho_mod, macho, optimize, &.{})).step);
    }

    // UEFI: every OS-free test (the qemu rigs' root: simd, BLAS, every
    // interface at every element type, and LAPACK), also those the UEFI rigs
    // (`qemu-test-<arch>-uefi`) filter out, built with their runner, not run
    const uefi_step = b.step("rstd-check-uefi", "Build rstd's OS-free tests (simd, BLAS, LAPACK) for UEFI (x86_64, aarch64) without running them");
    for ([_]std.Target.Cpu.Arch{ .x86_64, .aarch64 }) |arch| {
        const uefi = b.resolveTargetQuery(.{ .cpu_arch = arch, .os_tag = .uefi });
        const uefi_mod = b.createModule(.{ .root_source_file = b.path("src/_test_qemu.zig"), .target = uefi, .optimize = optimize });
        try blaskern.addBlasKernels(b, uefi_mod, uefi, .all, null);
        try qemutest.addUefiImports(b, uefi_mod);
        uefi_step.dependOn(&b.addTest(.{
            .name = b.fmt("rstd-check-uefi-{t}", .{arch}),
            .root_module = uefi_mod,
            .test_runner = .{ .path = b.path("test/qemu/runner.zig"), .mode = .simple },
        }).step);
    }

    // zfls vendors rstd and builds it with Zig 0.17.0 for these targets
    // (design.md, Requirements); every core's kernels, built, not run
    const zig017_step = b.step("rstd-check-zig017", "Build rstd's unit tests for zfls's targets without running them; run with Zig 0.17.0");
    if (sysinfo.zig_version.major != 0 or sysinfo.zig_version.minor != 17)
        zig017_step.dependOn(&b.addFail("rstd-check-zig017 checks Zig 0.17.0; run it with that zig").step);
    for (zfls_targets) |query| {
        const t = b.resolveTargetQuery(query);
        const t_mod = try addRstdModule(b, null, t, optimize, .all, false, false);
        const t_test = try addRstdTest(b, t_mod, t, optimize, test_filter);
        // emitted, so codegen and the link run too, as in zfls's builds
        _ = t_test.getEmittedBin();
        zig017_step.dependOn(&t_test.step);
    }

    // the same tests that need no OS, on emulated machines
    try qemutest.addQemuTests(b, optimize);
    try wasmtest.addWasmTests(b, optimize);
}
