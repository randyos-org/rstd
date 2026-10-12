//! `wasm-test` step: runs the same OS-free tests as `qemu-test` as a WASI
//! module under Node (`tools/wasm_run.js`), for the WebAssembly kernels, which
//! no emulator here runs. The module is built with SIMD128 on, so the SIMD
//! core's kernels are the ones exercised.

const std = @import("std");
const rstdbuild = @import("../src/buildutils/__root__.zig");
const Build = rstdbuild.Build;
const Optimize = rstdbuild.Optimize;
const blaskern = @import("blaskern.zig");

/// Same selection as the emulated machines.
const default_filters: []const []const u8 = &.{ "linalg.blas.blas1.axpy", "simd." };

pub fn addWasmTests(b: *Build, optimize: Optimize) void {
    const step = b.step("wasm-test", "Run rstd's hardware-independent tests as a WASI module under Node");
    const target = b.resolveTargetQuery(.{
        .cpu_arch = .wasm32,
        .os_tag = .wasi,
        .cpu_features_add = std.Target.wasm.featureSet(&.{.simd128}),
    });
    const test_mod = b.createModule(.{
        .root_source_file = b.path("src/_test_qemu.zig"),
        .target = target,
        .optimize = optimize,
    });
    blaskern.addBlasKernels(b, test_mod, target, .all);
    const test_exe = b.addTest(.{
        .name = "rstd-wasm",
        .root_module = test_mod,
        .filters = default_filters,
    });

    const run = b.addSystemCommand(&.{"node"});
    run.addFileArg(b.path("tools/wasm_run.js"));
    run.addArtifactArg(test_exe);
    run.has_side_effects = true;
    step.dependOn(&run.step);
}
