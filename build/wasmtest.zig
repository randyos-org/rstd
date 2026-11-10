//! `wasm-test` step: the rigs' OS-free tests (`src/_test_qemu.zig`) as WASI modules under Node (`tools/wasm_run.js`), for the wasm
//! kernels no emulator here runs: one module with SIMD128 for the SIMD core, one without for the portable path.

const std = @import("std");
const rstdbuild = @import("../src/buildutils/__root__.zig");
const Build = rstdbuild.Build;
const Optimize = rstdbuild.Optimize;
const blaskern = @import("blaskern.zig");
const machines = @import("qemu/machines.zig");

/// The rigs' selection, plus the files WASI can also run.
const default_filters: []const []const u8 = machines.default_filters ++ &[_][]const u8{
    "linalg.blas.abi.",
    "linalg.blas.threading.",
};

const Variant = struct { name: []const u8, features: []const std.Target.wasm.Feature };

const variants = [_]Variant{
    .{ .name = "simd128", .features = &.{.simd128} },
    .{ .name = "portable", .features = &.{} },
};

pub fn addWasmTests(b: *Build, optimize: Optimize) !void {
    const step = b.step("wasm-test", "Run rstd's OS-free tests as WASI modules under Node, with and without SIMD128");
    const check = b.step("wasm-check", "Build wasm-test's modules without running them");
    const filters = b.option([]const []const u8, "wasm-test-filter", "Only run wasm tests whose name contains this (repeatable)") orelse default_filters;
    // one module compiles at a time: two at once can pass a CI runner's memory (issue #268)
    var prev: ?*Build.Step = null;
    for (variants) |v| {
        const one = b.step(b.fmt("wasm-test-{s}", .{v.name}), b.fmt("Run wasm-test's {s} module only", .{v.name}));
        const one_check = b.step(b.fmt("wasm-check-{s}", .{v.name}), b.fmt("Build wasm-test's {s} module without running it", .{v.name}));
        const target = b.resolveTargetQuery(.{
            .cpu_arch = .wasm32,
            .os_tag = .wasi,
            .cpu_features_add = std.Target.wasm.featureSet(v.features),
        });
        const test_mod = b.createModule(.{
            .root_source_file = b.path("src/_test_qemu.zig"),
            .target = target,
            .optimize = optimize,
        });
        try blaskern.addBlasKernels(b, test_mod, target, .all, null);
        const test_exe = b.addTest(.{
            .name = b.fmt("rstd-wasm-{s}", .{v.name}),
            .root_module = test_mod,
            .filters = filters,
        });
        if (prev) |p| test_exe.step.dependOn(p);
        prev = &test_exe.step;

        const run = b.addSystemCommand(&.{"node"});
        run.setName(b.fmt("run rstd-wasm-{s}", .{v.name}));
        run.addFileArg(b.path("tools/wasm_run.js"));
        run.addArtifactArg2(test_exe, .{});
        run.has_side_effects = true;
        one.dependOn(&run.step);
        step.dependOn(&run.step);
        one_check.dependOn(&test_exe.step);
        check.dependOn(&test_exe.step);
    }
}
