const std = @import("std");
pub const Build = std.Build;
pub const Step = Build.Step;
pub const BuildDir = Step.WriteFile;
pub const InstallDir = Step.InstallDir;
pub const OptimizeMode = std.builtin.OptimizeMode;
pub const TargetQuery = std.Target.Query;
pub const RunStep = Step.Run;
pub const CompileStep = Step.Compile;
pub const Module = Build.Module;
pub const BuildOptions = Step.Options;
pub const ResolvedTarget = Build.ResolvedTarget;

pub fn getOptimize(b: *Build) OptimizeMode {
    return b.standardOptimizeOption(.{});
}

pub fn addBuildDir(b: *Build) *BuildDir {
    return b.addWriteFiles();
}

/// install compile's docs under zig-out/docs/<name>/, target_step depends on it
pub fn addCompileStepDocsToStep(b: *Build, target_step: *Step, compile: *CompileStep, name: []const u8) void {
    const install = b.addInstallDirectory(.{
        .source_dir = compile.getEmittedDocs(),
        .install_dir = .prefix,
        .install_subdir = b.fmt("docs/{s}", .{name}),
    });
    target_step.dependOn(&install.step);
}

/// common docs step handler
pub const Docs = struct {
    step: *Step,

    /// install compile's docs, shared "docs" step depends on it. Only for
    /// things always expected to build -- other archs should call
    /// addCompileStepDocsToStep against their own step so a broken stub
    /// can't break the build.
    pub fn addCompileStepDocs(docs: Docs, b: *Build, compile: *CompileStep, name: []const u8) void {
        addCompileStepDocsToStep(b, docs.step, compile, name);
    }
};

pub fn addDocs(b: *Build) Docs {
    const docs_step = b.step("docs", "Generate and install API documentation");
    // also runs as part of default zig build/install, not just explicit request
    b.getInstallStep().dependOn(docs_step);

    const docs_port = b.option(u16, "docs-port", "Port for `zig build run-docs` to serve documentation on") orelse 3000;
    const run_docs_step = b.step("run-docs", "Serve the generated documentation locally");
    const server_cmd = b.addSystemCommand(&.{
        "python", "-m", "http.server", b.fmt("{d}", .{docs_port}), "-b", "127.0.0.1", "-d", "zig-out/docs/",
    });
    server_cmd.step.dependOn(docs_step);
    run_docs_step.dependOn(&server_cmd.step);

    return .{ .step = docs_step };
}

pub fn addBuildOptionsModule(b: *Build) *BuildOptions {
    return b.addOptions();
}

pub fn addBuildOptsModToModule(build_options: *BuildOptions, module: *Module) void {
    module.addOptions("build_options", build_options);
}

/// convert a given
pub fn installDirLazyPath(b: *Build, instdir: *InstallDir) Build.LazyPath {
    const dest_rel_path = instdir.options.install_subdir;
    return switch (instdir.options.install_dir) {
        .prefix => .{ .relative = .{ .base = .install_prefix, .sub_path = dest_rel_path } },
        .lib => .{ .relative = .{ .base = .install_lib, .sub_path = dest_rel_path } },
        .bin => .{ .relative = .{ .base = .install_bin, .sub_path = dest_rel_path } },
        .header => .{ .relative = .{ .base = .install_include, .sub_path = dest_rel_path } },
        .custom => |p| .{ .relative = .{
            .base = .install_prefix,
            .sub_path = if (dest_rel_path.len == 0) p else b.pathJoin(&.{ p, dest_rel_path }),
        } },
    };
}
