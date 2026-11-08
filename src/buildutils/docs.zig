const rstdbuild = @import("__root__.zig");
const Build = rstdbuild.Build;
const Step = rstdbuild.Step;
const CompileStep = rstdbuild.CompileStep;
const Module = rstdbuild.Module;
const ResolvedTarget = rstdbuild.ResolvedTarget;
const console = rstdbuild.console;
const cloneModuleWithTarget = rstdbuild.cloneModuleWithTarget;

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
    /// things always expected to build; other archs should call
    /// addCompileStepDocsToStep against their own step so a broken stub
    /// can't break the build.
    pub fn addCompileStepDocs(docs: Docs, b: *Build, compile: *CompileStep, name: []const u8) void {
        addCompileStepDocsToStep(b, docs.step, compile, name);
    }

    /// `mod`'s docs, through an unlinked object named `name` built for `target`.
    /// A target-less module (one registered for consumers) is cloned with one.
    pub fn addModuleDocs(docs: Docs, b: *Build, mod: *Module, target: ?ResolvedTarget, name: []const u8) !void {
        const docs_mod = try cloneModuleWithTarget(b, mod, target);
        // Docs only analyze the source; the BLAS kernel objects would make
        // this object link others, which the self-hosted COFF linker rejects.
        docs_mod.link_objects = .empty;
        docs.addCompileStepDocs(b, b.addObject(.{ .name = name, .root_module = docs_mod }), name);
    }
};

pub fn addDocs(b: *Build) Docs {
    const docs_step = b.step("docs", "Generate and install API documentation");
    // also runs as part of default zig build/install, not just explicit request
    b.getInstallStep().dependOn(docs_step);

    const docs_port = b.option(u16, "docs-port", "Port for `zig build run-docs` to serve documentation on") orelse 3000;
    const run_docs_step = b.step("run-docs", "Serve the generated documentation locally");
    const server_cmd = console(b, &.{
        "python", "-m", "http.server", b.fmt("{d}", .{docs_port}), "-b", "127.0.0.1", "-d", "zig-out/docs/",
    });
    server_cmd.step.dependOn(docs_step);
    run_docs_step.dependOn(&server_cmd.step);

    return .{ .step = docs_step };
}
