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
pub const LazyPath = Build.LazyPath;
pub const Dependency = Build.Dependency;

pub fn getOptimize(b: *Build) OptimizeMode {
    return b.standardOptimizeOption(.{});
}

pub fn getTarget(b: *Build) ResolvedTarget {
    return b.standardTargetOptions(.{});
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

pub fn addBuildOptionsModule(b: *Build, pkg_name: []const u8, version: []const u8) *BuildOptions {
    const options = b.addOptions();
    options.addOption([]const u8, "package_name", pkg_name);
    options.addOption([]const u8, "version", version);
    return options;
}

pub const LogLevel = std.log.Level;
const LogLevelInt = @typeInfo(LogLevel).@"enum".tag_type;
const ScopeNames: type = []const []const u8;
pub const ScopeLevel = struct {
    scope: []const u8,
    level: LogLevel,
};

pub fn addLogScopeOptions(
    b: *Build,
    build_options: *BuildOptions,
    default_scope_levels: []const ScopeLevel,
) !void {
    const no_log_scope_args = b.option(
        ScopeNames,
        "no-log-scope",
        "Restrict a scope to error-level-only logging, on top of the default list (repeatable; pass with no value to clear the list, including defaults)",
    ) orelse &.{};
    const log_scope_args = b.option(
        ScopeNames,
        "log-scope",
        "Clear a scope's level override, restoring it to the global default log level (repeatable)",
    ) orelse &.{};
    const log_scope_level_args = b.option(
        ScopeNames,
        "log-scope-level",
        "Set a scope's log level explicitly, as \"scope=level\" (level one of err/warn/info/debug; repeatable)",
    ) orelse &.{};

    const reset_defaults = for (no_log_scope_args) |scope| {
        if (scope.len == 0) break true;
    } else false;

    // Values are @backingInt(LogLevel); logging.zig converts back with @enumFromInt.
    var levels: std.StringArrayHashMapUnmanaged(LogLevelInt) = .empty;
    if (!reset_defaults) {
        for (default_scope_levels) |dsl| try levels.put(b.allocator, dsl.scope, @backingInt(dsl.level));
    }
    for (no_log_scope_args) |scope| {
        if (scope.len == 0) continue;
        try levels.put(b.allocator, scope, @backingInt(LogLevel.err));
    }
    for (log_scope_args) |scope| {
        _ = levels.swapRemove(scope);
    }
    for (log_scope_level_args) |arg| {
        const eq = std.mem.indexOfScalar(u8, arg, '=') orelse
            std.debug.panic("-Dlog-scope-level value \"{s}\" is missing '=' (expected \"scope=level\")", .{arg});
        const scope = arg[0..eq];
        const level_name = arg[eq + 1 ..];
        const level = std.meta.stringToEnum(LogLevel, level_name) orelse
            std.debug.panic("-Dlog-scope-level value \"{s}\" has unknown level \"{s}\" (expected err, warn, info, or debug)", .{ arg, level_name });
        try levels.put(b.allocator, scope, @backingInt(level));
    }

    build_options.addOption([]const []const u8, "logger_scope_names", levels.keys());
    build_options.addOption([]const LogLevelInt, "logger_scope_name_levels", levels.values());
}

pub fn addBuildOptsModToModule(build_options: *BuildOptions, module: *Module) void {
    module.addOptions("build_options", build_options);
}

/// Shallow-clones `existing` (same root source file, imports, and every
/// other option already set on it) into a new `Module` bound to `target`.
///
/// A `Compile` step (`addObject`/`addExecutable`/`addTest`/...) requires its
/// root module to have a known target, but a module meant to be consumed via
/// `addImport` - like the one returned by `dep.module(...)` - is usually left
/// target-less on purpose so each consumer can pick its own. This lets a
/// standalone `Compile` step (docs, tests, ...) reuse that module's contents
/// without redeclaring `root_source_file` and re-adding every import/option.
pub fn cloneModuleWithTarget(b: *Build, existing: *Module, target: ?ResolvedTarget) !*Module {
    const clone = try b.allocator.create(Module);
    clone.* = existing.*;
    clone.resolved_target = target;
    // The struct copy above also copies the memoized import-graph cache
    // (keyed off `existing`'s own identity as "root"). Reset it so the
    // clone - a distinct root - recomputes its own instead of silently
    // reusing one rooted at `existing`.
    clone.cached_graph = .{ .modules = &.{}, .names = &.{} };
    return clone;
}

/// convert a given
pub fn installDirLazyPath(b: *Build, instdir: *InstallDir) LazyPath {
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

pub const ZigCLink = struct {
    lib: *CompileStep,
    incl: ?LazyPath = null,

    pub fn addLinkToModule(self: *const ZigCLink, mod: *Module) void {
        mod.linkLibrary(self.lib);
        if (self.incl) |incl| mod.addIncludePath(incl);
    }
};

pub const ZigCBind = struct {
    mod: *Module,
    clink: ZigCLink,
};

pub const ZigModule = struct {
    mod: *Module,
    name: []const u8,
    /// Set when `mod` is a thin Zig wrapper over a compiled C/C++ library
    /// (freetype, harfbuzz) - the underlying library still has to be
    /// linked into whatever consumes this module, since importing the
    /// wrapper module alone only pulls in type/function declarations.
    clink: ?ZigCLink = null,
    dep: ?*Dependency = null,

    pub fn addImportToModule(self: *const ZigModule, mod: *Module) void {
        mod.addImport(self.name, self.mod);
        if (self.clink) |clink| clink.addLinkToModule(mod);
    }
};

/// add a FreeDesktop `.desktop` file and icon to the install step for Linux targets
pub fn addLinuxDesktopFile(
    b: *Build,
    target: ResolvedTarget,
    icon_png: LazyPath,
    name: []const u8,
    display_name: []const u8,
    version: []const u8,
) !void {
    if (target.result.os.tag != .linux) return;

    const desktop_name = try std.fmt.allocPrint(b.allocator, "{s}.desktop", .{name});
    const desktop_contents = try std.fmt.allocPrint(b.allocator,
        \\[Desktop Entry]
        \\Type=Application
        \\Version=1.0
        \\Name={s}
        \\Comment={s} v{s}
        \\Exec={s}
        \\Icon={s}
        \\Terminal=false
        \\Categories=Utility;
        \\
    , .{ display_name, display_name, version, name, name });
    const build_dir = addBuildDir(b);
    const desktop_path = build_dir.add(desktop_name, desktop_contents);

    const icon_name = try std.fmt.allocPrint(b.allocator, "{s}.png", .{name});
    const install_step = b.getInstallStep();
    install_step.dependOn(&b.addInstallFileWithDir(
        desktop_path,
        .prefix,
        b.pathJoin(&.{ "share/applications", desktop_name }),
    ).step);
    install_step.dependOn(&b.addInstallFileWithDir(
        icon_png,
        .prefix,
        b.pathJoin(&.{ "share/icons/hicolor/256x256/apps", icon_name }),
    ).step);
}

pub const Icon = struct {
    const Self = @This();
    rc_build_path: LazyPath,
    png_build_path: LazyPath,
    png_name: []const u8,
    rgba_build_path: LazyPath,
    rgba_name: []const u8,

    pub fn addToModule(self: *const Self, mod: *Module) void {
        mod.addWin32ResourceFile(.{ .file = self.rc_build_path });
        mod.addAnonymousImport(self.png_name, .{ .root_source_file = self.png_build_path });
        mod.addAnonymousImport(self.rgba_name, .{ .root_source_file = self.rgba_build_path });
    }

    pub fn addGenIconStep(
        b: *Build,
        optimize: OptimizeMode,
        in_png_path: LazyPath,
        png_name: []const u8,
        ico_name: []const u8,
        rc_name: []const u8,
        rgba_name: []const u8,
    ) !Icon {
        const target = b.graph.host;
        const rstd = try addRstd(b, target, null);

        const mod = b.createModule(.{
            .root_source_file = rstd.dep.?.path("src/rgl/gen_icon.zig"),
            .target = target,
            .optimize = optimize,
        });
        rstd.addImportToModule(mod);

        const exe = b.addExecutable(.{
            .name = "genicon",
            .root_module = mod,
        });
        const run = b.addRunArtifact(exe);
        run.addFileArg(in_png_path);

        // record the paths of the generated output files so the build graph knows about them
        const png_build_path = run.addOutputFileArg(png_name);
        const ico_build_path = run.addOutputFileArg(ico_name);
        const rgba_build_path = run.addOutputFileArg(rgba_name);

        // The .rc ICON directive references "icon.ico" by a bare name, resolved
        // relative to the .rc file's own directory - so the generated icon.ico
        // needs to sit next to a (likewise generated) icon.rc, not just exist
        // somewhere else in the cache under its own path.
        const rc_build_dir = addBuildDir(b);
        const rc_build_path = rc_build_dir.add(rc_name, b.fmt("1 ICON \"{s}\"\n", .{ico_name}));
        _ = rc_build_dir.addCopyFile(ico_build_path, ico_name);

        // `zig build icon` - not needed for the app build itself (that already
        // depends on genicon's outputs above), but installs copies to zig-out/
        // so a human can open them after tweaking the drawing code.
        const icon_step = b.step("icon", "Regenerate icon.png/icon.ico into zig-out/ for inspection");
        icon_step.dependOn(&b.addInstallFileWithDir(png_build_path, .prefix, png_name).step);
        icon_step.dependOn(&b.addInstallFileWithDir(ico_build_path, .prefix, ico_name).step);
        icon_step.dependOn(&b.addInstallFileWithDir(rgba_build_path, .prefix, rgba_name).step);

        return .{
            .rc_build_path = rc_build_path,
            .png_build_path = png_build_path,
            .png_name = png_name,
            .rgba_build_path = rgba_build_path,
            .rgba_name = rgba_name,
        };
    }
};

pub fn addRstd(
    b: *Build,
    target: ?ResolvedTarget,
    docs: ?Docs,
) !ZigModule {
    const rstd_dep = b.dependency("rstd", .{
        .target = target,
    });
    // rstd_mod already has its own "build_options" import, wired up by
    // rstd's own build.zig when it configured itself - don't re-add one
    // here. `Module.addOptions`/`Options.createModule()` isn't memoized, so
    // doing that would stamp out a second Module rooted at the exact same
    // generated file the caller's own module wraps too, which the compiler
    // now rejects ("file exists in modules 'build_options' and
    // 'build_options0'").
    const rstd_mod = rstd_dep.module("rstd");

    if (docs) |d| {
        // We want the rstd docs to be built separately from the other packages,
        // so we will make a separate compile step for them and attach the docs.
        // The registered "rstd" module has no target (consumers set their own via
        // addImport), so clone it with one bound rather than redeclaring
        // root_source_file/imports/options.
        const rstd_docs_mod = try cloneModuleWithTarget(b, rstd_mod, target);
        // add the standalone rstd build step, just generates an unlinked object file
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

pub fn addWin32(b: *Build) ZigModule {
    return .{ .mod = b.dependency("win32", .{}).module("win32"), .name = "win32" };
}

pub fn addWin32ToModIfWinTgt(b: *Build, target: ResolvedTarget, mod: *Module) void {
    if (target.result.os.tag == .windows) {
        mod.link_libc = true;
        const w32 = addWin32(b);
        w32.addImportToModule(mod);
    }
}
