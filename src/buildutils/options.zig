const std = @import("std");

const rstdbuild = @import("__root__.zig");
const Build = rstdbuild.Build;
const Module = rstdbuild.Module;
const BuildOptions = rstdbuild.BuildOptions;
const DepModule = @import("deps.zig").DepModule;

/// Creates a package's `build_options`. Pass the package's rstd module
/// (from `addRstd`) to make rstd read these options too: rstd code such as
/// `rstd.Logger` imports `build_options` and expects the application's
/// settings (e.g. `addLogScopeOptions`), falling back to rstd's own when no
/// parent provides them.
pub fn addBuildOptionsModule(b: *Build, pkg_name: []const u8, version: []const u8, rstd: ?*const DepModule) !*BuildOptions {
    const options = b.addOptions();
    options.addOption([]const u8, "package_name", pkg_name);
    options.addOption([]const u8, "version", version);
    if (rstd) |r| r.mod.addImport("build_options", try optionsModule(options));
    return options;
}

/// One module per options step. `Options.createModule` (and so
/// `Module.addOptions`) wraps the generated file in a new module on every
/// call, and two such modules in one compilation are rejected ("file exists
/// in modules 'build_options' and 'build_options0'"), so everything that
/// imports a given options step must share this module.
var options_modules: std.array_hash_map.Auto(*BuildOptions, *Module) = .empty;

pub fn optionsModule(options: *BuildOptions) !*Module {
    const gop = try options_modules.getOrPut(options.step.owner.allocator, options);
    if (!gop.found_existing) gop.value_ptr.* = options.createModule();
    return gop.value_ptr.*;
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
    var levels: std.array_hash_map.String(LogLevelInt) = .empty;
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
        const eq = std.mem.findScalar(u8, arg, '=') orelse
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

pub fn addBuildOptsModToModule(build_options: *BuildOptions, module: *Module) !void {
    module.addImport("build_options", try optionsModule(build_options));
}
