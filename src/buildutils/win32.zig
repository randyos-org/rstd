const rstdbuild = @import("__root__.zig");
const Build = rstdbuild.Build;
const Module = rstdbuild.Module;
const ResolvedTarget = rstdbuild.ResolvedTarget;
const DepModule = rstdbuild.DepModule;

/// zigwin32, resolved through rstd's own build.zig.zon so packages that use
/// rstd get rstd's pinned zigwin32 without declaring it themselves. A package
/// that does declare `win32` (rstd itself, or one that wants its own pin)
/// gets that one instead.
pub fn addWin32(b: *Build) DepModule {
    return .{ .mod = win32Owner(b).dependency("win32", .{}).module("win32"), .name = "win32" };
}

pub fn addWin32ToModIfWinTgt(b: *Build, target: ResolvedTarget, mod: *Module) void {
    if (target.result.os.tag == .windows) {
        mod.link_libc = true;
        const w32 = addWin32(b);
        w32.addImportToModule(mod);
    }
}

/// The builder whose build.zig.zon provides `win32`: the caller's own when it
/// declares one, otherwise rstd's. Dependencies are cached by package hash and
/// options, so any rstd instance yields the same zigwin32 module as the one
/// rstd's own module imports.
fn win32Owner(b: *Build) *Build {
    return if (rstdbuild.declaresDependency(b, "win32")) b else b.dependency("rstd", .{}).builder;
}
