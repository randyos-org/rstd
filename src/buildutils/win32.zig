const rstdbuild = @import("__root__.zig");
const Build = rstdbuild.Build;
const ResolvedTarget = rstdbuild.ResolvedTarget;
const Module = rstdbuild.Module;
const DepModule = @import("deps.zig").DepModule;

pub fn addWin32(b: *Build) DepModule {
    return .{ .mod = b.dependency("win32", .{}).module("win32"), .name = "win32" };
}

pub fn addWin32ToModIfWinTgt(b: *Build, target: ResolvedTarget, mod: *Module) void {
    if (target.result.os.tag == .windows) {
        mod.link_libc = true;
        const w32 = addWin32(b);
        w32.addImportToModule(mod);
    }
}
