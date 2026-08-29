// const std = @import("std");
const rstdbuild = @import("__root__.zig");
const Build = rstdbuild.Build;
const Module = rstdbuild.Module;
const CompileStep = rstdbuild.CompileStep;
const LazyPath = rstdbuild.LazyPath;
const Dependency = rstdbuild.Dependency;

pub const DepCLink = struct {
    lib: *CompileStep,
    incl: ?LazyPath = null,

    pub fn addLinkToModule(self: *const DepCLink, mod: *Module) void {
        mod.linkLibrary(self.lib);
        if (self.incl) |incl| mod.addIncludePath(incl);
    }
};

pub const DepCBind = struct {
    mod: *Module,
    clink: DepCLink,
};

pub const DepModule = struct {
    mod: *Module,
    name: []const u8,
    /// Set when `mod` is a thin Zig wrapper over a compiled C/C++ library
    /// (freetype, harfbuzz); the underlying library still has to be
    /// linked into whatever consumes this module, since importing the
    /// wrapper module alone only pulls in type/function declarations.
    clink: ?DepCLink = null,
    dep: ?*Dependency = null,

    pub fn addImportToModule(self: *const DepModule, mod: *Module) void {
        mod.addImport(self.name, self.mod);
        if (self.clink) |clink| clink.addLinkToModule(mod);
    }
};
