const rstdbuild = @import("__root__.zig");
const Build = rstdbuild.Build;
const Module = rstdbuild.Module;
const ResolvedTarget = rstdbuild.ResolvedTarget;
const DepModule = rstdbuild.DepModule;
const Docs = rstdbuild.Docs;

/// Which cores get a BLAS kernel object: `target`, only the one the build
/// target maps to; `all`, every core it can run (for runtime dispatch);
/// `none`, no kernel objects, for consumers that leave BLAS unused (the
/// portable kernels stay installed).
pub const BlasCores = enum { target, all, none };

pub const RstdOptions = struct {
    _rstd_mod: ?*Module = null,
    /// The consumer's optimize mode; picks `blas_cores`' default (`target` in
    /// debug, else `all`). Unset, rstd treats the build as debug.
    optimize: ?rstdbuild.Optimize = null,
    blas_cores: ?BlasCores = null,
    /// RenderDoc capture support in `rgl.renderdoc`; fetches renderdoc and
    /// translate-c only when set.
    renderdoc: bool = false,
};

pub fn addRstd(
    b: *Build,
    target: ?ResolvedTarget,
    docs: ?Docs,
    vulkan: bool,
    options: RstdOptions,
) !DepModule {
    const rstd_dep = if (options._rstd_mod) |_| null else b.dependency("rstd", .{
        .target = target,
        .vulkan = vulkan,
        .optimize = options.optimize,
        .@"blas-cores" = options.blas_cores,
        .renderdoc = options.renderdoc,
    });
    // rstd_mod already has its own "build_options" import, wired up by
    // rstd's own build.zig when it configured itself; don't re-add one
    // here. `Module.addOptions`/`Options.createModule()` isn't memoized, so
    // doing that would stamp out a second Module rooted at the exact same
    // generated file the caller's own module wraps too, which the compiler
    // now rejects ("file exists in modules 'build_options' and
    // 'build_options0'"). Using `cloneRootForTarget` rather than
    // `cloneModuleWithTarget` here also matters: `rmod` is rstd's own native
    // `rstd_mod`, which by this point carries `win32`/`renderdoc` (added for
    // the native target and pulling in libc); fine for docs/tests, but
    // fatal for a spirv32-vulkan `target` (no libc), which is the case this
    // branch (`_rstd_mod` reuse, from `wireRglShaders`) exists for.
    const rstd_mod = if (options._rstd_mod) |rmod| rstdbuild.cloneRootForTarget(b, rmod, target) else rstd_dep.?.module("rstd");

    if (docs) |d| try d.addModuleDocs(b, rstd_mod, target, "rstd");

    return .{
        .mod = rstd_mod,
        .name = "rstd",
        .dep = rstd_dep,
    };
}
