const rstdbuild = @import("__root__.zig");
const docs_mod = @import("docs.zig");
const module = @import("module.zig");

const Build = rstdbuild.Build;
const ResolvedTarget = rstdbuild.ResolvedTarget;
const Docs = docs_mod.Docs;
const DepModule = @import("deps.zig").DepModule;
const Module = rstdbuild.Module;

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
    // the native target and pulling in libc) -- fine for docs/tests, but
    // fatal for a spirv32-vulkan `target` (no libc), which is the case this
    // branch (`_rstd_mod` reuse, from `wireRglShaders`) exists for.
    const rstd_mod = if (options._rstd_mod) |rmod| module.cloneRootForTarget(b, rmod, target) else rstd_dep.?.module("rstd");

    if (docs) |d| {
        // The registered "rstd" module has no target (consumers set their
        // own via addImport), so clone it with one bound rather than
        // redeclaring root_source_file/imports/options.
        const rstd_docs_mod = try module.cloneModuleWithTarget(b, rstd_mod, target);
        // Docs only analyze the source; the BLAS kernel objects would make
        // this object link others, which the self-hosted COFF linker rejects.
        rstd_docs_mod.link_objects = .empty;
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
