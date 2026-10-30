//! Per-core BLAS kernel objects.
//!
//! OpenBLAS (`DYNAMIC_ARCH`) compiles its kernel tree once per core type with
//! that core's flags. This does the same: one LLVM-built object per group of
//! cores whose flags are equal, from `src/linalg/blas/_kernel.zig`, linked into
//! the `rstd` module so the main library can stay on the target's baseline CPU
//! and any backend.
//!
//! What each core is compiled with is not written here: it is
//! `src/linalg/blas/cpu/_corespec.zig`, shared with the runtime. This takes
//! the target architecture's cores with kernels (`corespec.kernelCores`),
//! keeps those the target can run, groups them by equal spec, and builds one
//! object per group.
//!
//! `-Dblas-cores=target` (the debug default) builds only the group the target
//! maps to: the host's detected core for a native target, else the group with
//! the most features the target has. A detected core not linked at runtime
//! uses that one instead (`asm_registry.install`).
//!
//! Supporting another architecture means filling in its `kernels/<arch>.zig`
//! (a stub exists for each) and listing its cores in `corespec.kernelCores`;
//! its specs are already in `_corespec.zig`.

const std = @import("std");
const sysinfo = @import("builtin");
const rstdbuild = @import("../src/buildutils/__root__.zig");
const Build = rstdbuild.Build;
const Module = rstdbuild.Module;
const Compile = rstdbuild.CompileStep;
const LazyPath = rstdbuild.LazyPath;
const ResolvedTarget = rstdbuild.ResolvedTarget;

// Only `cpu/`'s light files (`_corespec.zig` and the `_detect/<arch>.zig` detectors
// it names): the registries or `cpu/<arch>.zig` would make this script's
// compile depend on every kernel source.
const corespec = @import("../src/linalg/blas/cpu/_corespec.zig");
const abi = @import("../src/linalg/blas/abi.zig");

/// Cores with equal specs, served by one object.
const Group = struct {
    spec: corespec.Spec,
    query: std.Target.Query,
    cores: std.ArrayList([]const u8) = .empty,
};

var force_llvm: ?bool = null;
var relapack = false;
var option_read = false;

/// The kernel objects built when `-Dblas-cores` is not given: the target's
/// core in debug, else all; none where no kernel takes the target's calling
/// convention (`abi.zig`), so a consumer's build links the portable kernels
/// instead of stopping. Asking for kernels there still stops the build.
pub fn defaultCores(target: ResolvedTarget, optimize: rstdbuild.Optimize) rstdbuild.BlasCores {
    if (abi.variantOf(&target.result) == .unsupported) return .none;
    return if (optimize == .debug) .target else .all;
}

/// Builds the per-core kernel objects for `target` and links them into `mod`:
/// for every core the target can run (`.all`), or only the core the target
/// maps to (`.target`), or none (`.none`). It also gives `mod` the
/// `blaskern_linked` module, which tells `asm_registry.install` what got linked,
/// and `lapack_options` (`-Drelapack`), since every module with the BLAS has
/// LAPACK too.
pub fn addBlasKernels(b: *Build, mod: *Module, target: ResolvedTarget, cores: rstdbuild.BlasCores) !void {
    if (!option_read) {
        option_read = true;
        force_llvm = b.option(bool, "blas-llvm", "Force the BLAS kernel objects on (true) or off (false) LLVM; unset uses it only where the kernels need it");
        relapack = b.option(bool, "relapack", "LAPACK's potrf, getrf, trtri and lauum from ReLAPACK's recursive algorithms instead of OpenBLAS's, as OpenBLAS's BUILD_RELAPACK with RELAPACK_REPLACE") orelse false;
    }
    const lapack_opts = b.addOptions();
    lapack_opts.addOption(bool, "relapack", relapack);
    mod.addImport("lapack_options", lapack_opts.createModule());
    // the core type depends on the architecture, so it becomes compile-time here
    switch (target.result.cpu.arch) {
        inline else => |a| try addForArch(b, mod, target, cores, a),
    }
}

fn addForArch(
    b: *Build,
    mod: *Module,
    target: ResolvedTarget,
    cores: rstdbuild.BlasCores,
    comptime arch_tag: std.Target.Cpu.Arch,
) !void {
    const arch = target.result.cpu.arch;
    var groups: std.ArrayList(Group) = .empty;

    // the converted assembly has a path per calling convention (`abi.zig`)
    switch (abi.variantOf(&target.result)) {
        .variant => if (cores == .none) return addLinked(b, mod, &.{}, null),
        .none => return addLinked(b, mod, &.{}, null),
        .unsupported => |why| if (cores == .none) return addLinked(b, mod, &.{}, null) else std.process.fatal(
            "no BLAS kernels for {s}-{s}-{s}: {s}; pass -Dblas-cores=none for the portable ones",
            .{ @tagName(arch), @tagName(target.result.os.tag), @tagName(target.result.abi), why },
        ),
    }

    for (corespec.kernelCores(arch_tag)) |core| {
        // its kernels must run on the target CPU
        if (!corespec.runsOn(target.result.cpu, core)) continue;
        const spec = corespec.of(core);
        const g = for (groups.items) |*g| {
            if (g.spec.eql(spec)) break g;
        } else blk: {
            try groups.append(b.allocator, .{ .spec = spec, .query = groupQuery(target, spec) });
            break :blk &groups.items[groups.items.len - 1];
        };
        try g.cores.append(b.allocator, @tagName(core));
        try checkIsa(b, target, g.query, core);
    }

    var fallback: ?[]const u8 = null;
    if (cores == .target) {
        const pick = pickTarget(b, target, arch_tag, groups.items);
        if (pick) |p| {
            const g = groups.items[p.group];
            groups.clearRetainingCapacity();
            groups.appendAssumeCapacity(g);
            fallback = p.core;
        } else groups.clearRetainingCapacity();
    }

    var linked: std.ArrayList([]const u8) = .empty;
    for (groups.items) |g| try linked.appendSlice(b.allocator, g.cores.items);
    addLinked(b, mod, linked.items, fallback);

    for (groups.items) |g| {
        const opts = b.addOptions();
        opts.addOption([]const []const u8, "cores", g.cores.items);
        opts.addOption(bool, "fuse_fma_set", g.spec.fuse_fma != null);
        opts.addOption(bool, "fuse_fma", g.spec.fuse_fma orelse false);

        // LLVM only where the kernels need it; `-Dblas-llvm` overrides either way.
        const llvm = force_llvm orelse needsLlvm(arch, g.cores.items);
        const obj = b.addObject(.{
            .name = b.fmt("rstd-blaskernel-{s}", .{g.cores.items[0]}),
            .root_module = b.createModule(.{
                .root_source_file = b.path("src/linalg/blas/_kernel.zig"),
                .target = b.resolveTargetQuery(g.query),
                // Kernels are always built for speed, whatever the consumer uses.
                .optimize = .fast,
                .pic = objectPic(mod, target),
                .imports = &.{.{ .name = "blaskernel_options", .module = opts.createModule() }},
            }),
            .use_llvm = llvm,
            // the self-hosted linker writes the object otherwise, and pads and
            // misplaces SPARC symbols (calls land in zeroed `.text`); LLD cannot
            // write Mach-O, so macOS objects keep the self-hosted linker
            .use_lld = llvm and target.result.ofmt != .macho,
        });
        // `zig cc` ignores `-ofmt`, so an explicit object format goes unchecked
        if (llvm and g.query.ofmt == null) mod.addObjectFile(try checkedObject(b, obj, g.query)) else mod.addObject(obj);
    }
}

/// `obj`'s object, linkable only once `asmCheck` passes.
fn checkedObject(b: *Build, obj: *Compile, query: std.Target.Query) !LazyPath {
    const wf = b.addWriteFiles();
    wf.step.dependOn(try asmCheck(b, obj, query));
    return wf.addCopyFile(obj.getEmittedBin(), obj.out_filename);
}

/// A step that fails where `comp`'s inline asm has a line LLVM rejects: it
/// compiles `comp`'s LLVM bitcode with clang. Zig 0.18 drops such lines, yet
/// exits 0 and caches the output (issues.md, "Inline assembly errors exit 0").
/// `comp` must use LLVM, and `query` name no `ofmt`, which `zig cc` ignores.
pub fn asmCheck(b: *Build, comp: *Compile, query: std.Target.Query) !*Build.Step {
    const triple = try query.zigTriple(b.allocator);
    const check = b.addSystemCommand(&.{ b.graph.zig_exe, "cc", "-c", "-Wno-override-module", "-target", triple });
    const target = b.resolveTargetQuery(query).result;
    // Functions carry their own CPU; module-level asm needs `-mcpu`. Not on
    // x86, whose assembler ignores CPU features.
    if (!target.cpu.arch.isX86())
        check.addArg(b.fmt("-mcpu={s}", .{try query.serializeCpuAlloc(b.allocator)}));
    // as Zig builds it; `zig cc` makes freestanding mips PIC, which its asm
    // fails, and non-abicalls MIPS code (`objectPic`) cannot be PIC
    if (comp.root_module.pic == false or (comp.root_module.pic == null and target.os.tag == .freestanding))
        check.addArg("-fno-PIC");
    // wasm bitcode already has its `wasm-feature-*` flags, which clang adds again
    // (issues.md, "wasm bitcode gets its feature flags twice")
    if (target.cpu.arch.isWasm()) check.addArgs(&.{ "-Xclang", "-disable-llvm-verifier" });
    // a function per section, so calls are relocations the linker resolves; in
    // one section a large image's branches outrun their range (+-32 MiB on
    // PowerPC and 32-bit arm), which the assembler reports as an error
    check.addArg("-ffunction-sections");
    // an `SDKROOT` for another platform (Swift's for Windows) makes `zig cc`
    // warn for Darwin targets, and stderr fails the step; no SDK is read here
    check.removeEnvironmentVariable("SDKROOT");
    check.setName(b.fmt("check inline asm of {s}", .{comp.name}));
    check.addFileArg(comp.getEmittedLlvmBc());
    check.addArg("-o");
    _ = check.addOutputFileArg(b.fmt("asmcheck-{s}.o", .{comp.name}));
    return &check.step;
}

/// The objects' PIC where it must match the image's: MIPS marks PIC code
/// abicalls, and `ld.lld` warns when it meets the rest. An object is PIC by
/// default, an executable only where it links dynamically (glibc); a shared
/// library sets `pic` on its module. Elsewhere Zig's default, which links
/// into anything.
fn objectPic(mod: *Module, target: ResolvedTarget) ?bool {
    const t = &target.result;
    if (!t.cpu.arch.isMips()) return null;
    return mod.pic orelse switch (t.os.tag) {
        .freestanding => false,
        .linux => mod.link_libc == true and t.isGnuLibC(),
        else => null,
    };
}

/// Tells `asm_registry.install` which cores' objects are linked.
fn addLinked(b: *Build, mod: *Module, linked: []const []const u8, fallback: ?[]const u8) void {
    const opts = b.addOptions();
    opts.addOption([]const []const u8, "cores", linked);
    opts.addOption(?[]const u8, "fallback", fallback);
    mod.addImport("blaskern_linked", opts.createModule());
}

/// The target query a group's object is built with (`corespec.objectQuery`).
fn groupQuery(target: ResolvedTarget, spec: corespec.Spec) std.Target.Query {
    return corespec.objectQuery(target.query, target.result.os, target.result.cpu, spec);
}

/// Fails the build where `query`'s object needs a feature that a CPU running
/// `core` may lack and `target` does not promise (`corespec.excess`).
fn checkIsa(b: *Build, target: ResolvedTarget, query: std.Target.Query, core: anytype) !void {
    const arch = target.result.cpu.arch;
    const object = b.resolveTargetQuery(query).result.cpu.features;
    const ex = corespec.excess(core, arch, object, target.result.cpu.features);
    if (ex.isEmpty()) return;
    var names: std.ArrayList(u8) = .empty;
    for (arch.allFeaturesList(), 0..) |f, i| {
        if (ex.isEnabled(@intCast(i))) try names.print(b.allocator, " {s}", .{f.name});
    }
    std.process.fatal("the BLAS kernel object serving {s} on {s} uses{s}, which some {s} CPUs lack (`machine` in cpu/_corespec.zig)", .{
        @tagName(core), target.result.cpu.model.name, names.items, @tagName(core),
    });
}

/// The group, and the core in it, that `target` maps to. A native target asks
/// the host, as the runtime would (`corespec.Detector`); otherwise it is the
/// group with the most features that are all the target's.
fn pickTarget(
    b: *Build,
    target: ResolvedTarget,
    comptime arch_tag: std.Target.Cpu.Arch,
    groups: []const Group,
) ?struct { group: usize, core: []const u8 } {
    if (comptime arch_tag == sysinfo.target.cpu.arch) {
        if (target.query.isNativeCpu()) {
            const name = @tagName(corespec.objectCore(arch_tag, corespec.Detector(arch_tag).detect()));
            for (groups, 0..) |g, i| {
                for (g.cores.items) |c| {
                    if (std.mem.eql(u8, c, name)) return .{ .group = i, .core = c };
                }
            }
        }
    }
    var best: ?usize = null;
    var best_count: usize = 0;
    for (groups, 0..) |g, i| {
        const features = b.resolveTargetQuery(g.query).result.cpu.features;
        if (!target.result.cpu.features.isSuperSetOf(features)) continue;
        if (best == null or features.count() > best_count) {
            best = i;
            best_count = features.count();
        }
    }
    const i = best orelse return null;
    return .{ .group = i, .core = groups[i].cores.items[0] };
}

/// Whether the kernels of `cores` can only be built with LLVM, found by building
/// every architecture with `-Dblas-llvm=false`:
///
///   * wasm32/wasm64: yes. The self-hosted backend builds them but emits invalid
///     wasm for the SIMD kernel (node: "type error in branch, expected v128, got
///     i32" in `Real(f32).unit`).
///   * x86_64: yes for every core that has kernels. The self-hosted backend's own
///     assembler rejects the register-based memory operands (`(%rdx)`) of the
///     converted assembly and, past SSE, the VEX/EVEX mnemonics (`vzeroupper`).
///   * aarch64, riscv64 and the rest: yes; they have no self-hosted backend that
///     covers inline assembly and vectors, and riscv64 also fails on `zbb` in std.
fn needsLlvm(arch: std.Target.Cpu.Arch, cores: []const []const u8) bool {
    _ = arch;
    _ = cores;
    return true;
}
