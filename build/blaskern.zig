//! Per-core BLAS kernel objects.
//!
//! OpenBLAS (`DYNAMIC_ARCH`) compiles its kernel tree once per core type with
//! that core's flags. This does the same: one LLVM-built object per group of
//! cores whose flags are equal, from `src/linalg/blas/_kernel.zig`, linked into
//! the `rstd` module so the main library can stay on the target's baseline CPU
//! and any backend.
//!
//! What each core is compiled with is not written here: it is
//! `src/linalg/blas/cpu/_corespec.zig`, shared with the runtime. This walks
//! the target architecture's core types, keeps those a kernel registry
//! (`src/linalg/blas/kernels/<arch>.zig`) has kernels for, groups them by
//! equal spec, and builds one object per group.
//!
//! The architecture's core types and kernel registry come from
//! `cpu.GetCpuArch(arch)`. Supporting another architecture means filling in its
//! `kernels/<arch>.zig` (a stub exists for each); its specs are already in
//! `_corespec.zig`.

const std = @import("std");
const rstd = @import("../src/__root__.zig");
const rstdbuild = rstd.buildutils;
const Build = rstdbuild.Build;
const Module = rstdbuild.Module;
const ResolvedTarget = rstdbuild.ResolvedTarget;

const cpu = @import("../src/linalg/blas/cpu/__root__.zig");
const corespec = @import("../src/linalg/blas/cpu/_corespec.zig");

/// Cores with equal specs, served by one object.
const Group = struct {
    spec: corespec.Spec,
    cores: std.ArrayList([]const u8) = .empty,
};

var force_llvm: ?bool = null;
var option_read = false;

/// Builds the per-core kernel objects for `target` and links them into `mod`.
pub fn addBlasKernels(b: *Build, mod: *Module, target: ResolvedTarget) void {
    if (!option_read) {
        option_read = true;
        force_llvm = b.option(bool, "blas-llvm", "Force the BLAS kernel objects on (true) or off (false) LLVM; unset uses it only where the kernels need it");
    }
    // `GetCpuArch` picks a type, so the runtime architecture becomes a
    // compile-time one here; each arm is one instantiation per architecture family.
    switch (target.result.cpu.arch) {
        inline else => |a| {
            const arch = cpu.GetCpuArch(a);
            addForArch(b, mod, target, a, arch.CoreType, arch.kernel);
        },
    }
}

fn addForArch(b: *Build, mod: *Module, target: ResolvedTarget, comptime arch_tag: std.Target.Cpu.Arch, comptime Core: type, comptime registry: type) void {
    const arch = target.result.cpu.arch;
    if (!registry.supports(arch)) return;
    var groups: std.ArrayList(Group) = .empty;

    inline for (@typeInfo(Core).@"enum".field_names) |name| {
        const core = @field(Core, name);
        // its kernels must run on the target CPU
        const runs = if (@hasDecl(registry, "runsOn")) registry.runsOn(target.result.cpu, core) else true;
        if (comptime registry.hasKernels(arch_tag, core)) {
            if (runs) {
                const spec = corespec.of(core);
                const g = for (groups.items) |*g| {
                    if (g.spec.eql(spec)) break g;
                } else blk: {
                    groups.append(b.allocator, .{ .spec = spec }) catch @panic("OOM");
                    break :blk &groups.items[groups.items.len - 1];
                };
                g.cores.append(b.allocator, name) catch @panic("OOM");
            }
        }
    }

    for (groups.items) |g| {
        var query = target.query;
        query.cpu_features_sub = .empty;
        if (g.spec.model) |name| {
            query.cpu_model = .{ .explicit = findModel(arch, name) };
            query.cpu_features_add = .empty;
        } else {
            query.cpu_model = .baseline;
            query.cpu_features_add = featureSet(arch, g.spec.features);
        }

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
                .target = b.resolveTargetQuery(query),
                // Kernels are always built for speed, whatever the consumer uses.
                .optimize = .fast,
                .imports = &.{.{ .name = "blaskernel_options", .module = opts.createModule() }},
            }),
            .use_llvm = llvm,
            // the self-hosted linker writes the object otherwise, and pads and
            // misplaces SPARC symbols (calls land in zeroed `.text`)
            .use_lld = llvm,
        });
        mod.addObject(obj);
    }
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

fn findModel(arch: std.Target.Cpu.Arch, name: []const u8) *const std.Target.Cpu.Model {
    for (arch.allCpuModels()) |m| {
        if (std.mem.eql(u8, m.name, name)) return m;
    }
    std.debug.panic("no CPU model '{s}' for {s}", .{ name, @tagName(arch) });
}

fn featureSet(arch: std.Target.Cpu.Arch, names: []const []const u8) std.Target.Cpu.Feature.Set {
    var set: std.Target.Cpu.Feature.Set = .empty;
    const all = arch.allFeaturesList();
    for (names) |name| {
        const i = for (all, 0..) |feat, k| {
            if (std.mem.eql(u8, feat.name, name)) break k;
        } else std.debug.panic("no CPU feature '{s}' for {s}", .{ name, @tagName(arch) });
        set.addFeature(@intCast(i));
    }
    return set;
}
