//! The per-core kernels that are assembly upstream, for any routine at any
//! BLAS level. `index` is one architecture's generated `asm/__root__.zig` for
//! a routine: which converted file has which kernel, each file saying which
//! cores (and build variants) it serves. Each kernel is the source's own
//! function, naked, so it is exported and installed as is. Also each
//! architecture's registry (`Arch`), and where the `@Vector` kernels go (`Vector`).

const std = @import("std");
const sysinfo = @import("builtin");
const util = @import("util.zig");
const Vtable = @import("Vtable.zig");
const abi = @import("abi.zig");
const select = @import("cpu/select.zig");
const coretypes = @import("cpu/coretypes.zig");
const vector = @import("vector.zig");
const corespec = @import("cpu/corespec.zig");
// this file, for names `Arch` shares
const registry = @This();

const Slot = util.Slot;

/// The build variant (`tools/blas_asm/arch/<arch>.zig`) the target's calling
/// convention needs (`abi.zig`). The compiled target's, so `has`, which the
/// build script asks about the target it is building for, does not depend on it.
pub const variant: []const u8 = abi.variantName(&sysinfo.target);

/// The cores whose kernel objects this build links, and the one standing in
/// for a detected core that is not linked (`build/blaskern.zig`).
const linked = @import("blaskern_linked");

/// Whether `core`'s kernel object is linked.
pub fn isLinked(core: anytype) bool {
    @setEvalBranchQuota(100_000);
    for (linked.cores) |name| {
        if (std.mem.eql(u8, name, @tagName(core))) return true;
    }
    return false;
}

/// The linked objects as `select.served` reads them: which cores have one, and
/// the build's pick for the target as the last resort.
pub fn LinkedCtx(comptime Core: type) type {
    return struct {
        pub const arch = sysinfo.target.cpu.arch;
        pub const fallback: ?Core = if (linked.fallback) |f| std.meta.stringToEnum(Core, f) else null;
        pub fn ok(_: @This(), core: Core) bool {
            return isLinked(core);
        }
    };
}

/// The core whose object serves `core`: itself if linked, else the closest linked
/// core that needs nothing `core` does not (`select.served`), else the build's
/// pick for the target. Null means the portable kernels.
pub fn served(core: anytype) ?@TypeOf(core) {
    const Core = @TypeOf(core);
    return select.served(Core, core, LinkedCtx(Core){});
}

/// Points the vtable at `core`'s kernels from `modules` (an arch's list in
/// `kernels/<arch>.zig`), or at those of the core that serves it (`served`).
/// With none, the portable kernels stay.
pub fn install(vt: *Vtable, comptime modules: anytype, core: anytype) void {
    switch (served(core) orelse return) {
        inline else => |c| {
            if (comptime !isLinked(c)) return;
            inline for (modules) |m| m.install(vt, c);
        },
    }
}

/// An architecture's kernel registry over `modules` (each with `has`,
/// `exportCore` and `install`), for the targets in `serves`. With
/// `by_object`, a core another's object serves has no kernels of its own
/// (`corespec.objectCore`; 32-bit x86).
pub fn Arch(comptime C: type, comptime serves: []const std.Target.Cpu.Arch, comptime modules: anytype, comptime by_object: bool) type {
    return struct {
        pub const Core = C;

        /// Whether this registry has kernels for `arch`.
        pub fn supports(arch: std.Target.Cpu.Arch) bool {
            return std.mem.findScalar(std.Target.Cpu.Arch, serves, arch) != null;
        }

        /// Whether any routine has a kernel for `core`.
        pub fn hasKernels(comptime arch: std.Target.Cpu.Arch, comptime core: C) bool {
            if (comptime by_object and corespec.objectCore(arch, core) != core) return false;
            inline for (modules) |m| {
                if (comptime m.has(core)) return true;
            }
            return false;
        }

        /// Exports every kernel `core` has.
        pub fn exportCore(comptime core: C) void {
            inline for (modules) |m| m.exportCore(core);
        }

        /// Points the vtable at the kernels `core` has (`registry.install`);
        /// with no modules, the portable ones stay.
        pub fn install(vt: *Vtable, core: C) void {
            if (comptime modules.len == 0) return;
            registry.install(vt, modules, core);
        }
    };
}

/// One row of a generated `asm/__root__.zig`: a converted file, and the slot
/// its kernel fills.
pub const Entry = struct {
    slot: Slot,
    module: type,

    /// Whether any build of `core` uses this file.
    pub fn serves(comptime e: Entry, comptime core: []const u8) bool {
        inline for (e.module.keys) |key| {
            if (comptime std.mem.eql(u8, @tagName(key.core), core)) return true;
        }
        return false;
    }

    /// The file's function for `core` built as `v`, if it has one; a merged
    /// file's takes the key's `parts`.
    pub fn kernel(comptime e: Entry, comptime core: []const u8, comptime v: []const u8) ?type {
        inline for (e.module.keys) |key| {
            if (comptime std.mem.eql(u8, @tagName(key.core), core) and std.mem.eql(u8, key.variant, v)) {
                return if (@hasField(@TypeOf(key), "parts")) e.module.Body(key.parts) else e.module.Body;
            }
        }
        return null;
    }
};

/// A kernel upstream builds that the generator does not convert, and the
/// cores it builds it for: a row of an `asm/__root__.zig`'s `skipped`
/// (`tools/blas_asm/roles.zig`'s `skip`).
pub const Skipped = struct { slot: Slot, cores: []const @EnumLiteral() };

/// Exports `f` as `core`'s kernel for `slot`, under `util.kernelSymbol`.
pub fn exportKernel(comptime f: anytype, comptime slot: Slot, comptime core: anytype) void {
    exportNamed(f, util.kernelSymbol(slot, core));
}

/// The kernel `core`'s object exports for `slot` (`exportKernel`).
pub fn externKernel(comptime slot: Slot, comptime core: anytype) util.SlotFn(slot) {
    return externNamed(slot, util.kernelSymbol(slot, core));
}

/// Exports `f` as `name`, hidden in the per-core object.
fn exportNamed(comptime f: anytype, comptime name: []const u8) void {
    @export(f, .{ .name = name, .visibility = .hidden });
}

/// The `slot` kernel a per-core object exports as `name`.
fn externNamed(comptime slot: Slot, comptime name: []const u8) util.SlotFn(slot) {
    return @extern(util.SlotFn(slot), .{ .name = name, .visibility = .hidden });
}

/// Points each of `slots` that `core` has (`hasSlot`) at the kernel the
/// per-core object exports for it, for a routine whose kernels are not all
/// converted assembly (`Kernels.install` does this from `index`).
pub fn installSlots(vt: *Vtable, comptime core: anytype, comptime slots: []const Slot, comptime hasSlot: anytype) void {
    inline for (slots) |slot| {
        if (comptime hasSlot(core, slot)) @field(vt, @tagName(slot)) = externKernel(slot, core);
    }
}

/// What a hand port's kernel is written in: `asm` blocks, `@Vector`s, or
/// plain Zig.
pub const Written = enum { assembly, vector, scalar };

/// A hand-written kernel module a routine's `<arch>/kernel.zig` composes
/// beside its generated indexes: on `cores`, `module.entry(slot)` fills
/// `slots`; with `arch`, only in that OpenBLAS architecture's build
/// (`kernel/<arch>`). `tools/blas_asm/survey.zig` reads a `ports` table as
/// text, so keep its fields literal.
pub fn Port(comptime Core: type) type {
    return struct {
        module: type,
        cores: []const Core,
        slots: []const Slot,
        written: Written,
        arch: ?[]const u8 = null,
    };
}

/// The kernels of `ports` in OpenBLAS architecture `arch`'s build.
pub fn Ports(comptime Core: type, comptime arch: []const u8, comptime ports: []const Port(Core)) type {
    return struct {
        fn find(comptime core: Core, comptime slot: Slot) ?Port(Core) {
            for (ports) |p| {
                if (p.arch) |a| if (!std.mem.eql(u8, a, arch)) continue;
                if (std.mem.indexOfScalar(Core, p.cores, core) != null and std.mem.indexOfScalar(Slot, p.slots, slot) != null) return p;
            }
            return null;
        }

        /// Whether a port fills `slot` on `core`.
        pub fn hasSlot(comptime core: Core, comptime slot: Slot) bool {
            return comptime find(core, slot) != null;
        }

        /// Whether a port fills any slot on `core`.
        pub fn has(comptime core: Core) bool {
            for (ports) |p| {
                if (p.arch) |a| if (!std.mem.eql(u8, a, arch)) continue;
                if (std.mem.indexOfScalar(Core, p.cores, core) != null) return true;
            }
            return false;
        }

        /// Exports the port's kernel for `slot` on `core`, if one fills it.
        pub fn exportOne(comptime slot: Slot, comptime core: Core) void {
            if (comptime find(core, slot)) |p| exportKernel(p.module.entry(slot), slot, core);
        }

        /// Exports every kernel the ports have on `core`, for a routine
        /// whose kernels are all ports.
        pub fn exportCore(comptime core: Core) void {
            inline for (ports, 0..) |p, i| inline for (p.slots) |slot| {
                // once a slot, from the first port that lists it
                if (comptime firstWith(slot) == i) exportOne(slot, core);
            };
        }

        fn firstWith(comptime slot: Slot) usize {
            for (ports, 0..) |p, i| {
                if (std.mem.indexOfScalar(Slot, p.slots, slot) != null) return i;
            }
            unreachable;
        }

        /// Points the vtable at every kernel the ports have on `core`.
        pub fn install(vt: *Vtable, comptime core: Core) void {
            inline for (ports) |p| installSlots(vt, core, p.slots, hasSlot);
        }
    };
}

/// Routine modules side by side (each with `has`, `hasSlot`, `exportCore`
/// and `install`; the first with `skips`), as a routine's converted assembly
/// and its ports: each exports and installs its own, in order.
pub fn Union(comptime Core: type, comptime parts: anytype) type {
    return struct {
        pub fn hasSlot(comptime core: Core, comptime slot: Slot) bool {
            inline for (parts) |m| {
                if (comptime m.hasSlot(core, slot)) return true;
            }
            return false;
        }

        pub fn has(comptime core: Core) bool {
            inline for (parts) |m| {
                if (comptime m.has(core)) return true;
            }
            return false;
        }

        pub fn skips(comptime core: Core, comptime slot: Slot) bool {
            return parts[0].skips(core, slot);
        }

        pub fn exportCore(comptime core: Core) void {
            inline for (parts) |m| m.exportCore(core);
        }

        pub fn install(vt: *Vtable, comptime core: Core) void {
            inline for (parts) |m| m.install(vt, core);
        }
    };
}

pub fn Kernels(comptime Core: type, comptime index: type) type {
    return struct {
        /// `e`'s function for `core` on this target, if it has one and the
        /// `@Vector` kernels do not take its slot (`displaces`).
        fn kernel(comptime e: Entry, comptime core: Core) ?type {
            if (comptime vector.routed(e.slot)) @compileError(@tagName(e.slot) ++ " is routed, yet converted; tools/blas_asm/roles.zig's skip and vector.routed disagree");
            if (comptime displaces(core, e.slot)) return null;
            return e.kernel(@tagName(core), variant);
        }

        /// Whether any of these kernels serves `core`, built any way.
        pub fn has(comptime core: Core) bool {
            @setEvalBranchQuota(1_000_000);
            inline for (index.entries) |e| {
                if (comptime e.serves(@tagName(core))) return true;
            }
            return false;
        }

        /// Whether `slot` has a kernel on `core` for this target.
        pub fn hasSlot(comptime core: Core, comptime slot: Slot) bool {
            inline for (index.entries) |e| {
                if (comptime e.slot == slot and kernel(e, core) != null) return true;
            }
            return false;
        }

        /// Whether upstream builds `slot` on `core` from assembly the
        /// generator left out (`index.skipped`).
        pub fn skips(comptime core: Core, comptime slot: Slot) bool {
            if (!@hasDecl(index, "skipped")) return false;
            inline for (index.skipped) |s| {
                if (s.slot != slot) continue;
                inline for (s.cores) |c| {
                    if (comptime std.mem.eql(u8, @tagName(c), @tagName(core))) return true;
                }
            }
            return false;
        }

        /// Exports the kernel `slot` has on `core`, if it has one.
        pub fn exportOne(comptime slot: Slot, comptime core: Core) void {
            inline for (index.entries) |e| {
                if (comptime e.slot == slot) {
                    if (comptime kernel(e, core)) |K| exportKernel(&K.entry, e.slot, core);
                }
            }
        }

        /// Exports every kernel `core` has. Call from a container-level `comptime`
        /// block of the object's root file.
        pub fn exportCore(comptime core: Core) void {
            inline for (index.entries) |e| {
                if (comptime kernel(e, core)) |K| {
                    exportKernel(&K.entry, e.slot, core);
                } else if (comptime e.serves(@tagName(core)) and !hasSlot(core, e.slot) and !displaces(core, e.slot)) {
                    @compileError("no '" ++ variant ++ "' build of " ++ @tagName(e.slot) ++ " for " ++ @tagName(core) ++ "; tools/blas_asm/arch and abi.zig disagree");
                }
            }
        }

        /// Points the vtable at the kernels `core` has.
        pub fn install(vt: *Vtable, comptime core: Core) void {
            inline for (index.entries) |e| {
                if (comptime kernel(e, core) != null) @field(vt, @tagName(e.slot)) = externKernel(e.slot, core);
            }
        }
    };
}

// ---- where the `@Vector` kernels go -------------------------------------------

/// Where the `@Vector` kernels take slots from the converted ones beyond `vector.routed`
/// (`displaces`): every slot of theirs but `keep` on arm64's `vortex` and
/// `vortexm4`, where they mostly match or beat the converted kernels
/// (`tools/blasbench/results/results-m6.md`); `sum` is unmeasured there. `tools/blas_asm/survey.zig` reads it
/// as text.
pub const displacing: struct { cores: []const coretypes.Arm64CoreType, keep: []const Slot } = .{
    .cores = &.{ .vortex, .vortexm4 },
    .keep = &.{ .srot_k, .drot_k, .csrot_k, .zdrot_k, .ssum_k, .dsum_k, .csum_k, .zsum_k },
};

/// The cores whose vectors flush `f32` (`util.vectorIeeeOn` on their specs,
/// as a test checks): ARMv7's NEON ones. Their objects build no `s` or `c`
/// kernels (`vector.kernels`), so those take no such slot there.
/// `tools/blas_asm/survey.zig` reads it as text.
pub const flushing: struct { cores: []const coretypes.ArmCoreType } = .{
    .cores = &.{ .cortexa9, .cortexa15 },
};

/// Whether `slot` is over `f32` and `core` one of `flushing`'s.
fn flushes(comptime core: anytype, comptime slot: Slot) bool {
    if (@TypeOf(core) != coretypes.ArmCoreType) return false;
    if (std.mem.indexOfScalar(coretypes.ArmCoreType, flushing.cores, core) == null) return false;
    inline for (vector.every) |k| {
        if (k.slot == slot) return k.F == f32;
    }
    return false;
}

test "flushing names the arm cores whose specs flush f32 lanes" {
    // soft-float ABIs keep floats out of every vector unit
    if (sysinfo.target.abi.float() == .soft) return error.SkipZigTest;
    inline for (comptime corespec.kernelCores(.arm)) |core| {
        const listed = std.mem.indexOfScalar(coretypes.ArmCoreType, flushing.cores, core) != null;
        try std.testing.expectEqual(!util.vectorIeeeOn(f32, select.specCpu(core)), listed);
    }
}

/// Whether the `@Vector` kernels take `slot` from `core`'s converted ones, beyond
/// `vector.routed`, so its object leaves those out (`displacing`). On Apple's
/// systems Accelerate is to serve instead, and the converted kernels stay.
pub fn displaces(comptime core: anytype, comptime slot: Slot) bool {
    if (sysinfo.target.os.tag.isDarwin()) return false;
    if (@TypeOf(core) != coretypes.Arm64CoreType) return false;
    return std.mem.indexOfScalar(coretypes.Arm64CoreType, displacing.cores, core) != null and
        std.mem.indexOfScalar(Slot, displacing.keep, slot) == null and fills(slot);
}

/// Whether one of the `@Vector` kernels fills `slot` on this target.
pub fn fills(comptime slot: Slot) bool {
    inline for (vector.kernels) |k| {
        if (k.slot == slot) return true;
    }
    return false;
}

/// Where a registry puts the `@Vector` kernels (`vector.zig`), its `vector_rule`, which
/// `tools/blas_asm/survey.zig` reads as text.
pub fn Rule(comptime Core: type) type {
    return union(enum) {
        /// The `routed` slots, on `cores` or every core when null; a `peels`
        /// slot only where a converted module `skips` its assembly.
        routed: ?[]const Core,
        /// The `routed` slots; off Apple, whose Accelerate is to serve the
        /// rest, also those they `displace` and those no converted module
        /// fills (`hasSlot`), mostly C upstream, for the portable loop.
        unfilled,
        /// `unfilled` on `cores` only: those with a vector unit where others
        /// have none (riscv64). The rest take none of these kernels.
        unfilled_on: []const Core,
    };
}

/// Whether `rule` puts the `@Vector` kernels in `slot` on `core`, over the
/// registry's `converted` modules.
pub fn takes(comptime Core: type, comptime rule: Rule(Core), comptime converted: anytype, comptime core: Core, comptime slot: Slot) bool {
    if (flushes(core, slot)) return false;
    switch (rule) {
        .routed => |cores| {
            if (!vector.routed(slot)) return false;
            if (cores) |cs| if (std.mem.indexOfScalar(Core, cs, core) == null) return false;
            if (!vector.peels(slot)) return true;
            inline for (converted) |m| {
                if (@hasDecl(m, "skips") and m.skips(core, slot)) return true;
            }
            return false;
        },
        .unfilled => return unfilled(Core, converted, core, slot),
        .unfilled_on => |cores| return std.mem.indexOfScalar(Core, cores, core) != null and unfilled(Core, converted, core, slot),
    }
}

/// `Rule.unfilled`'s choice for `slot` on `core`.
fn unfilled(comptime Core: type, comptime converted: anytype, comptime core: Core, comptime slot: Slot) bool {
    if (vector.routed(slot)) return true;
    if (sysinfo.target.os.tag.isDarwin()) return false;
    if (displaces(core, slot)) return true;
    inline for (converted) |m| {
        if (@hasDecl(m, "hasSlot") and m.hasSlot(core, slot)) return false;
    }
    return true;
}

/// The `@Vector` kernels four vectors a step, in a per-core object for the
/// cores and slots `rule` picks (`takes`). Like the gemm kernels, they never
/// make a core need an object of its own (`has`).
pub fn Vector(comptime Core: type, comptime rule: Rule(Core), comptime converted: anytype) type {
    return struct {
        /// Whether these kernels take `slot` on `core`.
        pub fn uses(comptime core: Core, comptime slot: Slot) bool {
            return takes(Core, rule, converted, core, slot);
        }

        pub fn has(comptime core: Core) bool {
            _ = core;
            return false;
        }

        pub fn exportCore(comptime core: Core) void {
            inline for (vector.kernels) |k| {
                if (comptime uses(core, k.slot)) exportNamed(vector.kernel(k, 4), vector.symbol(k.slot, core));
            }
        }

        pub fn install(vt: *Vtable, comptime core: Core) void {
            inline for (vector.kernels) |k| {
                if (comptime uses(core, k.slot)) @field(vt, @tagName(k.slot)) = externNamed(k.slot, vector.symbol(k.slot, core));
            }
        }
    };
}
