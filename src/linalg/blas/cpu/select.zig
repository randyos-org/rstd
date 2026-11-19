//! Choosing the core that serves a CPU, as pure decisions. A detector says which
//! core the CPU is (`detect/<arch>.zig`); this finds the core to run when that
//! one has no kernel object in the build or the CPU cannot execute it. Light,
//! like `corespec.zig`, so tests feed it synthetic values on any host.
//!
//! The reference does this with `#define gotoblas_X gotoblas_Y` chains
//! (`dynamic.c`), which follow no rule and are wrong for some builds (Barcelona
//! to Nehalem, SSE4.1 code on a CPU without SSSE3). Here the rule is the ISA: a
//! stand-in must need no feature the wanted core does not, and among those the
//! one with the most features wins.

const std = @import("std");
const coretypes = @import("coretypes.zig");
const corespec = @import("corespec.zig");
const x86d = @import("detect/x86.zig");

const Set = std.Target.Cpu.Feature.Set;

/// A detector's answer: the core, and a warning for the user when the CPU or OS
/// fell short of the core the model table first picked (`FALLBACK_VERBOSE`).
pub fn Decision(comptime Core: type) type {
    return struct {
        core: Core,
        note: ?[]const u8 = null,
    };
}

/// The architecture whose feature names a core type's specs use.
fn archOf(comptime Core: type) std.Target.Cpu.Arch {
    return switch (Core) {
        coretypes.X86CoreType => .x86_64,
        coretypes.Arm64CoreType => .aarch64,
        coretypes.ArmCoreType => .arm,
        coretypes.PowerCoreType => .powerpc64le,
        coretypes.RiscvCoreType => .riscv64,
        coretypes.MipsCoreType => .mips64,
        coretypes.LoongarchCoreType => .loongarch64,
        coretypes.S390xCoreType => .s390x,
        coretypes.SparcCoreType => .sparc64,
        coretypes.WasmCoreType => .wasm32,
        coretypes.GenericCoreType => .x86_64,
        else => @compileError("no feature namespace for " ++ @typeName(Core)),
    };
}

/// `spec`'s ISA features and what they imply: the model's, then the ones
/// added, less the ones taken out.
fn specSet(comptime arch: std.Target.Cpu.Arch, comptime spec: corespec.Spec) Set {
    var set: Set = .empty;
    if (spec.model) |name| {
        for (arch.allCpuModels()) |m| {
            if (std.mem.eql(u8, m.name, name)) set.addFeatureSet(m.features);
        }
    }
    const all = arch.allFeaturesList();
    for (spec.features) |name| {
        for (all, 0..) |f, i| {
            if (std.mem.eql(u8, f.name, name)) set.addFeature(@intCast(i));
        }
    }
    set.populateDependencies(all);
    for (all, 0..) |f, i| {
        if (corespec.notIsa(f.name)) set.removeFeature(@intCast(i));
        for (spec.without) |name| if (std.mem.eql(u8, f.name, name)) set.removeFeature(@intCast(i));
    }
    return set;
}

/// Every core's feature set, indexed by the enum's value.
pub fn Sets(comptime Core: type) type {
    return struct {
        pub const all: [@typeInfo(Core).@"enum".field_names.len]Set = blk: {
            @setEvalBranchQuota(10_000_000);
            var out: [@typeInfo(Core).@"enum".field_names.len]Set = undefined;
            for (@typeInfo(Core).@"enum".field_names, 0..) |n, i| {
                out[i] = specSet(archOf(Core), corespec.of(@field(Core, n)));
            }
            break :blk out;
        };
    };
}

/// A CPU with `core`'s spec's features and no model's, to ask what its
/// object's vectors compute with (`util.vectorIeeeOn`).
pub fn specCpu(comptime core: anytype) std.Target.Cpu {
    return comptime blk: {
        @setEvalBranchQuota(1_000_000);
        const arch = archOf(@TypeOf(core));
        break :blk .{ .arch = arch, .model = std.Target.Cpu.Model.generic(arch), .features = specSet(arch, corespec.of(core)) };
    };
}

/// Whether `a`'s kernels need no feature `b`'s do not.
pub fn subsetOf(comptime Core: type, a: Core, b: Core) bool {
    const sets = Sets(Core).all;
    return sets[@backingInt(b)].isSuperSetOf(sets[@backingInt(a)]);
}

/// What `closest` may stand in: only a core needing no feature the wanted one
/// lacks (so it runs wherever that does), or, failing that, the richest accepted.
pub const Mode = enum { subset, any };

/// The core to run in place of `want`: among those `ctx.ok` accepts, one that
/// needs no feature `want` does not (the one with the most features), else, in
/// mode `any`, the one with the most features at all. Ties go to the earlier core
/// in the enum. Null if none qualifies.
pub fn closest(comptime Core: type, want: Core, ctx: anytype, mode: Mode) ?Core {
    const sets = Sets(Core).all;
    var best: ?Core = null;
    var best_sub = false;
    var best_count: usize = 0;
    for (std.enums.values(Core)) |c| {
        if (!ctx.ok(c)) continue;
        if (c == want) return c;
        if (staticOnly(c)) continue;
        const sub = subsetOf(Core, c, want);
        if (mode == .subset and !sub) continue;
        const n: usize = sets[@backingInt(c)].count();
        if (best == null or (sub and !best_sub) or (sub == best_sub and n > best_count)) {
            best = c;
            best_sub = sub;
            best_count = n;
        }
    }
    return best;
}

/// The core whose kernels serve `core`: the core whose object it takes
/// (`corespec.objectCore`) if `linked` has that object, else the closest one
/// that covers, else `linked.fallback` (the build's pick for the target, which
/// the target CPU runs); null for the portable kernels. `linked` has `ok(core)`
/// (whether the object is linked) and `fallback` and `arch` declarations.
pub fn served(comptime Core: type, core: Core, linked: anytype) ?Core {
    const want = corespec.objectCore(@TypeOf(linked).arch, core);
    if (linked.ok(want)) return want;
    return closest(Core, want, linked, .subset) orelse @TypeOf(linked).fallback;
}

/// The reference's static-only `TARGET`s that no detector returns: only
/// `OPENBLAS_CORETYPE` names them, never a fallback (their tuning is for
/// another class of CPU).
pub fn staticOnly(core: anytype) bool {
    return std.mem.findScalar(@TypeOf(core), coretypes.staticOnly(@TypeOf(core)), core) != null;
}

/// What `OPENBLAS_CORETYPE` can name: every core but the placeholder.
pub fn selectable(core: anytype) bool {
    return !std.mem.eql(u8, @tagName(core), "unknown");
}

/// `name` as a core, matched case-insensitively. Null for no such core and for
/// the placeholder.
pub fn byName(comptime Core: type, name: []const u8) ?Core {
    inline for (@typeInfo(Core).@"enum".field_names) |n| {
        if (std.ascii.eqlIgnoreCase(n, name)) {
            const c = @field(Core, n);
            return if (selectable(c)) c else null;
        }
    }
    return null;
}

/// `d`, unless the CPU cannot run its core: then the closest core it can (the
/// reference trusts its model tables; a table can name a core for a CPU whose OS
/// or hypervisor hides what the core's kernels use).
pub fn guard(comptime D: type, d: Decision(D.CoreType), caps: D.Caps) Decision(D.CoreType) {
    if (D.runnable(d.core, caps)) return d;
    // the detector's own line of descent first (same class, same tuning)
    if (@hasDecl(D, "parent")) {
        var p = D.parent(d.core);
        while (p) |k| : (p = D.parent(k)) {
            if (D.runnable(k, caps)) return .{ .core = k, .note = d.note };
        }
    }
    const stand = closest(D.CoreType, d.core, RunnableCtx(D){ .caps = caps }, .any) orelse return d;
    return .{ .core = stand, .note = d.note };
}

/// `D.probe()`, or under `hwcap_mask` (`LD_HWCAP_MASK`, all bits if null)
/// where `D` takes one (`probeMasked`, s390x).
pub fn probe(comptime D: type, hwcap_mask: ?u64) D.Caps {
    if (comptime @hasDecl(D, "probeMasked")) return D.probeMasked(hwcap_mask orelse ~@as(u64, 0));
    return D.probe();
}

/// The running CPU's core: `D`'s decision on its probe.
pub fn detect(comptime D: type, hwcap_mask: ?u64) D.CoreType {
    return D.decide(probe(D, hwcap_mask)).core;
}

/// Whether this machine runs `core`'s kernels, for tests: `D.runsHere`
/// where `D` has one (mips), else `D.runnable` on its probe.
pub fn canRun(comptime D: type, core: D.CoreType) bool {
    if (comptime @hasDecl(D, "runsHere")) return D.runsHere(core);
    return D.runnable(core, D.probe());
}

/// What `choose` decided: the core the CPU or the user wants, and the core
/// whose kernels serve it (null: the portable kernels).
pub fn Chosen(comptime Core: type) type {
    return struct { core: Core, served: ?Core };
}

/// The cores the CPU can run, for standing in for a forced one.
fn RunnableCtx(comptime D: type) type {
    return struct {
        caps: D.Caps,
        pub fn ok(ctx: @This(), c: D.CoreType) bool {
            return selectable(c) and D.runnable(c, ctx.caps);
        }
    };
}

/// `RunnableCtx`, within what the detected core `within` needs where `D.runnable`
/// cannot see the ISA level (`level_unseen`): there `runnable` passes an ARMv8.2
/// core on an ARMv8.0 CPU, but the detector's part table does not (T18).
fn ForcedCtx(comptime D: type) type {
    return struct {
        caps: D.Caps,
        within: D.CoreType,
        pub fn ok(ctx: @This(), c: D.CoreType) bool {
            if (!(RunnableCtx(D){ .caps = ctx.caps }).ok(c)) return false;
            if (comptime !@hasDecl(D, "level_unseen")) return true;
            return subsetOf(D.CoreType, c, ctx.within);
        }
    };
}

/// The whole choice, as `gotoblas_dynamic_init` makes it, for a detector `D`
/// (`detect/<arch>.zig`: `CoreType`, `Caps`, `decide`, `runnable`). `forced` is
/// `OPENBLAS_CORETYPE`. A name no core has is reported and detection decides, a
/// forced core the CPU cannot run is replaced by the closest one it can (`ForcedCtx`), and a
/// core with no kernels in the build is served by the closest one that has
/// (`served`). `log` has `warn(level, fmt, args)`.
pub fn choose(comptime D: type, caps: D.Caps, forced: ?[]const u8, linked: anytype, log: anytype) Chosen(D.CoreType) {
    const Core = D.CoreType;
    var d: Decision(Core) = undefined;
    var from_user = false;
    if (forced) |name| {
        if (byName(Core, name)) |c| {
            from_user = true;
            if (D.runnable(c, caps)) {
                d = .{ .core = c };
            } else {
                const det = D.decide(caps).core;
                const alt = closest(Core, c, ForcedCtx(D){ .caps = caps, .within = det }, .any) orelse det;
                log.warn(1, "core {s} cannot run on this CPU, using {s}", .{ @tagName(c), @tagName(alt) });
                d = .{ .core = alt };
            }
        } else {
            log.warn(1, "Core not found: {s}", .{name});
        }
    }
    if (!from_user) {
        d = D.decide(caps);
        if (d.note) |n| log.warn(1, "{s}", .{n});
    }
    const s = served(Core, d.core, linked);
    if (s == null) {
        log.warn(if (from_user) 1 else 2, "no kernels for {s} in this build, using the portable ones", .{@tagName(d.core)});
    } else if (s.? != d.core) {
        log.warn(if (from_user) 1 else 2, "no kernels for {s} in this build, using {s}", .{ @tagName(d.core), @tagName(s.?) });
    }
    return .{ .core = d.core, .served = s };
}

/// Test doubles for `choose`.
pub const testing = struct {
    /// A build with every core of `Core` linked and no fallback.
    pub fn AllLinked(comptime target: std.Target.Cpu.Arch, comptime Core: type) type {
        return struct {
            pub const arch: std.Target.Cpu.Arch = target;
            pub const fallback: ?Core = null;
            pub fn ok(_: @This(), _: Core) bool {
                return true;
            }
        };
    }

    /// A `log` that drops every warning.
    pub const NoLog = struct {
        pub fn warn(_: NoLog, _: i32, comptime _: []const u8, _: anytype) void {}
    };
};
