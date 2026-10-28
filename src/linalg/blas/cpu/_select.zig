//! Choosing the core that serves a CPU, as pure decisions. A detector says which
//! core the CPU is (`_detect/<arch>.zig`); this finds the core to run when that
//! one has no kernel object in the build or the CPU cannot execute it. Light,
//! like `_corespec.zig`, so tests feed it synthetic values on any host.
//!
//! The reference does this with `#define gotoblas_X gotoblas_Y` chains
//! (`dynamic.c`), which follow no rule and are wrong for some builds (Barcelona
//! to Nehalem, SSE4.1 code on a CPU without SSSE3). Here the rule is the ISA: a
//! stand-in must need no feature the wanted core does not, and among those the
//! one with the most features wins.

const std = @import("std");
const coretypes = @import("_coretypes.zig");
const corespec = @import("_corespec.zig");

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

/// The core whose kernels serve `core`: itself if `linked` has its object, else the
/// closest one it covers, else `linked.fallback` (the build's pick for the
/// target, which the target CPU runs); null for the portable kernels. `linked`
/// has `ok(core)` (whether the object is linked) and a `fallback` declaration.
pub fn served(comptime Core: type, core: Core, linked: anytype) ?Core {
    if (linked.ok(core)) return core;
    return closest(Core, core, linked, .subset) orelse @TypeOf(linked).fallback;
}

/// The reference's static-only `TARGET`s that no detector returns: only
/// `OPENBLAS_CORETYPE` names them, never a fallback (their tuning is for
/// another class of CPU).
pub fn staticOnly(core: anytype) bool {
    return switch (@TypeOf(core)) {
        coretypes.X86CoreType => switch (core) {
            .p5, .p6, .viac3, .yonah => true,
            else => false,
        },
        coretypes.Arm64CoreType => switch (core) {
            .cortexa76, .cortexa510, .cortexa710, .cortexx1, .cortexx2, .ampere1, .ft2000 => true,
            else => false,
        },
        coretypes.RiscvCoreType => switch (core) {
            .c910v, .u74, .x280 => true,
            else => false,
        },
        coretypes.MipsCoreType => switch (core) {
            .mips24k, .mips1004k, .p5600, .i6400, .i6500, .p6600 => true,
            else => false,
        },
        else => false,
    };
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

/// The whole choice, as `gotoblas_dynamic_init` makes it, for a detector `D`
/// (`_detect/<arch>.zig`: `CoreType`, `Caps`, `decide`, `runnable`). `forced` is
/// `OPENBLAS_CORETYPE`. A name no core has is reported and detection decides, a
/// forced core the CPU cannot run is replaced by the closest one it can, and a
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
                const alt = closest(Core, c, RunnableCtx(D){ .caps = caps }, .any) orelse D.decide(caps).core;
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

// ---- tests -----------------------------------------------------------------

const X86 = coretypes.X86CoreType;

const TestCtx = struct {
    linked: []const X86 = &.{},
    cannot: []const X86 = &.{},

    pub fn ok(ctx: TestCtx, c: X86) bool {
        return std.mem.findScalar(X86, ctx.linked, c) != null and
            std.mem.findScalar(X86, ctx.cannot, c) == null;
    }
};

test "subsetOf follows the ISA, not the order of the enum" {
    try std.testing.expect(subsetOf(X86, .prescott, .haswell));
    try std.testing.expect(subsetOf(X86, .haswell, .skylakex));
    try std.testing.expect(subsetOf(X86, .skylakex, .cooperlake));
    try std.testing.expect(subsetOf(X86, .cooperlake, .sapphirerapids));
    try std.testing.expect(!subsetOf(X86, .haswell, .sandybridge));
    // Barcelona has SSE3 only; Nehalem needs SSSE3 and SSE4.1
    try std.testing.expect(!subsetOf(X86, .nehalem, .barcelona));
    // FMA4 is not in the Intel cores
    try std.testing.expect(!subsetOf(X86, .bulldozer, .haswell));
    try std.testing.expect(subsetOf(X86, .bulldozer, .piledriver));
}

test "closest takes the richest core the wanted one covers" {
    const all = [_]X86{ .prescott, .nehalem, .sandybridge, .haswell, .zen, .skylakex };
    try std.testing.expectEqual(@as(?X86, .skylakex), closest(X86, .skylakex, TestCtx{ .linked = &all }, .any));
    try std.testing.expectEqual(@as(?X86, .haswell), closest(X86, .cooperlake, TestCtx{ .linked = all[0..5] }, .any));
    try std.testing.expectEqual(@as(?X86, .sandybridge), closest(X86, .haswell, TestCtx{ .linked = all[0..3] }, .any));
    try std.testing.expectEqual(@as(?X86, .prescott), closest(X86, .barcelona, TestCtx{ .linked = all[0..4] }, .any));
    // nothing linked, nothing to run
    try std.testing.expectEqual(@as(?X86, null), closest(X86, .haswell, TestCtx{}, .any));
}

test "closest skips a core the CPU cannot run" {
    const linked = [_]X86{ .prescott, .nehalem, .sandybridge, .haswell };
    const ctx: TestCtx = .{ .linked = &linked, .cannot = &.{ .haswell, .sandybridge } };
    try std.testing.expectEqual(@as(?X86, .nehalem), closest(X86, .haswell, ctx, .any));
}

test "closest falls back to the richest core when none is covered" {
    // FMA4 cores: the only linked ones are Intel's; the richer of those runs
    const linked = [_]X86{ .prescott, .sandybridge };
    try std.testing.expectEqual(@as(?X86, .sandybridge), closest(X86, .bulldozer, TestCtx{ .linked = &linked }, .any));
}

test "byName is case insensitive and never names the placeholder" {
    try std.testing.expectEqual(@as(?X86, .sapphirerapids), byName(X86, "SapphireRapids"));
    try std.testing.expectEqual(@as(?X86, .opteron_sse3), byName(X86, "OPTERON_SSE3"));
    try std.testing.expectEqual(@as(?X86, null), byName(X86, "unknown"));
    try std.testing.expectEqual(@as(?X86, null), byName(X86, "haswell2"));
}

test "every core type has a feature set and the sets nest along the AVX-512 line" {
    inline for (.{
        coretypes.X86CoreType,       coretypes.Arm64CoreType, coretypes.ArmCoreType,
        coretypes.PowerCoreType,     coretypes.RiscvCoreType, coretypes.MipsCoreType,
        coretypes.LoongarchCoreType, coretypes.S390xCoreType, coretypes.SparcCoreType,
        coretypes.WasmCoreType,
    }) |C| _ = Sets(C).all;
}

const Log = struct {
    lines: [8][96]u8 = undefined,
    lens: [8]usize = @splat(0),
    n: usize = 0,
    max_level: i32 = 0,

    pub fn warn(l: *Log, level: i32, comptime fmt: []const u8, args: anytype) void {
        l.max_level = @max(l.max_level, level);
        if (l.n == l.lines.len) return;
        l.lens[l.n] = if (std.mem.print(&l.lines[l.n], fmt, args)) |s| s.len else |_| 0;
        l.n += 1;
    }
    fn has(l: *const Log, needle: []const u8) bool {
        for (0..l.n) |i| if (std.mem.find(u8, l.lines[i][0..l.lens[i]], needle) != null) return true;
        return false;
    }
};

fn LinkedList(comptime list: []const X86) type {
    return struct {
        pub const fallback: ?X86 = null;
        pub fn ok(_: @This(), c: X86) bool {
            return std.mem.findScalar(X86, list, c) != null;
        }
    };
}

const x86d = @import("_detect/x86.zig");

fn intelCaps(eax1: u32, level: u8) x86d.Caps {
    return .{ .vendor = .intel, .eax1 = eax1, .features = .init(.{
        .sse3 = true,
        .ssse3 = true,
        .sse4_1 = true,
        .avx = level >= 1,
        .fma = level >= 2,
        .avx2 = level >= 2,
        .avx512 = level >= 3,
        .avx512_bf16 = level >= 4,
        .amx_bf16 = level >= 5,
    }) };
}

test "choose: detection with every core linked is the detector's answer" {
    const all = comptime std.enums.values(X86)[1..];
    var log: Log = .{};
    const r = choose(x86d, intelCaps(0x306c3, 2), null, LinkedList(all){}, &log);
    try std.testing.expectEqual(X86.haswell, r.core);
    try std.testing.expectEqual(@as(?X86, .haswell), r.served);
    try std.testing.expectEqual(@as(usize, 0), log.n);
}

test "choose: a core with no object is served by the closest one it covers, and says so" {
    const list = [_]X86{ .prescott, .nehalem, .sandybridge, .haswell };
    var log: Log = .{};
    const r = choose(x86d, intelCaps(0x806f8, 5), null, LinkedList(&list){}, &log);
    try std.testing.expectEqual(X86.sapphirerapids, r.core);
    try std.testing.expectEqual(@as(?X86, .haswell), r.served);
    try std.testing.expect(log.has("using haswell"));
    try std.testing.expectEqual(@as(i32, 2), log.max_level);
}

test "choose: nothing linked means the portable kernels" {
    var log: Log = .{};
    const r = choose(x86d, intelCaps(0x306c3, 2), null, LinkedList(&.{}){}, &log);
    try std.testing.expectEqual(@as(?X86, null), r.served);
    try std.testing.expect(log.has("portable"));
}

test "choose: a forced core is honored when the CPU can run it" {
    const list = [_]X86{ .prescott, .nehalem, .sandybridge, .haswell };
    var log: Log = .{};
    const r = choose(x86d, intelCaps(0x306c3, 2), "NeHaLeM", LinkedList(&list){}, &log);
    try std.testing.expectEqual(X86.nehalem, r.core);
    try std.testing.expectEqual(@as(?X86, .nehalem), r.served);
    try std.testing.expectEqual(@as(usize, 0), log.n);
}

test "choose: a forced core the CPU cannot run is replaced, never executed" {
    const list = [_]X86{ .prescott, .nehalem, .sandybridge, .haswell, .skylakex };
    var log: Log = .{};
    // Haswell has no AVX-512
    const r = choose(x86d, intelCaps(0x306c3, 2), "skylakex", LinkedList(&list){}, &log);
    try std.testing.expectEqual(X86.haswell, r.core);
    try std.testing.expect(log.has("cannot run"));
    try std.testing.expectEqual(@as(i32, 1), log.max_level);
}

test "choose: a 64-bit build refuses the 32-bit-only cores by name" {
    const all = comptime std.enums.values(X86)[1..];
    for ([_][]const u8{ "katmai", "northwood", "p6", "yonah", "athlon" }) |name| {
        var log: Log = .{};
        var caps = intelCaps(0x306c3, 2);
        caps.is64 = true;
        const r = choose(x86d, caps, name, LinkedList(all){}, &log);
        try std.testing.expect(x86d.runnable(r.core, caps));
        try std.testing.expect(log.has("cannot run"));
    }
}

test "fallbacks never pick a static-only core; the user may name one" {
    const A = coretypes.Arm64CoreType;
    const Ctx = struct {
        pub fn ok(_: @This(), c: A) bool {
            return c != .neoversev1 and c != .armv8sve and c != .a64fx;
        }
    };
    // a V1 short of SVE support: Cortex-A510's and X2's sets are SVE subsets, but static
    const s = closest(A, .neoversev1, Ctx{}, .any).?;
    try std.testing.expect(!staticOnly(s));
    try std.testing.expectEqual(@as(?A, .cortexx2), closest(A, .cortexx2, struct {
        pub fn ok(_: @This(), _: A) bool {
            return true;
        }
    }{}, .any));
}

test "choose: an unknown name is reported and detection decides" {
    const list = [_]X86{ .prescott, .haswell };
    var log: Log = .{};
    const r = choose(x86d, intelCaps(0x306c3, 2), "haswel", LinkedList(&list){}, &log);
    try std.testing.expectEqual(X86.haswell, r.core);
    try std.testing.expect(log.has("Core not found: haswel"));
    // the placeholder is no name either
    var log2: Log = .{};
    _ = choose(x86d, intelCaps(0x306c3, 2), "unknown", LinkedList(&list){}, &log2);
    try std.testing.expect(log2.has("Core not found"));
}

/// The cores of `arch` with objects, as a `served` context.
fn Objects(comptime arch: std.Target.Cpu.Arch) type {
    const Core = corespec.CoreType(arch);
    return struct {
        pub const fallback: ?Core = null;
        pub fn ok(_: @This(), c: Core) bool {
            return std.mem.findScalar(Core, corespec.kernelCores(arch), c) != null;
        }
    };
}

/// Every core of `arch` (but the placeholder) is served by an object that needs
/// nothing it lacks, or by none at all.
fn checkServed(comptime arch: std.Target.Cpu.Arch, comptime must_have_object: bool) !void {
    const Core = corespec.CoreType(arch);
    for (std.enums.values(Core)) |c| {
        if (!selectable(c)) continue;
        const s = served(Core, c, Objects(arch){});
        if (s) |k| {
            try std.testing.expect(subsetOf(Core, k, c));
        } else if (must_have_object and Sets(Core).all[@backingInt(c)].count() != 0) {
            std.debug.print("{s}: {t} has no core to run it\n", .{ @tagName(arch), c });
            return error.TestUnexpectedResult;
        }
    }
}

test "every detectable core has an object that serves it" {
    try checkServed(.x86_64, true);
    try checkServed(.x86, true);
    try checkServed(.aarch64, true);
    try checkServed(.arm, true);
}

test "no core is ever served by one that needs more" {
    try checkServed(.powerpc64le, false);
    try checkServed(.riscv64, false);
    try checkServed(.s390x, false);
    try checkServed(.mips64, false);
    try checkServed(.loongarch64, false);
    try checkServed(.sparc64, false);
    try checkServed(.wasm32, false);
}
