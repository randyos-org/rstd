//! `select.zig`'s tests.

const std = @import("std");
const coretypes = @import("../coretypes.zig");
const corespec = @import("../corespec.zig");
const x86d = @import("../detect/x86.zig");
const select = @import("../select.zig");

const Sets = select.Sets;
const subsetOf = select.subsetOf;
const closest = select.closest;
const served = select.served;
const staticOnly = select.staticOnly;
const selectable = select.selectable;
const byName = select.byName;
const choose = select.choose;

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
    inline for (coretypes.all) |C| _ = Sets(C).all;
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
        pub const arch: std.Target.Cpu.Arch = .x86_64;
        pub const fallback: ?X86 = null;
        pub fn ok(_: @This(), c: X86) bool {
            return std.mem.findScalar(X86, list, c) != null;
        }
    };
}

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

/// The cores of `of` with objects, as a `served` context.
fn Objects(comptime of: std.Target.Cpu.Arch) type {
    const Core = corespec.CoreType(of);
    return struct {
        pub const arch = of;
        pub const fallback: ?Core = null;
        pub fn ok(_: @This(), c: Core) bool {
            return std.mem.findScalar(Core, corespec.kernelCores(of), c) != null;
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
