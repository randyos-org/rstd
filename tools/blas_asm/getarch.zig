//! Reads `getarch.c`: what a core's `config.h` defines when `Makefile.system`
//! builds its `DYNAMIC_ARCH` kernel (`-DFORCE_<entry>`). `gen.zig` takes the
//! `-D` macros for each core's sources, `params.zig` for its `param.h` probe.

const std = @import("std");
const map = @import("openblas.zig");

/// An entry's `ARCHCONFIG` `-D` macros (subarch such as `PENTIUM3`, cache
/// sizes, `HAVE_*`, the DTB), then `CORE_<CORENAME>`.
pub const Config = struct {
    defines: []const []const u8,
    /// `DTB_DEFAULT_ENTRIES`, where `ARCHCONFIG` sets it.
    dtb: ?i64,
};

var configs: std.StringHashMapUnmanaged(?Config) = .empty;

/// `entry`'s config, preprocessing `getarch.c` for `target` (the triple that
/// takes the architecture's branch) with `extra` (`-DNO_AVX`, ...); null
/// where the entry sets no `ARCHCONFIG`. Preprocessing reads any entry for
/// any target. Each is preprocessed once.
pub fn config(target: []const u8, entry: []const u8, extra: []const []const u8) !?Config {
    const a = map.arena;
    var key: std.ArrayList(u8) = .empty;
    try key.print(a, "{s} {s}", .{ target, entry });
    for (extra) |e| try key.print(a, " {s}", .{e});
    if (configs.get(key.items)) |c| return c;

    var argv: std.ArrayList([]const u8) = .empty;
    try argv.appendSlice(a, &.{ "zig", "cc", "-E", "-dM", "-w", "-DDYNAMIC_ARCH", try a.print("-DFORCE_{s}", .{entry}) });
    try argv.appendSlice(a, extra);
    try argv.appendSlice(a, &.{ "-target", target, "-I", map.root, try map.path(&.{ map.root, "getarch.c" }) });
    const cc = try std.process.run(a, map.io, .{ .argv = argv.items });
    const ok = switch (cc.term) {
        .exited => |c| c == 0,
        else => false,
    };
    if (!ok or cc.stderr.len > 0) map.die("preprocessing getarch.c for {s} ({s}) failed:\n{s}", .{ entry, target, cc.stderr });

    const found: ?Config = if (try macroStrings(cc.stdout, "ARCHCONFIG")) |archconfig| blk: {
        const corename = try macroStrings(cc.stdout, "CORENAME") orelse map.die("FORCE_{s} has no CORENAME", .{entry});
        var defines: std.ArrayList([]const u8) = .empty;
        var dtb: ?i64 = null;
        var toks = std.mem.tokenizeAny(u8, archconfig, " \t");
        while (toks.next()) |t| {
            if (!std.mem.startsWith(u8, t, "-D")) continue;
            try defines.append(a, t);
            const dtb_key = "-DDTB_DEFAULT_ENTRIES=";
            if (std.mem.startsWith(u8, t, dtb_key)) dtb = try std.fmt.parseInt(i64, t[dtb_key.len..], 10);
        }
        try defines.append(a, try a.print("-DCORE_{s}", .{corename}));
        break :blk .{ .defines = defines.items, .dtb = dtb };
    } else null;
    try configs.put(a, key.items, found);
    return found;
}

/// The concatenated string literals of `#define <name> "..." "..."` in `-dM` output.
fn macroStrings(dump: []const u8, name: []const u8) !?[]const u8 {
    var lines = std.mem.splitScalar(u8, dump, '\n');
    while (lines.next()) |raw| {
        const line = std.mem.trimEnd(u8, raw, " \r");
        if (!std.mem.startsWith(u8, line, "#define ")) continue;
        const rest = line["#define ".len..];
        if (!std.mem.startsWith(u8, rest, name) or rest.len == name.len or rest[name.len] != ' ') continue;
        var out: std.ArrayList(u8) = .empty;
        var it = std.mem.splitScalar(u8, rest[name.len..], '"');
        var i: usize = 0;
        while (it.next()) |piece| : (i += 1) {
            if (i % 2 == 1) try out.appendSlice(map.arena, piece);
        }
        return out.items;
    }
    return null;
}

/// The entry (`FORCE_<it>`) whose `CORENAME` is `corename`, where exactly one
/// is: `P5`'s is `FORCE_P2`, the Pentium II target.
pub fn forceOf(corename: []const u8) !?[]const u8 {
    const text = map.read(try map.path(&.{ map.root, "getarch.c" })) orelse map.die("cannot read getarch.c", .{});
    var force: ?[]const u8 = null;
    var found: ?[]const u8 = null;
    var lines = std.mem.splitScalar(u8, text, '\n');
    while (lines.next()) |raw| {
        const line = std.mem.trim(u8, raw, " \t\r");
        if (std.mem.startsWith(u8, line, "#ifdef FORCE_")) {
            force = std.mem.trim(u8, line["#ifdef FORCE_".len..], " \t");
        } else if (std.mem.startsWith(u8, line, "#define CORENAME")) {
            const q = std.mem.trim(u8, line["#define CORENAME".len..], " \t\"");
            if (force == null or !std.mem.eql(u8, q, corename)) continue;
            if (found != null) return null;
            found = force;
        }
    }
    return found;
}

/// `Makefile.system`'s "fallbacks for 32bit": a `BINARY=32` build of the first
/// `TARGET` runs `getarch.c` for the second, so it gets that core's `config.h`,
/// `param.h` set and kernels. Its `TARGET_CORE` list sends `ZEN` to `NEHALEM`
/// instead; `dynamic.c` and `kernel/x86/KERNEL.ZEN` follow `TARGET`'s.
const binary32_targets = [_][2][]const u8{
    .{ "HASWELL", "NEHALEM" },        .{ "SKYLAKEX", "NEHALEM" },      .{ "COOPERLAKE", "NEHALEM" },
    .{ "SAPPHIRERAPIDS", "NEHALEM" }, .{ "SANDYBRIDGE", "NEHALEM" },   .{ "BULLDOZER", "BARCELONA" },
    .{ "PILEDRIVER", "BARCELONA" },   .{ "STEAMROLLER", "BARCELONA" }, .{ "EXCAVATOR", "BARCELONA" },
    .{ "ZEN", "BARCELONA" },          .{ "ARMV8", "ARMV7" },           .{ "POWER8", "POWER6" },
    .{ "POWER9", "POWER6" },          .{ "POWER10", "POWER6" },
};

/// The entry a 32-bit build of the upper-case target `core` takes.
pub fn target32(core: []const u8) []const u8 {
    for (binary32_targets) |t| {
        if (std.mem.eql(u8, t[0], core)) return t[1];
    }
    return core;
}
