//! Extracts OpenBLAS's per-core gemm blocking parameters from `param.h` into
//! `src/linalg/blas/blas3/gemm/params_data.zig`.
//!
//!     OPENBLAS=<checkout> zig run tools/blas_asm/params.zig
//!
//! For each core it compiles and runs a small C program that includes `param.h`
//! with what the core's `config_kernel.h` defines (its `getarch.c` entry's
//! macros, see `Config`, and its architecture's), and prints the
//! values the header gives: per precision `P`, `Q`, `R` and the unroll sizes,
//! plus the offsets, alignment and the tuning knobs that land in the vtable. A
//! `P` or `R` that the header leaves to a run-time variable (`sgemm_p`,
//! `sgemm_r`, which `driver/others/parameter.c` fills in) prints as -1 and
//! becomes `.dynamic` / `.from_buffer`. The rules for the dynamic ones are in
//! `params.zig`, not here. Each table also takes the `BUFFER_SIZE` its core's
//! build sees in `common_<arch>.h`, preprocessed with the same macros.

const std = @import("std");
const Io = std.Io;
const map = @import("kernel_map.zig");

/// One architecture family: the header macros it needs and the cores it has.
const Family = struct {
    /// The generated function's name.
    name: []const u8,
    /// Enum type in `cpu/_coretypes.zig`.
    core_type: []const u8,
    /// What `c_check` puts in `config.h` for the target (the probe itself is
    /// built for the host): architecture, word size, byte order.
    defines: []const []const u8,
    /// The `zig cc -target` triple that makes `getarch.c` take this architecture's
    /// branch, with libc headers to include (wasm's block is not gated).
    target: []const u8,
    /// Cores with no `getarch.c` entry of their own: the entries to take, most
    /// specific last.
    renames: []const [2][]const u8 = &.{},
    /// What `common_<arch>.h` needs to preprocess for a core with no table.
    base: []const []const u8 = &.{},
};

const families = [_]Family{
    .{ .name = "x86_64", .core_type = "X86CoreType", .target = "x86_64-linux-gnu", .defines = &.{ "-DARCH_X86_64", "-DBINARY64", "-D__64BIT__" } },
    .{ .name = "x86", .core_type = "X86CoreType", .target = "x86-linux-gnu", .defines = &.{ "-DARCH_X86", "-DBINARY32", "-D__32BIT__" } },
    .{ .name = "arm64", .core_type = "Arm64CoreType", .target = "aarch64-linux-gnu", .defines = &.{ "-DARCH_ARM64", "-DBINARY64", "-D__64BIT__", "-DARMV8" } },
    .{ .name = "arm", .core_type = "ArmCoreType", .target = "arm-linux-gnueabihf", .defines = &.{ "-DARCH_ARM", "-DBINARY32", "-D__32BIT__" }, .renames = &.{ .{ "neon", "ARMV7" }, .{ "vfpv4", "ARMV7,CORTEXA15" } }, .base = &.{"-DARMV5"} },
    // `PPC970`'s unrolls follow the byte order
    .{ .name = "power64", .core_type = "PowerCoreType", .target = "powerpc64le-linux-gnu", .defines = &.{ "-DARCH_POWER", "-DBINARY64", "-D__64BIT__", "-DOS_LINUX" } },
    .{ .name = "power64be", .core_type = "PowerCoreType", .target = "powerpc64-linux-musl", .defines = &.{ "-DARCH_POWER", "-DBINARY64", "-D__64BIT__", "-DOS_LINUX", "-U__BYTE_ORDER__", "-D__BYTE_ORDER__=__ORDER_BIG_ENDIAN__" } },
    .{ .name = "power32", .core_type = "PowerCoreType", .target = "powerpc-linux-musleabihf", .defines = &.{ "-DARCH_POWER", "-DBINARY32", "-D__32BIT__", "-DOS_LINUX", "-U__BYTE_ORDER__", "-D__BYTE_ORDER__=__ORDER_BIG_ENDIAN__" } },
    .{ .name = "power32le", .core_type = "PowerCoreType", .target = "powerpc-linux-musleabihf", .defines = &.{ "-DARCH_POWER", "-DBINARY32", "-D__32BIT__", "-DOS_LINUX" } },
    .{ .name = "riscv64", .core_type = "RiscvCoreType", .target = "riscv64-linux-gnu", .defines = &.{ "-DARCH_RISCV64", "-DBINARY64", "-D__64BIT__" } },
    .{ .name = "s390x", .core_type = "S390xCoreType", .target = "s390x-linux-gnu", .defines = &.{ "-DARCH_ZARCH", "-DBINARY64", "-D__64BIT__" } },
    .{ .name = "mips64", .core_type = "MipsCoreType", .target = "mips64el-linux-gnuabi64", .defines = &.{ "-DARCH_MIPS64", "-DBINARY64", "-D__64BIT__" } },
    .{ .name = "mips32", .core_type = "MipsCoreType", .target = "mips-linux-gnueabihf", .defines = &.{ "-DARCH_MIPS", "-DBINARY32", "-D__32BIT__" } },
    .{ .name = "loongarch64", .core_type = "LoongarchCoreType", .target = "loongarch64-linux-gnu", .defines = &.{ "-DARCH_LOONGARCH64", "-DBINARY64", "-D__64BIT__" } },
    .{ .name = "sparc", .core_type = "SparcCoreType", .target = "sparc-linux-gnu", .defines = &.{ "-DARCH_SPARC", "-DBINARY32", "-D__32BIT__" } },
    .{ .name = "wasm", .core_type = "WasmCoreType", .target = "x86_64-linux-gnu", .defines = &.{ "-DARCH_WASM", "-DBINARY32", "-D__32BIT__" } },
    .{ .name = "generic", .core_type = "GenericCoreType", .target = "x86_64-linux-gnu", .defines = &.{ "-DARCH_X86_64", "-DBINARY64", "-D__64BIT__" } },
};

const probe_head =
    \\#include <stdio.h>
    \\typedef long BLASLONG;
    \\typedef unsigned long BLASULONG;
    \\long sgemm_q=-1,dgemm_q=-1,cgemm_q=-1,zgemm_q=-1,qgemm_q=-1,xgemm_q=-1;
    \\long shgemm_p=-1,shgemm_q=-1,shgemm_r=-1,bgemm_p=-1,bgemm_q=-1,bgemm_r=-1,sbgemm_p=-1,sbgemm_q=-1,sbgemm_r=-1;
    \\long cgemm3m_p=-1,cgemm3m_q=-1,cgemm3m_r=-1,zgemm3m_p=-1,zgemm3m_q=-1,zgemm3m_r=-1,xgemm3m_p=-1,xgemm3m_q=-1,xgemm3m_r=-1;
    \\long sgemm_p=-1,dgemm_p=-1,cgemm_p=-1,zgemm_p=-1,qgemm_p=-1,xgemm_p=-1;
    \\long sgemm_r=-1,dgemm_r=-1,cgemm_r=-1,zgemm_r=-1,qgemm_r=-1,xgemm_r=-1;
    \\#include "param.h"
    \\#define V(name, x) printf("%s=%ld\n", name, (long)(x))
    \\#define N(name) printf("%s=?\n", name)
    \\int main(void) {
    \\
;

/// A precision's field in `Table` and the head of its header macros.
const Prec = struct { field: []const u8, head: []const u8 };

const precs = [_]Prec{
    .{ .field = "s", .head = "SGEMM" },
    .{ .field = "d", .head = "DGEMM" },
    .{ .field = "q", .head = "QGEMM" },
    .{ .field = "c", .head = "CGEMM" },
    .{ .field = "z", .head = "ZGEMM" },
    .{ .field = "x", .head = "XGEMM" },
    .{ .field = "sh", .head = "SHGEMM" },
    .{ .field = "b", .head = "BGEMM" },
    .{ .field = "sb", .head = "SBGEMM" },
    .{ .field = "c3m", .head = "CGEMM3M" },
    .{ .field = "z3m", .head = "ZGEMM3M" },
    .{ .field = "x3m", .head = "XGEMM3M" },
};

/// The values one precision prints, keyed `<field>.<what>`.
fn precisionProbe(a: std.mem.Allocator, w: *std.ArrayList(u8), pr: Prec) !void {
    for ([_][]const u8{ "P", "Q", "R", "UNROLL_M", "UNROLL_N", "UNROLL_MN", "P_BASE" }) |what| {
        try w.print(a, "#ifdef {s}_DEFAULT_{s}\n  V(\"{s}.{s}\", {s}_DEFAULT_{s});\n#else\n  N(\"{s}.{s}\");\n#endif\n", .{ pr.head, what, pr.field, what, pr.head, what, pr.field, what });
    }
    try w.print(a, "#ifdef {s}_ALIGN_K\n  V(\"{s}.ALIGN_K\", {s}_ALIGN_K);\n#else\n  N(\"{s}.ALIGN_K\");\n#endif\n", .{ pr.head, pr.field, pr.head, pr.field });
}

const globals = [_][2][]const u8{
    .{ "offset_a", "GEMM_DEFAULT_OFFSET_A" },
    .{ "offset_b", "GEMM_DEFAULT_OFFSET_B" },
    .{ "align", "GEMM_DEFAULT_ALIGN" },
    .{ "switch_ratio", "SWITCH_RATIO" },
    .{ "divide_rate", "GEMM_DIVIDE_RATE" },
    .{ "divide_limit", "GEMM_DIVIDE_LIMIT" },
    .{ "preferred_size", "GEMM_PREFERRED_SIZE" },
    .{ "symv_p", "SYMV_P" },
};

/// `common_<arch>.h`, which `common.h` includes for `fam`'s `ARCH_<ARCH>`.
fn commonHeader(a: std.mem.Allocator, fam: Family) ![]const u8 {
    for (fam.defines) |d| {
        if (!std.mem.startsWith(u8, d, "-DARCH_")) continue;
        const name = try a.dupe(u8, d["-DARCH_".len..]);
        for (name) |*c| c.* = std.ascii.toLower(c.*);
        return a.print("common_{s}.h", .{name});
    }
    map.die("{s} has no ARCH_ macro", .{fam.name});
}

/// `BUFFER_SIZE` as `fam`'s `common_<arch>.h` defines it under `defines`:
/// preprocessed for the target, then its `(<n> << <shift>)` evaluated. Not
/// under `DYNAMIC_ARCH`: only `common_power.h:923` branches on it, to the
/// 256 MiB its dynamic cores (POWER6 to 10) have anyway; the other POWER
/// cores, never in a dynamic build, keep their static sizes.
fn bufferSize(a: std.mem.Allocator, io: std.Io, checkout: []const u8, fam: Family, defines: []const []const u8) !i64 {
    const src = ".zig-cache/params_buffer.c";
    try Io.Dir.cwd().writeFile(io, .{ .sub_path = src, .data = try a.print("#include \"{s}\"\nrstd_buffer_size BUFFER_SIZE\n", .{try commonHeader(a, fam)}) });
    var argv: std.ArrayList([]const u8) = .empty;
    try argv.appendSlice(a, &.{ "zig", "cc", "-E", "-P", "-w", "-target", fam.target, "-I", checkout });
    try argv.appendSlice(a, fam.defines);
    try argv.appendSlice(a, defines);
    try argv.append(a, src);
    const cc = try std.process.run(a, io, .{ .argv = argv.items });
    const ok = switch (cc.term) {
        .exited => |c| c == 0,
        else => false,
    };
    if (!ok) map.die("preprocessing {s} for {s} failed:\n{s}", .{ src, fam.name, cc.stderr });
    const key = "rstd_buffer_size ";
    const at = std.mem.find(u8, cc.stdout, key) orelse map.die("no BUFFER_SIZE for {s}", .{fam.name});
    const line = std.mem.sliceTo(cc.stdout[at + key.len ..], '\n');
    var it = std.mem.tokenizeAny(u8, line, "() \t\r");
    const n = it.next() orelse map.die("BUFFER_SIZE of {s}: {s}", .{ fam.name, line });
    const op = it.next() orelse map.die("BUFFER_SIZE of {s}: {s}", .{ fam.name, line });
    const shift = it.next() orelse map.die("BUFFER_SIZE of {s}: {s}", .{ fam.name, line });
    if (!std.mem.eql(u8, op, "<<") or it.next() != null) map.die("BUFFER_SIZE of {s}: {s}", .{ fam.name, line });
    const base = try std.fmt.parseInt(i64, std.mem.trimEnd(u8, n, "UuLl"), 10);
    return base << try std.fmt.parseInt(u6, shift, 10);
}

fn probeSource(a: std.mem.Allocator) ![]const u8 {
    var w: std.ArrayList(u8) = .empty;
    try w.appendSlice(a, probe_head);
    for (precs) |pr| try precisionProbe(a, &w, pr);
    for (globals) |g| {
        try w.print(a, "#ifdef {s}\n  V(\"{s}\", {s});\n#else\n  N(\"{s}\");\n#endif\n", .{ g[1], g[0], g[1], g[0] });
    }
    try w.appendSlice(a, "  return 0;\n}\n");
    return w.items;
}

/// The enum's field names, read from `cpu/_coretypes.zig`.
fn coreNames(a: std.mem.Allocator, text: []const u8, type_name: []const u8) ![]const []const u8 {
    const head = try a.print("pub const {s} = enum", .{type_name});
    const start = std.mem.find(u8, text, head) orelse map.die("no enum {s}", .{type_name});
    const open = std.mem.findScalarPos(u8, text, start, '{').?;
    const close = std.mem.findPos(u8, text, open, "};") orelse std.mem.findScalarPos(u8, text, open, '}').?;
    var out: std.ArrayList([]const u8) = .empty;
    var lines = std.mem.splitScalar(u8, text[open + 1 .. close], '\n');
    while (lines.next()) |raw| {
        var line = std.mem.trim(u8, raw, " \t\r");
        if (line.len == 0 or std.mem.startsWith(u8, line, "//")) continue;
        if (std.mem.find(u8, line, "//")) |c| line = std.mem.trim(u8, line[0..c], " \t");
        // `a, b, c,` on one line
        var it = std.mem.tokenizeAny(u8, line, ", ");
        while (it.next()) |name| {
            if (name.len > 0 and (std.ascii.isAlphabetic(name[0]) or name[0] == '@')) try out.append(a, try a.dupe(u8, name));
        }
    }
    return out.items;
}

const Values = std.array_hash_map.String(?i64);

fn parseValues(a: std.mem.Allocator, text: []const u8) !Values {
    var v: Values = .empty;
    var lines = std.mem.splitScalar(u8, text, '\n');
    while (lines.next()) |line| {
        const eq = std.mem.findScalar(u8, line, '=') orelse continue;
        const key = try a.dupe(u8, std.mem.trim(u8, line[0..eq], " \r"));
        const val = std.mem.trim(u8, line[eq + 1 ..], " \r");
        try v.put(a, key, if (std.mem.eql(u8, val, "?")) null else try std.fmt.parseInt(i64, val, 10));
    }
    return v;
}

fn upper(a: std.mem.Allocator, s: []const u8) ![]const u8 {
    const out = try a.dupe(u8, s);
    for (out) |*c| c.* = std.ascii.toUpper(c.*);
    return out;
}

fn size(a: std.mem.Allocator, w: *std.ArrayList(u8), v: ?i64, dynamic: []const u8) !void {
    if (v) |x| {
        if (x < 0) try w.appendSlice(a, dynamic) else try w.print(a, ".{{ .fixed = {d} }}", .{x});
    } else try w.appendSlice(a, dynamic);
}

/// Writes a switch's prongs, one per distinct value, with every core that has
/// it (in table order) on the prong's left side.
fn writeProngs(a: std.mem.Allocator, out: *std.ArrayList(u8), names: []const []const u8, values: []const []const u8) !void {
    for (values, 0..) |v, i| {
        const seen = for (values[0..i]) |p| {
            if (std.mem.eql(u8, p, v)) break true;
        } else false;
        if (seen) continue;
        try out.appendSlice(a, "        ");
        var first = true;
        for (names, values) |n, w| {
            if (!std.mem.eql(u8, w, v)) continue;
            if (!first) try out.appendSlice(a, ", ");
            try out.print(a, ".{s}", .{n});
            first = false;
        }
        try out.print(a, " => {s},\n", .{v});
    }
}

/// What `getarch.c` puts in a core's `config_kernel.h` when `Makefile.system`
/// builds its `DYNAMIC_ARCH` kernel (`-DFORCE_<core>`): `ARCHCONFIG`'s `-D`
/// macros (subarch such as `PENTIUM3`, cache sizes, `HAVE_*`, the DTB) plus
/// `CORE_<CORENAME>`. `param.h` keys some cores on these, not on `CORE_`. A
/// 32-bit build takes its target's fallback entry (`map.target32`) and adds
/// `NO_AVX`, `NO_AVX2` and `NO_AVX512`, as `Makefile.system` does.
const Config = struct { defines: []const []const u8, dtb: ?i64 };

fn binary32(fam: Family) bool {
    for (fam.defines) |d| {
        if (std.mem.eql(u8, d, "-DBINARY32")) return true;
    }
    return false;
}

/// The config of the first entry found among: the core's own name, `macro`'s
/// parts (last first), the entry whose `CORENAME` is the core's
/// (`forceOf`). Preprocessing is enough to read it, and works for any target.
fn coreConfig(a: std.mem.Allocator, io: std.Io, checkout: []const u8, getarch: []const u8, fam: Family, core: []const u8, macro: []const u8) !?Config {
    const b32 = binary32(fam);
    var cands: std.ArrayList([]const u8) = .empty;
    const own = try upper(a, core);
    try cands.append(a, if (b32) map.target32(own) else own);
    var parts = std.mem.splitBackwardsScalar(u8, macro, ',');
    while (parts.next()) |part| try cands.append(a, part);
    if (forceOf(getarch, cands.items[0])) |f| try cands.append(a, f);
    for (cands.items) |name| {
        var argv: std.ArrayList([]const u8) = .empty;
        try argv.appendSlice(a, &.{ "zig", "cc", "-E", "-dM", "-w", "-DDYNAMIC_ARCH", try a.print("-DFORCE_{s}", .{name}) });
        if (b32) try argv.appendSlice(a, &.{ "-DNO_AVX", "-DNO_AVX2", "-DNO_AVX512" });
        try argv.appendSlice(a, &.{ "-target", fam.target, "-I", checkout, try map.path(&.{ checkout, "getarch.c" }) });
        const cc = try std.process.run(a, io, .{ .argv = argv.items });
        const ok = switch (cc.term) {
            .exited => |c| c == 0,
            else => false,
        };
        if (!ok or cc.stderr.len > 0) map.die("preprocessing getarch.c for {s} ({s}) failed:\n{s}", .{ name, fam.name, cc.stderr });
        const archconfig = try macroStrings(a, cc.stdout, "ARCHCONFIG") orelse continue;
        const corename = try macroStrings(a, cc.stdout, "CORENAME") orelse map.die("FORCE_{s} has no CORENAME", .{name});
        var defines: std.ArrayList([]const u8) = .empty;
        var dtb: ?i64 = null;
        var toks = std.mem.tokenizeAny(u8, archconfig, " \t");
        while (toks.next()) |t| {
            if (!std.mem.startsWith(u8, t, "-D")) continue;
            try defines.append(a, t);
            const key = "-DDTB_DEFAULT_ENTRIES=";
            if (std.mem.startsWith(u8, t, key)) dtb = try std.fmt.parseInt(i64, t[key.len..], 10);
        }
        try defines.append(a, try a.print("-DCORE_{s}", .{corename}));
        return .{ .defines = defines.items, .dtb = dtb };
    }
    return null;
}

/// The concatenated string literals of `#define <name> "..." "..."` in `-dM` output.
fn macroStrings(a: std.mem.Allocator, dump: []const u8, name: []const u8) !?[]const u8 {
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
            if (i % 2 == 1) try out.appendSlice(a, piece);
        }
        return out.items;
    }
    return null;
}

/// The `getarch.c` entry (`FORCE_<it>`) whose `CORENAME` is `corename`, where
/// exactly one is: `P5`'s is `FORCE_P2`, the Pentium II target.
fn forceOf(getarch: []const u8, corename: []const u8) ?[]const u8 {
    var force: ?[]const u8 = null;
    var found: ?[]const u8 = null;
    var lines = std.mem.splitScalar(u8, getarch, '\n');
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

/// The `getarch.c` entries of `name`, its upper-case name unless renamed; also
/// its `param.h` macros where no entry is found.
fn coreMacro(a: std.mem.Allocator, fam: Family, name: []const u8) ![]const u8 {
    for (fam.renames) |r| {
        if (std.mem.eql(u8, r[0], name)) return r[1];
    }
    return upper(a, name);
}

pub fn main(init: std.process.Init) !void {
    const a = init.arena.allocator();
    const checkout = init.environ_map.get("OPENBLAS") orelse
        map.die("set OPENBLAS to the OpenBLAS checkout", .{});
    try map.init(a, init.io, checkout);

    const io = init.io;
    const types_text = try Io.Dir.cwd().readFileAlloc(io, "src/linalg/blas/cpu/_coretypes.zig", a, .limited(1 << 20));
    const getarch = map.read(try map.path(&.{ checkout, "getarch.c" })) orelse map.die("cannot read getarch.c", .{});
    const src_path = ".zig-cache/params_probe.c";
    try Io.Dir.cwd().writeFile(io, .{ .sub_path = src_path, .data = try probeSource(a) });

    var out: std.ArrayList(u8) = .empty;
    try out.appendSlice(a,
        \\//! Per-core gemm blocking parameters, extracted from OpenBLAS `param.h`
        \\//! (`BUFFER_SIZE` from `common_<arch>.h`). Generated by
        \\//! `tools/blas_asm/params.zig`; re-derive it rather than editing it.
        \\//! `.dynamic` and `.from_buffer` mark the values the header leaves to
        \\//! `driver/others/parameter.c` (see `params.zig`).
        \\
        \\const coretypes = @import("../../cpu/_coretypes.zig");
        \\const params = @import("params.zig");
        \\const Table = params.Table;
        \\const Shape = params.Shape;
        \\
        \\
    );

    for (families) |fam| {
        const names = try coreNames(a, types_text, fam.core_type);
        try out.print(a, "/// `BUFFER_SIZE` of a core with no table of its own.\npub const {s}_buffer_size: i64 = {d};\n\n", .{ fam.name, try bufferSize(a, io, checkout, fam, fam.base) });
        var dtb: std.ArrayList(u8) = .empty;
        defer dtb.deinit(a);
        try dtb.print(a, "/// `DTB_DEFAULT_ENTRIES`; null where `getarch.c` has no entry for the core.\npub fn {s}_dtb(core: coretypes.{s}) ?i32 {{\n    return switch (core) {{\n", .{ fam.name, fam.core_type });
        const macros = try a.alloc([]const u8, names.len);
        for (names, macros) |name, *m| m.* = try coreMacro(a, fam, name);
        const configs = try a.alloc(?Config, names.len);
        for (names, macros, configs) |name, macro, *c| c.* = try coreConfig(a, io, checkout, getarch, fam, name, macro);
        const dtb_vals = try a.alloc([]const u8, names.len);
        for (configs, dtb_vals) |c, *v| {
            v.* = if (c) |cfg| if (cfg.dtb) |n| try a.print("{d}", .{n}) else "null" else "null";
        }
        try writeProngs(a, &dtb, names, dtb_vals);
        try dtb.appendSlice(a, "    };\n}\n\n");
        try out.print(a, "pub fn {s}(core: coretypes.{s}) ?Table {{\n    return switch (core) {{\n", .{ fam.name, fam.core_type });
        const table_vals = try a.alloc([]const u8, names.len);
        for (names, macros, configs, table_vals) |name, macro, config, *tv| {
            var val: std.ArrayList(u8) = .empty;
            var argv: std.ArrayList([]const u8) = .empty;
            try argv.appendSlice(a, &.{ "zig", "cc", "-w", "-DDYNAMIC_ARCH", "-I", try map.path(&.{checkout}), "-I", ".zig-cache" });
            try argv.appendSlice(a, fam.defines);
            const core_defines = if (config) |c| c.defines else defs: {
                var list: std.ArrayList([]const u8) = .empty;
                var defs = std.mem.tokenizeScalar(u8, macro, ',');
                while (defs.next()) |d| try list.append(a, try a.print("-D{s}", .{d}));
                break :defs list.items;
            };
            try argv.appendSlice(a, core_defines);
            try argv.appendSlice(a, &.{ src_path, "-o", ".zig-cache/params_probe.exe" });
            const cc = try std.process.run(a, io, .{ .argv = argv.items });
            const built = switch (cc.term) {
                .exited => |c| c == 0,
                else => false,
            };
            if (!built) {
                tv.* = "null";
                std.debug.print("{s}/{s}: header does not build\n", .{ fam.name, name });
                continue;
            }
            const run = try std.process.run(a, io, .{ .argv = &.{".zig-cache/params_probe.exe"} });
            const vals = try parseValues(a, run.stdout);
            const have_p = vals.get("s.P") orelse null;
            const have_q = vals.get("s.Q") orelse null;
            if (have_p == null or have_q == null) {
                tv.* = "null";
                std.debug.print("{s}/{s}: header defines no blocking\n", .{ fam.name, name });
                continue;
            }
            try val.appendSlice(a, ".{\n");
            for (precs) |pr| {
                const key = struct {
                    fn k(al: std.mem.Allocator, f: []const u8, w: []const u8) ![]const u8 {
                        return al.print("{s}.{s}", .{ f, w });
                    }
                }.k;
                const um = (vals.get(try key(a, pr.field, "UNROLL_M")) orelse null);
                const pp = (vals.get(try key(a, pr.field, "P")) orelse null);
                if (um == null and pp == null) continue;
                try val.print(a, "            .{s} = Shape{{ .unroll_m = {d}, .unroll_n = {d}, .unroll_mn = {d}, .align_k = {d}, .p_base = {d}, .p = ", .{
                    pr.field,
                    um orelse 0,
                    (vals.get(try key(a, pr.field, "UNROLL_N")) orelse null) orelse 0,
                    (vals.get(try key(a, pr.field, "UNROLL_MN")) orelse null) orelse 0,
                    (vals.get(try key(a, pr.field, "ALIGN_K")) orelse null) orelse 0,
                    (vals.get(try key(a, pr.field, "P_BASE")) orelse null) orelse 0,
                });
                try size(a, &val, pp, ".dynamic");
                try val.appendSlice(a, ", .q = ");
                try size(a, &val, (vals.get(try key(a, pr.field, "Q")) orelse null), ".dynamic");
                try val.appendSlice(a, ", .r = ");
                try size(a, &val, (vals.get(try key(a, pr.field, "R")) orelse null), ".from_buffer");
                try val.appendSlice(a, " },\n");
            }
            for (globals) |g| {
                if ((vals.get(g[0]) orelse null)) |x| {
                    const field = if (std.mem.eql(u8, g[0], "align")) "align_" else g[0];
                    try val.print(a, "            .{s} = {d},\n", .{ field, x });
                }
            }
            try val.print(a, "            .buffer_size = {d},\n", .{try bufferSize(a, io, checkout, fam, core_defines)});
            try val.appendSlice(a, "        }");
            tv.* = val.items;
        }
        try writeProngs(a, &out, names, table_vals);
        try out.appendSlice(a, "    };\n}\n\n");
        try out.appendSlice(a, dtb.items);
    }

    const dest = "src/linalg/blas/blas3/gemm/params_data.zig";
    try Io.Dir.cwd().writeFile(io, .{ .sub_path = dest, .data = out.items });
    _ = try std.process.run(a, io, .{ .argv = &.{ "zig", "fmt", dest } });
    std.debug.print("wrote {s}\n", .{dest});
}
