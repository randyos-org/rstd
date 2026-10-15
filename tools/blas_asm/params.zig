//! Extracts OpenBLAS's per-core gemm blocking parameters from `param.h` into
//! `src/linalg/blas/blas3/gemm/params_data.zig`.
//!
//!     OPENBLAS=<checkout> zig run tools/blas_asm/params.zig
//!
//! For each core it compiles and runs a small C program that includes `param.h`
//! with the core's macro (and its architecture's) defined, and prints the
//! values the header gives: per precision `P`, `Q`, `R` and the unroll sizes,
//! plus the offsets, alignment and the tuning knobs that land in the vtable. A
//! `P` or `R` that the header leaves to a run-time variable (`sgemm_p`,
//! `sgemm_r`, which `driver/others/parameter.c` fills in) prints as -1 and
//! becomes `.dynamic` / `.from_buffer`. The rules for the dynamic ones are in
//! `params.zig`, not here.

const std = @import("std");
const Io = std.Io;
const map = @import("kernel_map.zig");

/// One architecture family: the header macros it needs and the cores it has.
const Family = struct {
    /// The generated function's name.
    name: []const u8,
    /// Enum type in `cpu/_coretypes.zig`.
    core_type: []const u8,
    defines: []const []const u8,
    /// Cores whose macro is not just the upper-cased enum name.
    renames: []const [2][]const u8 = &.{},
};

/// The old x86 cores are `CORE_<NAME>` in `param.h`.
const x86_renames = [_][2][]const u8{
    .{ "katmai", "CORE_KATMAI" },
    .{ "coppermine", "CORE_COPPERMINE" },
    .{ "northwood", "CORE_NORTHWOOD" },
    .{ "prescott", "CORE_PRESCOTT" },
    .{ "banias", "CORE_BANIAS" },
    .{ "opteron_sse3", "OPTERON" },
};

const families = [_]Family{
    .{ .name = "x86_64", .core_type = "X86CoreType", .defines = &.{ "-DARCH_X86_64", "-DBINARY64" }, .renames = &x86_renames },
    .{ .name = "x86", .core_type = "X86CoreType", .defines = &.{ "-DARCH_X86", "-DBINARY32" }, .renames = &x86_renames },
    .{ .name = "arm64", .core_type = "Arm64CoreType", .defines = &.{ "-DARCH_ARM64", "-DBINARY64", "-DARMV8" } },
    .{ .name = "arm", .core_type = "ArmCoreType", .defines = &.{ "-DARCH_ARM", "-DBINARY32" }, .renames = &.{ .{ "neon", "ARMV7" }, .{ "vfpv4", "ARMV7,CORTEXA15" } } },
    .{ .name = "power64", .core_type = "PowerCoreType", .defines = &.{ "-DARCH_POWER", "-DBINARY64", "-DOS_LINUX" } },
    .{ .name = "power32", .core_type = "PowerCoreType", .defines = &.{ "-DARCH_POWER", "-DBINARY32", "-DOS_LINUX" } },
    .{ .name = "riscv64", .core_type = "RiscvCoreType", .defines = &.{ "-DARCH_RISCV64", "-DBINARY64" } },
    .{ .name = "s390x", .core_type = "S390xCoreType", .defines = &.{ "-DARCH_ZARCH", "-DBINARY64" } },
    .{ .name = "mips64", .core_type = "MipsCoreType", .defines = &.{ "-DARCH_MIPS64", "-DBINARY64" } },
    .{ .name = "mips32", .core_type = "MipsCoreType", .defines = &.{ "-DARCH_MIPS", "-DBINARY32" } },
    .{ .name = "loongarch64", .core_type = "LoongarchCoreType", .defines = &.{ "-DARCH_LOONGARCH64", "-DBINARY64" } },
    .{ .name = "sparc", .core_type = "SparcCoreType", .defines = &.{ "-DARCH_SPARC", "-DBINARY32" }, .renames = &.{ .{ "sparcv7", "SPARC,V7" }, .{ "sparc", "SPARC,V9" } } },
    .{ .name = "wasm", .core_type = "WasmCoreType", .defines = &.{ "-DARCH_WASM", "-DBINARY32" } },
    .{ .name = "generic", .core_type = "GenericCoreType", .defines = &.{ "-DARCH_X86_64", "-DBINARY64" }, .renames = &.{.{ "generic", "GENERIC" }} },
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
    .{ "dtb_entries", "DTB_DEFAULT_ENTRIES" },
    .{ "symv_p", "SYMV_P" },
};

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
    const head = try std.fmt.allocPrint(a, "pub const {s} = enum", .{type_name});
    const start = std.mem.indexOf(u8, text, head) orelse map.die("no enum {s}", .{type_name});
    const open = std.mem.indexOfScalarPos(u8, text, start, '{').?;
    const close = std.mem.indexOfPos(u8, text, open, "};") orelse std.mem.indexOfScalarPos(u8, text, open, '}').?;
    var out: std.ArrayList([]const u8) = .empty;
    var lines = std.mem.splitScalar(u8, text[open + 1 .. close], '\n');
    while (lines.next()) |raw| {
        var line = std.mem.trim(u8, raw, " \t\r");
        if (line.len == 0 or std.mem.startsWith(u8, line, "//")) continue;
        if (std.mem.indexOf(u8, line, "//")) |c| line = std.mem.trim(u8, line[0..c], " \t");
        // `a, b, c,` on one line
        var it = std.mem.tokenizeAny(u8, line, ", ");
        while (it.next()) |name| {
            if (name.len > 0 and (std.ascii.isAlphabetic(name[0]) or name[0] == '@')) try out.append(a, try a.dupe(u8, name));
        }
    }
    return out.items;
}

const Values = std.StringArrayHashMapUnmanaged(?i64);

fn parseValues(a: std.mem.Allocator, text: []const u8) !Values {
    var v: Values = .empty;
    var lines = std.mem.splitScalar(u8, text, '\n');
    while (lines.next()) |line| {
        const eq = std.mem.indexOfScalar(u8, line, '=') orelse continue;
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

pub fn main(init: std.process.Init) !void {
    const a = init.arena.allocator();
    const checkout = init.environ_map.get("OPENBLAS") orelse
        map.die("set OPENBLAS to the OpenBLAS checkout", .{});
    try map.init(a, init.io, checkout);

    const io = init.io;
    const types_text = try Io.Dir.cwd().readFileAlloc(io, "src/linalg/blas/cpu/_coretypes.zig", a, .limited(1 << 20));
    const src_path = ".zig-cache/params_probe.c";
    try Io.Dir.cwd().writeFile(io, .{ .sub_path = src_path, .data = try probeSource(a) });

    var out: std.ArrayList(u8) = .empty;
    try out.appendSlice(a,
        \\//! Per-core gemm blocking parameters, extracted from OpenBLAS `param.h`.
        \\//! Generated by `tools/blas_asm/params.zig`; re-derive it rather than
        \\//! editing it. `.dynamic` and `.from_buffer` mark the values the header
        \\//! leaves to `driver/others/parameter.c` (see `params.zig`).
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
        try out.print(a, "pub fn {s}(core: coretypes.{s}) ?Table {{\n    return switch (core) {{\n", .{ fam.name, fam.core_type });
        for (names) |name| {
            var macro = try upper(a, name);
            for (fam.renames) |r| {
                if (std.mem.eql(u8, r[0], name)) macro = r[1];
            }
            var argv: std.ArrayList([]const u8) = .empty;
            try argv.appendSlice(a, &.{ "zig", "cc", "-w", "-DDYNAMIC_ARCH", "-I", try map.path(&.{checkout}), "-I", ".zig-cache" });
            try argv.appendSlice(a, fam.defines);
            var defs = std.mem.tokenizeScalar(u8, macro, ',');
            while (defs.next()) |d| try argv.append(a, try std.fmt.allocPrint(a, "-D{s}", .{d}));
            try argv.appendSlice(a, &.{ src_path, "-o", ".zig-cache/params_probe.exe" });
            const cc = try std.process.run(a, io, .{ .argv = argv.items });
            const built = switch (cc.term) {
                .exited => |c| c == 0,
                else => false,
            };
            if (!built) {
                try out.print(a, "        .{s} => null,\n", .{name});
                std.debug.print("{s}/{s}: header does not build\n", .{ fam.name, name });
                continue;
            }
            const run = try std.process.run(a, io, .{ .argv = &.{".zig-cache/params_probe.exe"} });
            const vals = try parseValues(a, run.stdout);
            const have_p = vals.get("s.P") orelse null;
            const have_q = vals.get("s.Q") orelse null;
            if (have_p == null or have_q == null) {
                try out.print(a, "        .{s} => null,\n", .{name});
                std.debug.print("{s}/{s}: header defines no blocking\n", .{ fam.name, name });
                continue;
            }
            try out.print(a, "        .{s} => .{{\n", .{name});
            for (precs) |pr| {
                const key = struct {
                    fn k(al: std.mem.Allocator, f: []const u8, w: []const u8) ![]const u8 {
                        return std.fmt.allocPrint(al, "{s}.{s}", .{ f, w });
                    }
                }.k;
                const um = (vals.get(try key(a, pr.field, "UNROLL_M")) orelse null);
                const pp = (vals.get(try key(a, pr.field, "P")) orelse null);
                if (um == null and pp == null) continue;
                try out.print(a, "            .{s} = Shape{{ .unroll_m = {d}, .unroll_n = {d}, .unroll_mn = {d}, .align_k = {d}, .p_base = {d}, .p = ", .{
                    pr.field,
                    um orelse 0,
                    (vals.get(try key(a, pr.field, "UNROLL_N")) orelse null) orelse 0,
                    (vals.get(try key(a, pr.field, "UNROLL_MN")) orelse null) orelse 0,
                    (vals.get(try key(a, pr.field, "ALIGN_K")) orelse null) orelse 0,
                    (vals.get(try key(a, pr.field, "P_BASE")) orelse null) orelse 0,
                });
                try size(a, &out, pp, ".dynamic");
                try out.appendSlice(a, ", .q = ");
                try size(a, &out, (vals.get(try key(a, pr.field, "Q")) orelse null), ".dynamic");
                try out.appendSlice(a, ", .r = ");
                try size(a, &out, (vals.get(try key(a, pr.field, "R")) orelse null), ".from_buffer");
                try out.appendSlice(a, " },\n");
            }
            for (globals) |g| {
                if ((vals.get(g[0]) orelse null)) |x| {
                    const field = if (std.mem.eql(u8, g[0], "align")) "align_" else g[0];
                    try out.print(a, "            .{s} = {d},\n", .{ field, x });
                }
            }
            try out.appendSlice(a, "        },\n");
        }
        try out.appendSlice(a, "    };\n}\n\n");
    }

    const dest = "src/linalg/blas/blas3/gemm/params_data.zig";
    try Io.Dir.cwd().writeFile(io, .{ .sub_path = dest, .data = out.items });
    _ = try std.process.run(a, io, .{ .argv = &.{ "zig", "fmt", dest } });
    std.debug.print("wrote {s}\n", .{dest});
}
