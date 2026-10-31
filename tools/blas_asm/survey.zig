//! Lists, per architecture and core, which OpenBLAS source each kernel is
//! built from and what rstd serves in its place, so the gaps show.
//!
//!     OPENBLAS=<checkout> zig run tools/blas_asm/survey.zig [-- <routine> | l3 | l2 [<dest>]]
//!
//! With no argument, writes `blas1/<routine>/coverage.md` for every level-1
//! routine (the roles of `Makefile.L1`, grouped as `roles.dirFor` does), or
//! with `<routine>` for that one. With `l3`, the level-3 routines'
//! (`l3_families`), `gemm`'s among them. With `l2`, `blas2/coverage.md`
//! (`l2_families`), and which C sources carry assembly or intrinsics.
//! Run from the rstd root.
//!
//! Reads `kernel/<arch>/KERNEL.<core>`, falling back to `kernel/<arch>/KERNEL`
//! and then to the defaults in `kernel/Makefile.L*`, as the build does. A core
//! that assigns a kernel under a conditional lists every alternative, joined
//! by `|`; which one applies depends on the compiler or byte order. What rstd
//! serves is read from its generated indexes (`asm/`, `uk/`), the hand ports
//! of each `<arch>/kernel.zig` (`ports`), and its registries
//! (`kernels/<arch>.zig`, their `vector_rule`) and `vector.zig`.

const std = @import("std");
const Io = std.Io;
const asm_arch = @import("arch/__root__.zig");
const map = @import("kernel_map.zig");
const roles = @import("roles.zig");
const uk = @import("uk.zig");

/// An OpenBLAS `kernel/<name>` directory and where rstd serves it: the
/// routine subdirectory `dir`, the registry `kernels/<registry>.zig`, the
/// core enum `cpu/_coretypes.zig`'s `cores`. An empty `dir`: no rstd target.
const Arch = struct {
    name: []const u8,
    dir: []const u8 = "",
    registry: []const u8 = "",
    cores: []const u8 = "",
};

const arches = [_]Arch{
    .{ .name = "x86_64", .dir = "x86_64", .registry = "x86_64", .cores = "X86CoreType" },
    .{ .name = "x86", .dir = "x86", .registry = "x86", .cores = "X86CoreType" },
    .{ .name = "arm64", .dir = "arm64", .registry = "arm64", .cores = "Arm64CoreType" },
    .{ .name = "arm", .dir = "arm", .registry = "arm", .cores = "ArmCoreType" },
    .{ .name = "power", .dir = "power", .registry = "power", .cores = "PowerCoreType" },
    .{ .name = "mips64", .dir = "mips64", .registry = "mips", .cores = "MipsCoreType" },
    .{ .name = "mips", .dir = "mips64", .registry = "mips", .cores = "MipsCoreType" },
    .{ .name = "riscv64", .dir = "riscv64", .registry = "riscv", .cores = "RiscvCoreType" },
    .{ .name = "loongarch64", .dir = "loongarch64", .registry = "loongarch", .cores = "LoongarchCoreType" },
    .{ .name = "zarch", .dir = "s390x", .registry = "s390x", .cores = "S390xCoreType" },
    .{ .name = "sparc", .dir = "sparc", .registry = "sparc", .cores = "SparcCoreType" },
    .{ .name = "alpha" },
    .{ .name = "ia64" },
    .{ .name = "wasm", .dir = "wasm", .registry = "wasm", .cores = "WasmCoreType" },
};

const blas = "src/linalg/blas";

fn join(parts: []const []const u8) ![]const u8 {
    var out: std.ArrayList(u8) = .empty;
    for (parts, 0..) |p, i| {
        if (i > 0) try out.appendSlice(map.arena, " \\| ");
        try out.append(map.arena, '`');
        try out.appendSlice(map.arena, p);
        try out.append(map.arena, '`');
    }
    return out.items;
}

fn isAsm(name: []const u8) bool {
    return std.mem.endsWith(u8, name, ".S");
}

fn contains(list: []const []const u8, name: []const u8) bool {
    for (list) |l| {
        if (std.ascii.eqlIgnoreCase(l, name)) return true;
    }
    return false;
}

pub fn main(init: std.process.Init) !void {
    const arena = init.arena.allocator();
    const checkout = init.environ_map.get("OPENBLAS") orelse
        map.die("set OPENBLAS to the OpenBLAS checkout", .{});
    try map.init(arena, init.io, checkout);

    var args = try init.minimal.args.iterateAllocator(arena);
    _ = args.next();
    const mode = args.next();
    if (mode) |m| {
        if (std.mem.eql(u8, m, "l2")) return try l2(init, arena, args.next() orelse "src/linalg/blas/blas2/coverage.md");
        if (std.mem.eql(u8, m, "l3")) return try l3(init, arena);
    }
    try l1(init, arena, mode);
}

// ---- what upstream builds ---------------------------------------------------

/// What a core's source is: assembly, assembly the toolchain cannot
/// assemble, a C file whose block for that core (or the file itself, or what
/// it includes) has inline assembly or intrinsics, or plain C. A text
/// search, so a heuristic.
const Upstream = enum {
    assembly,
    unsupported,
    inline_asm,
    intrinsics,
    c,
    missing,

    fn tuned(u: Upstream) bool {
        return u == .assembly or u == .inline_asm or u == .intrinsics;
    }

    fn mark(u: Upstream) []const u8 {
        return switch (u) {
            .assembly, .c => "",
            .unsupported => " (unsupported)",
            .inline_asm => " +asm",
            .intrinsics => " +intrinsics",
            .missing => " (missing)",
        };
    }
};

/// A source's kind, and the CPU feature its intrinsics need (`needs`).
const Kind = struct { kind: Upstream, need: ?[]const u8 = null };

/// The feature each intrinsics header needs, as `cpu/_corespec.zig` names
/// it. Only one `machine` lists on every core that has it, so its absence
/// says the core lacks it: MSA, which the mips `KERNEL` files pick by the
/// compiler (`ifndef NO_MSA`), not the core.
const needs = [_][2][]const u8{.{ "msa.h", "msa" }};

fn textKind(t: []const u8) Kind {
    const need = for (needs) |n| {
        if (std.mem.find(u8, t, n[0]) != null) break n[1];
    } else null;
    return .{ .kind = if (hasAsm(t)) .inline_asm else if (hasIntrinsics(t)) .intrinsics else .c, .need = need };
}

/// `source`'s kind on `core`: a driver with a `#if defined(<core>)` chain
/// (`uk.readChain`) by the block it includes for `core`, or by its own text
/// when the chain has none for it; any other C file by itself and the `.c`
/// files it includes.
fn upstreamKind(arena: std.mem.Allocator, arch: []const u8, core: []const u8, source: []const u8) !Kind {
    const p = try map.path(&.{ map.root, "kernel", arch, source });
    const text = map.read(p) orelse return .{ .kind = .missing };
    if (isAsm(source)) return .{ .kind = if (assembles(arch, source)) .assembly else .unsupported };
    const chain = try uk.readChain(arena, text);
    if (chain.len > 0) {
        const dir = std.Io.Dir.path.dirname(p) orelse ".";
        for (chain) |link| {
            if (link.cores.len > 0 and !contains(link.cores, core)) continue;
            return textKind(map.read(try map.path(&.{ dir, link.file })) orelse return .{ .kind = .missing });
        }
        return textKind(text);
    }
    const info = try classifyC(arena, arch, source);
    var kind = textKind(text);
    for (info.includes.items) |i| {
        if (std.mem.endsWith(u8, i, "(asm)")) kind.kind = .inline_asm;
        if (kind.kind == .c and std.mem.endsWith(u8, i, "(intrinsics)")) kind.kind = .intrinsics;
    }
    return kind;
}

/// Whether the generator's `arch/<arch>.zig` can assemble `source`.
fn assembles(arch: []const u8, source: []const u8) bool {
    for (asm_arch.all) |a| {
        if (std.mem.eql(u8, a.name, arch)) return a.assembles(source);
    }
    return true;
}

/// The cell text of `sources`, each with its kind on `core`, and whether
/// any `core` can run is tuned.
fn upstreamCell(arena: std.mem.Allocator, ar: Arch, core: []const u8, sources: []const []const u8) !struct { text: []const u8, tuned: bool } {
    var out: std.ArrayList(u8) = .empty;
    var tuned = false;
    for (sources, 0..) |s, i| {
        const k = try upstreamKind(arena, ar.name, core, s);
        const lacks = if (k.need) |n| try lacksFeature(arena, ar, core, n) else false;
        tuned = tuned or (k.kind.tuned() and !lacks);
        if (i > 0) try out.appendSlice(arena, " \\| ");
        try out.print(arena, "`{s}`{s}", .{ s, k.kind.mark() });
        if (lacks) try out.print(arena, " (no {s})", .{k.need.?});
    }
    return .{ .text = out.items, .tuned = tuned };
}

/// The features `cpu/_corespec.zig`'s `machine` lists for each core, keyed
/// `<enum>.<core>`; read once.
var machine_features: ?std.array_hash_map.String([]const []const u8) = null;

/// Whether `machine` has an entry for `core` of `ar` without `feature`.
fn lacksFeature(arena: std.mem.Allocator, ar: Arch, core: []const u8, feature: []const u8) !bool {
    if (ar.cores.len == 0) return false;
    if (machine_features == null) machine_features = try readMachine(arena);
    const key = try arena.print("{s}.{s}", .{ ar.cores, try std.ascii.allocLowerString(arena, core) });
    const listed = machine_features.?.get(key) orelse return false;
    return !contains(listed, feature);
}

/// Reads `machine`'s prongs, one line each.
fn readMachine(arena: std.mem.Allocator) !std.array_hash_map.String([]const []const u8) {
    var out: std.array_hash_map.String([]const []const u8) = .empty;
    const text = map.read(blas ++ "/cpu/_corespec.zig") orelse map.die("cannot read cpu/_corespec.zig", .{});
    const at = std.mem.find(u8, text, "pub fn machine(") orelse map.die("no machine in cpu/_corespec.zig", .{});
    var enum_name: []const u8 = "";
    var lines = std.mem.splitScalar(u8, text[at..], '\n');
    while (lines.next()) |raw| {
        const line = std.mem.trim(u8, raw, " \t\r");
        if (std.mem.eql(u8, line, "}")) break;
        if (std.mem.startsWith(u8, line, "coretypes.")) {
            enum_name = line["coretypes.".len .. std.mem.findScalar(u8, line, ' ') orelse line.len];
            continue;
        }
        if (line.len == 0 or line[0] != '.') continue;
        const arrow = std.mem.find(u8, line, " => ") orelse continue;
        const value = line[arrow..];
        const features = if (field(value, "features")) |f| try quoted(arena, f[0 .. std.mem.findScalar(u8, f, '}') orelse f.len]) else &.{};
        var names = std.mem.tokenizeAny(u8, line[0..arrow], ", .");
        while (names.next()) |n| try out.put(arena, try arena.print("{s}.{s}", .{ enum_name, n }), features);
    }
    return out;
}

// ---- what rstd serves -------------------------------------------------------

/// A row of a generated index: a slot, the cores it names, and its block
/// kind (`uk/`) or nothing (`asm/`).
const Row = struct { slot: []const u8, cores: []const []const u8, block: []const u8 = "" };

/// Where rstd's per-core kernels for one routine on one architecture come
/// from, read from the generated indexes.
const Served = struct {
    /// The converted `.S` kernels (`asm/__root__.zig`'s `entries`).
    converted: []const Row = &.{},
    /// Upstream assembly left unconverted (its `skipped`).
    skipped: []const Row = &.{},
    /// The C drivers' blocks (`uk/__root__.zig`).
    blocks: []const Row = &.{},
    /// Ported files the routine's `kernel.zig` composes beside the indexes.
    hand: []const []const u8 = &.{},
    /// Its `ports` table, which says where those serve; null for none.
    ports: ?[]const Port = null,

    /// The port filling `slot` on `core` in `arch`'s build.
    fn port(s: Served, arch: []const u8, slot: []const u8, core: []const u8) ?Port {
        for (s.ports orelse return null) |p| {
            if (p.arch) |a| if (!std.mem.eql(u8, a, arch)) continue;
            if (contains(p.cores, core) and contains(p.slots, slot)) return p;
        }
        return null;
    }

    fn find(rows: []const Row, slot: []const u8, core: []const u8) ?Row {
        for (rows) |r| {
            if (std.mem.eql(u8, r.slot, slot) and contains(r.cores, core)) return r;
        }
        return null;
    }

    /// The core's block of the slot; none where two entries name it
    /// (`uk_block.Kernels`'s `found`).
    fn block(s: Served, slot: []const u8, core: []const u8) ?Row {
        var hit: ?Row = null;
        for (s.blocks) |r| {
            if (!std.mem.eql(u8, r.slot, slot) or !contains(r.cores, core)) continue;
            if (hit != null) return null;
            hit = r;
        }
        return hit;
    }
};

/// A row of a `kernel.zig`'s `ports` (`asm_registry.Port`).
const Port = struct {
    cores: []const []const u8,
    slots: []const []const u8,
    written: enum { assembly, vector, scalar },
    arch: ?[]const u8,
};

/// The names of the enum literals of `t`, a list starting `&.{`.
fn literals(arena: std.mem.Allocator, t: []const u8) ![]const []const u8 {
    if (!std.mem.startsWith(u8, t, "&.{")) map.die("not a literal list: {s}", .{t[0..@min(t.len, 40)]});
    const body = t[3 .. std.mem.findScalar(u8, t, '}') orelse t.len];
    var out: std.ArrayList([]const u8) = .empty;
    var it = std.mem.tokenizeAny(u8, body, ", \t\r\n");
    while (it.next()) |tok| {
        if (tok[0] != '.') map.die("not an enum literal: {s}", .{tok});
        try out.append(arena, tok[1..]);
    }
    return out.items;
}

/// What follows `.<name> = ` in `t`.
fn field(t: []const u8, comptime name: []const u8) ?[]const u8 {
    const tag = "." ++ name ++ " = ";
    const at = std.mem.find(u8, t, tag) orelse return null;
    return t[at + tag.len ..];
}

/// The entries of the `pub const ports` table of `text`, a `kernel.zig`.
fn readPorts(arena: std.mem.Allocator, text: []const u8) !?[]const Port {
    const tag = "pub const ports = [_]registry.Port(Core){";
    const at = std.mem.find(u8, text, tag) orelse return null;
    var out: std.ArrayList(Port) = .empty;
    var depth: usize = 0;
    var start: usize = 0;
    for (text[at + tag.len ..], at + tag.len..) |c, i| {
        if (c == '{') {
            if (depth == 0) start = i;
            depth += 1;
        } else if (c == '}') {
            if (depth == 0) break;
            depth -= 1;
            if (depth > 0) continue;
            const e = text[start..i];
            const written = field(e, "written") orelse map.die("a port without .written", .{});
            const name = written[1 .. std.mem.findAny(u8, written, ", }") orelse written.len];
            try out.append(arena, .{
                .cores = try literals(arena, field(e, "cores") orelse map.die("a port without .cores", .{})),
                .slots = try literals(arena, field(e, "slots") orelse map.die("a port without .slots", .{})),
                .written = std.meta.stringToEnum(@FieldType(Port, "written"), name) orelse map.die("unknown .written = .{s}", .{name}),
                .arch = if (field(e, "arch")) |a| (try quoted(arena, a))[0] else null,
            });
        }
    }
    return out.items;
}

/// How a registry routes the `@Vector` kernels: its `vector_rule`
/// (`vector.Rule`), or none.
const VectorRule = union(enum) {
    none,
    routed: ?[]const []const u8,
    unfilled,
    unfilled_on: []const []const u8,
};

fn readRule(arena: std.mem.Allocator, ar: Arch) !VectorRule {
    if (ar.registry.len == 0) return .none;
    const reg = map.read(try arena.print(blas ++ "/kernels/{s}.zig", .{ar.registry})) orelse return .none;
    const tag = "vector_rule: vector.Rule(Core) = ";
    const at = std.mem.find(u8, reg, tag) orelse {
        if (std.mem.find(u8, reg, "vector.Registry(") != null) map.die("kernels/{s}.zig has the @Vector kernels but no vector_rule", .{ar.registry});
        return .none;
    };
    const v = reg[at + tag.len .. std.mem.findScalarPos(u8, reg, at, ';') orelse reg.len];
    if (std.mem.eql(u8, v, ".unfilled")) return .unfilled;
    if (std.mem.eql(u8, v, ".{ .routed = null }")) return .{ .routed = null };
    if (field(v, "routed")) |r| return .{ .routed = try literals(arena, r) };
    if (field(v, "unfilled_on")) |r| return .{ .unfilled_on = try literals(arena, r) };
    map.die("kernels/{s}.zig: cannot read vector_rule = {s}", .{ ar.registry, v });
}

/// `vector.displacing`: the cores of enum `cores_enum` where the `@Vector`
/// kernels take every slot of theirs but `keep`; or `vector.flushing`: the
/// cores where they take no `s` or `c` slot (no `keep`).
const Displacing = struct { cores_enum: []const u8, cores: []const []const u8, keep: []const []const u8 = &.{} };

/// `vector.zig`'s `pub const <name>: struct { cores: ..., keep: ... }`.
fn readCores(arena: std.mem.Allocator, comptime name: []const u8, comptime keep: bool) !Displacing {
    const text = map.read(blas ++ "/vector.zig") orelse map.die("cannot read vector.zig", .{});
    const tag = "pub const " ++ name ++ ": struct { cores: []const coretypes.";
    const at = std.mem.find(u8, text, tag) orelse map.die("no " ++ name ++ " in vector.zig", .{});
    const rest = text[at + tag.len ..];
    const body = rest[0 .. std.mem.find(u8, rest, "};") orelse rest.len];
    const value = body[std.mem.find(u8, body, "} = .{") orelse map.die("cannot read vector." ++ name, .{}) ..];
    return .{
        .cores_enum = rest[0 .. std.mem.findAny(u8, rest, ", ") orelse 0],
        .cores = try literals(arena, field(value, "cores") orelse map.die("no " ++ name ++ ".cores", .{})),
        .keep = if (keep) try literals(arena, field(value, "keep") orelse map.die("no " ++ name ++ ".keep", .{})) else &.{},
    };
}

/// The quoted strings of `t`.
fn quoted(arena: std.mem.Allocator, t: []const u8) ![]const []const u8 {
    var out: std.ArrayList([]const u8) = .empty;
    var i: usize = 0;
    while (std.mem.findScalarPos(u8, t, i, '"')) |q| {
        const e = std.mem.findScalarPos(u8, t, q + 1, '"') orelse break;
        try out.append(arena, t[q + 1 .. e]);
        i = e + 1;
    }
    return out.items;
}

/// The rows of a generated index: each `.slot = .<slot>` up to the next.
/// An `asm/` entry's cores are the `keys` of the file it imports.
fn readIndex(arena: std.mem.Allocator, dir: []const u8, text: []const u8, entries: *std.ArrayList(Row), cored: *std.ArrayList(Row)) !void {
    const tag = ".slot = .";
    var at = std.mem.find(u8, text, tag);
    while (at) |a| {
        const start = a + tag.len;
        const next = std.mem.findPos(u8, text, start, tag);
        const body = text[start .. next orelse text.len];
        at = next;
        const slot = body[0 .. std.mem.findScalar(u8, body, ',') orelse continue];
        if (std.mem.find(u8, body, ".cores = &.{")) |c| {
            const list = body[c..][0 .. std.mem.findScalar(u8, body[c..], '}') orelse continue];
            var row: Row = .{ .slot = slot, .cores = try quoted(arena, list) };
            if (std.mem.find(u8, body, ".block = .")) |b| {
                const rest = body[b + ".block = .".len ..];
                row.block = rest[0 .. std.mem.findScalar(u8, rest, ',') orelse rest.len];
            }
            try cored.append(arena, row);
        } else if (std.mem.find(u8, body, "@import(\"")) |m| {
            const file = (try quoted(arena, body[m..]))[0];
            const mod = map.read(try map.path(&.{ dir, file })) orelse continue;
            const keys_at = std.mem.find(u8, mod, "const keys") orelse continue;
            const keys = mod[keys_at .. std.mem.findPos(u8, mod, keys_at, "};") orelse mod.len];
            var cores: std.ArrayList([]const u8) = .empty;
            var lines = std.mem.splitScalar(u8, keys, '\n');
            while (lines.next()) |l| {
                const q = try quoted(arena, l);
                if (q.len > 0) try cores.append(arena, q[0]);
            }
            try entries.append(arena, .{ .slot = slot, .cores = cores.items });
        }
    }
}

fn readServed(arena: std.mem.Allocator, ar: Arch, level: []const u8, routine: []const u8) !Served {
    var s: Served = .{};
    if (ar.dir.len == 0) return s;
    const reg = map.read(try arena.print(blas ++ "/kernels/{s}.zig", .{ar.registry})) orelse return s;
    const here = try arena.print("{s}/{s}/{s}", .{ level, routine, ar.dir });
    // an index the registry does not list serves nothing
    if (std.mem.find(u8, reg, try arena.print("../{s}/", .{here})) == null) return s;
    const root = try arena.print(blas ++ "/{s}", .{here});
    var converted: std.ArrayList(Row) = .empty;
    var skipped: std.ArrayList(Row) = .empty;
    var blocks: std.ArrayList(Row) = .empty;
    const asm_dir = try map.path(&.{ root, "asm" });
    if (map.read(try map.path(&.{ asm_dir, "__root__.zig" }))) |t| try readIndex(arena, asm_dir, t, &converted, &skipped);
    const uk_dir = try map.path(&.{ root, "uk" });
    var unused: std.ArrayList(Row) = .empty;
    if (map.read(try map.path(&.{ uk_dir, "__root__.zig" }))) |t| try readIndex(arena, uk_dir, t, &unused, &blocks);
    const kernel = map.read(try map.path(&.{ root, "kernel.zig" }));
    // the split indexes it imports (gemv's `uk/index_n.zig`, `uk/index_t.zig`)
    if (kernel) |t| for (try quoted(arena, t)) |file| {
        if (!std.mem.startsWith(u8, file, "uk/index")) continue;
        if (map.read(try map.path(&.{ root, file }))) |x| try readIndex(arena, uk_dir, x, &unused, &blocks);
    };
    s.converted = converted.items;
    s.skipped = skipped.items;
    s.blocks = blocks.items;
    // a sibling it imports that is neither an index nor the shared driver
    if (kernel) |t| {
        s.ports = try readPorts(arena, t);
        var hand: std.ArrayList([]const u8) = .empty;
        var i: usize = 0;
        while (std.mem.findPos(u8, t, i, "@import(\"")) |m| {
            i = m + 1;
            const file = (try quoted(arena, t[m..]))[0];
            if (std.mem.startsWith(u8, file, "../../") or std.mem.startsWith(u8, file, "asm/") or
                std.mem.startsWith(u8, file, "uk/") or std.mem.eql(u8, file, "../kernel.zig") or
                !std.mem.endsWith(u8, file, ".zig")) continue;
            try hand.append(arena, file);
        }
        s.hand = hand.items;
    }
    return s;
}

/// The names of `cpu/_coretypes.zig`'s `enum_name`.
fn coreNames(arena: std.mem.Allocator, enum_name: []const u8) ![]const []const u8 {
    var out: std.ArrayList([]const u8) = .empty;
    if (enum_name.len == 0) return out.items;
    const text = map.read(blas ++ "/cpu/_coretypes.zig") orelse map.die("cannot read cpu/_coretypes.zig; run from the rstd root", .{});
    const at = std.mem.find(u8, text, try arena.print("pub const {s} = enum {{", .{enum_name})) orelse map.die("no {s}", .{enum_name});
    const open = std.mem.findScalarPos(u8, text, at, '{').? + 1;
    const body = text[open..std.mem.findScalarPos(u8, text, open, '}').?];
    var lines = std.mem.splitScalar(u8, body, '\n');
    while (lines.next()) |raw| {
        const line = raw[0 .. std.mem.find(u8, raw, "//") orelse raw.len];
        var names = std.mem.tokenizeAny(u8, line, ", \t\r");
        while (names.next()) |n| try out.append(arena, n);
    }
    return out.items;
}

/// The slots `vector.zig` has kernels for: its `real_slots` and `complex_slots`.
fn vectorSlots(arena: std.mem.Allocator) ![]const []const u8 {
    const text = map.read(blas ++ "/vector.zig") orelse map.die("cannot read vector.zig", .{});
    const start = std.mem.find(u8, text, "const real_slots") orelse map.die("no real_slots in vector.zig", .{});
    const end = std.mem.findPos(u8, text, start, "/// A kernel:") orelse text.len;
    var out: std.ArrayList([]const u8) = .empty;
    var lines = std.mem.splitScalar(u8, text[start..end], '\n');
    while (lines.next()) |raw| {
        const line = std.mem.trim(u8, raw, " \t\r");
        const eq = std.mem.find(u8, line, " = .") orelse continue;
        if (!std.mem.endsWith(u8, line, ",")) continue;
        try out.append(arena, line[eq + " = .".len .. line.len - 1]);
    }
    return out.items;
}

/// What rstd serves for one slot on one core.
const Rstd = enum {
    converted,
    port_asm,
    port_vec,
    port,
    vector,
    portable,
    routed,
    extended,
    na,

    fn label(r: Rstd) []const u8 {
        return switch (r) {
            .converted => "asm",
            .port_asm => "C+asm",
            .port_vec => "C+vec",
            .port => "C",
            .vector => "vec",
            .portable, .extended => "portable",
            .routed => "routed",
            .na => "n/a",
        };
    }
};

/// A slot's `Rstd`, and whether a hand port may take it instead.
const Got = struct { r: Rstd, unsure: bool = false };

/// How the registry routes the `@Vector` kernels (`vector.takes`).
const Routing = struct { rule: VectorRule, displacing: Displacing, flushing: Displacing };

/// What `core` gets for `slot`, as the registry's `install` picks it: the
/// `@Vector` kernels where its `vector_rule` takes the slot, else a hand
/// port, else the converted assembly, else the ported driver, else the
/// portable kernels.
fn serve(ar: Arch, s: Served, routing: Routing, core: []const u8, role: roles.Role, fills: []const []const u8) Got {
    const slot = role.slot.?;
    if (role.extended) return .{
        .r = .extended,
    };
    const port = s.port(ar.name, slot, core);
    const indexed = Served.find(s.converted, slot, core) != null or s.block(slot, core) != null or port != null;
    const skip = roles.skipped(ar.name, slot);
    const routed = if (skip) |k| !k.cores or Served.find(s.skipped, slot, core) != null else false;
    const d = routing.displacing;
    const displaced = std.mem.eql(u8, ar.cores, d.cores_enum) and contains(d.cores, core) and !contains(d.keep, slot);
    const f = routing.flushing;
    const flushed = std.mem.eql(u8, ar.cores, f.cores_enum) and contains(f.cores, core) and (role.prec == .s or role.prec == .c);
    const vec = contains(fills, slot) and !flushed and switch (routing.rule) {
        .none => false,
        .routed => |cores| routed and (cores == null or contains(cores.?, core)),
        .unfilled => routed or !indexed or displaced,
        .unfilled_on => |cores| contains(cores, core) and (routed or !indexed or displaced),
    };
    const unsure = s.hand.len > 0 and s.ports == null and !indexed;
    if (vec) return .{ .r = .vector, .unsure = unsure };
    if (port) |p| return .{ .r = switch (p.written) {
        .assembly => .port_asm,
        .vector => .port_vec,
        .scalar => .port,
    } };
    if (Served.find(s.converted, slot, core) != null) return .{
        .r = .converted,
    };
    if (s.block(slot, core)) |b| {
        const r: Rstd = if (std.mem.eql(u8, b.block, "assembly")) .port_asm else if (std.mem.eql(u8, b.block, "vector")) .port_vec else .port;
        return .{
            .r = r,
        };
    }
    return .{ .r = if (routed) .routed else .portable, .unsure = unsure };
}

/// How a page reads, after its first paragraph.
const legend =
    \\Upstream: a `.S` file is assembly, `(unsupported)` where LLVM cannot
    \\assemble it (the generator's `arch/<arch>.zig` leaves it out); `+asm` and
    \\`+intrinsics` mark C whose block for that core (or the file, or a `.c` it
    \\includes) has inline assembly or intrinsics, a text search; unmarked C is
    \\plain C. `(no <feature>)` marks intrinsics the core cannot run
    \\(`cpu/_corespec.zig`'s `machine`), which upstream picks by the compiler.
    \\Neither counts as tuned.
    \\
    \\rstd:
    \\
    \\- `asm`: the converted `.S` (`asm/`);
    \\- `C+asm`, `C+vec`, `C`: the Zig port of the C upstream: around its
    \\  converted inline-assembly block (`uk/`) or in `asm` written for its
    \\  intrinsics, on `@Vector`s for its intrinsics, or in plain Zig (a hand
    \\  port is in its `kernel.zig`'s `ports`);
    \\- `vec`: the `@Vector` kernels (`vector.zig`);
    \\- `routed`: the portable kernels, by design (`vector.routed`);
    \\- `portable`: the portable kernels, **bold** where upstream is tuned (a gap);
    \\- `n/a`: rstd has no such core.
    \\
    \\The `Q` and `X` columns stay portable by design. On Apple's systems arm64
    \\takes `vec` only for routed slots. A core whose kernel object is not linked
    \\runs those of the core that serves it (`asm_registry.served`).
    \\
;

/// What a page's rstd cells add up to: the (core, role) pairs by label, and
/// the gaps, one row per core.
const Tally = struct {
    counts: std.array_hash_map.String(usize) = .empty,
    gaps: std.ArrayList(u8) = .empty,
    n_gaps: usize = 0,
    /// Gaps per architecture.
    by_arch: std.array_hash_map.String(usize) = .empty,

    /// The rstd half of a cell: what `core` gets for each role of `group`,
    /// **bold** where one is a gap; those roles go on `gap_roles`.
    fn cell(t: *Tally, arena: std.mem.Allocator, ar: Arch, s: Served, routing: Routing, core: []const u8, here: bool, group: []const roles.Role, tuned: bool, fills: []const []const u8, gap_roles: *std.ArrayList([]const u8)) ![]const u8 {
        var labels: std.ArrayList([]const u8) = .empty;
        var gap = false;
        for (group) |role| {
            const got = if (here) serve(ar, s, routing, core, role, fills) else Got{ .r = .na };
            const is_gap = got.r == .portable and tuned and !got.unsure;
            gap = gap or is_gap;
            const tag = try arena.print("{s}{s}", .{ got.r.label(), if (got.unsure) "?" else "" });
            const e = try t.counts.getOrPut(arena, tag);
            if (!e.found_existing) e.value_ptr.* = 0;
            e.value_ptr.* += 1;
            if (is_gap) {
                t.n_gaps += 1;
                const a = try t.by_arch.getOrPut(arena, ar.name);
                if (!a.found_existing) a.value_ptr.* = 0;
                a.value_ptr.* += 1;
                try gap_roles.append(arena, role.name);
            }
            if (!contains(labels.items, tag)) try labels.append(arena, tag);
        }
        const rstd = try std.mem.join(arena, "/", labels.items);
        return if (gap) try arena.print("**{s}**", .{rstd}) else rstd;
    }

    fn row(t: *Tally, arena: std.mem.Allocator, arch: []const u8, core: []const u8, gap_roles: []const []const u8) !void {
        if (gap_roles.len > 0) try t.gaps.print(arena, "| {s} | {s} | {s} |\n", .{ arch, core, try std.mem.join(arena, ", ", gap_roles) });
    }

    /// The `Gaps` and `Totals` sections.
    fn write(t: Tally, arena: std.mem.Allocator, w: *std.ArrayList(u8)) !void {
        try w.appendSlice(arena, "## Gaps\n\n");
        if (t.gaps.items.len == 0) {
            try w.appendSlice(arena, "None: no core runs the portable kernels where upstream is tuned.\n\n");
        } else {
            try w.appendSlice(arena, "The roles a core runs portable where upstream is tuned. By architecture:");
            for (t.by_arch.keys(), t.by_arch.values(), 0..) |k, v, i| try w.print(arena, "{s} {s} {d}", .{ if (i > 0) "," else "", k, v });
            try w.appendSlice(arena, ".\n\n| architecture | core | roles |\n|---|---|---|\n");
            try w.appendSlice(arena, t.gaps.items);
            try w.append(arena, '\n');
        }
        try w.appendSlice(arena, "## Totals\n\n(core, role) pairs by what rstd serves:");
        for (t.counts.keys(), t.counts.values(), 0..) |k, v, i| try w.print(arena, "{s} `{s}` {d}", .{ if (i > 0) "," else "", k, v });
        try w.print(arena, ". {d} gaps.\n", .{t.n_gaps});
    }
};

/// The line under an architecture's heading: no rstd target, or hand ports
/// this survey cannot place.
fn note(arena: std.mem.Allocator, w: *std.ArrayList(u8), ar: Arch, s: Served) !void {
    if (ar.dir.len == 0) {
        try w.appendSlice(arena, "rstd has no target here.\n\n");
    } else if (s.hand.len > 0 and s.ports == null) {
        var files: std.ArrayList(u8) = .empty;
        for (s.hand, 0..) |f, i| try files.print(arena, "{s}`{s}`", .{ if (i > 0) ", " else "", f });
        try w.print(arena, "`{s}/kernel.zig` also composes {s}, which this survey cannot read; `?` marks a slot one of them may take.\n\n", .{ ar.dir, files.items });
    }
}

// ---- level 1 ----------------------------------------------------------------

fn l1(init: std.process.Init, arena: std.mem.Allocator, only: ?[]const u8) !void {
    const text = map.read(try map.path(&.{ map.root, "kernel", "Makefile.L1" })).?;
    var routines: std.array_hash_map.String(std.ArrayList(roles.Role)) = .empty;
    for (try roles.parse(arena, text)) |r| {
        if (r.slot == null) continue;
        const dir = try roles.dirFor(arena, r.key);
        if (only) |o| if (!std.mem.eql(u8, o, dir)) continue;
        if (map.read(try arena.print(blas ++ "/blas1/{s}/__root__.zig", .{dir})) == null) continue;
        const e = try routines.getOrPut(arena, dir);
        if (!e.found_existing) e.value_ptr.* = .empty;
        try e.value_ptr.append(arena, r);
    }
    if (routines.count() == 0) map.die("no level-1 routine '{s}'", .{only orelse ""});
    const fills = try vectorSlots(arena);
    const displacing = try readCores(arena, "displacing", true);
    const flushing = try readCores(arena, "flushing", false);
    for (routines.keys(), routines.values()) |name, rs| try routine1(init, arena, name, rs.items, fills, displacing, flushing);
}

fn routine1(init: std.process.Init, arena: std.mem.Allocator, routine: []const u8, rs: []const roles.Role, fills: []const []const u8, displacing: Displacing, flushing: Displacing) !void {
    // one column per make variable, holding the roles built from it
    var columns: std.array_hash_map.String(std.ArrayList(roles.Role)) = .empty;
    const ordered = try arena.dupe(roles.Role, rs);
    std.mem.sort(roles.Role, ordered, {}, struct {
        // by precision, then the half-precision ones, the extended pair last
        fn rank(r: roles.Role) usize {
            const half = std.mem.startsWith(u8, r.key, "B") or std.mem.startsWith(u8, r.key, "SB");
            return @as(usize, if (r.extended) 2 else if (half) 1 else 0) * 4 + @backingInt(r.prec);
        }
        fn lt(_: void, a: roles.Role, b: roles.Role) bool {
            return rank(a) < rank(b);
        }
    }.lt);
    for (ordered) |r| {
        const e = try columns.getOrPut(arena, r.key);
        if (!e.found_existing) e.value_ptr.* = .empty;
        try e.value_ptr.append(arena, r);
    }

    var w: std.ArrayList(u8) = .empty;
    try w.print(arena,
        \\# `{s}` kernel coverage
        \\
        \\Generated by `tools/blas_asm/survey.zig` from OpenBLAS `kernel/` (`{s}`) and rstd's
        \\kernel indexes; re-derive it rather than editing it. Each cell is the
        \\source that core's `KERNEL.<core>` selects, falling back to the
        \\architecture's `KERNEL`, then to the `Makefile.L1` default (`|` separates
        \\alternatives chosen by a conditional), and after the arrow what rstd
        \\serves for that core.
        \\
        \\
    ++ legend ++ "\n", .{ routine, map.pin });

    var t: Tally = .{};

    for (arches) |ar| {
        const cores = try map.cores(ar.name) orelse continue;
        const known = try coreNames(arena, ar.cores);
        const s = try readServed(arena, ar, "blas1", routine);
        const routing: Routing = .{ .rule = try readRule(arena, ar), .displacing = displacing, .flushing = flushing };
        try w.print(arena, "## {s}\n\n", .{ar.name});
        try note(arena, &w, ar, s);
        try w.appendSlice(arena, "| core |");
        for (columns.values()) |group| {
            try w.append(arena, ' ');
            for (group.items, 0..) |r, i| try w.print(arena, "{s}{s}", .{ if (i > 0) "/" else "", r.name });
            try w.appendSlice(arena, " |");
        }
        try w.appendSlice(arena, "\n|---|");
        for (columns.keys()) |_| try w.appendSlice(arena, "---|");
        try w.append(arena, '\n');
        for (cores) |core| {
            const lower = try std.ascii.allocLowerString(arena, core);
            const here = ar.dir.len > 0 and contains(known, lower);
            var gap_roles: std.ArrayList([]const u8) = .empty;
            try w.print(arena, "| {s} |", .{core});
            for (columns.keys(), columns.values()) |key, group| {
                const up = try upstreamCell(arena, ar, core, try map.resolveKey(ar.name, core, key, "Makefile.L1"));
                const rstd = try t.cell(arena, ar, s, routing, lower, here, group.items, up.tuned, fills, &gap_roles);
                try w.print(arena, " {s} → {s} |", .{ up.text, rstd });
            }
            try w.append(arena, '\n');
            try t.row(arena, ar.name, core, gap_roles.items);
        }
        try w.append(arena, '\n');
    }

    try t.write(arena, &w);
    const dest = try arena.print(blas ++ "/blas1/{s}/coverage.md", .{routine});
    try Io.Dir.cwd().writeFile(init.io, .{ .sub_path = dest, .data = try alignTables(arena, w.items) });
    std.debug.print("wrote {s}: {d} gaps\n", .{ dest, t.n_gaps });
}

// ---- level 3 ----------------------------------------------------------------

/// A level-3 table column: a label and the make variable of its source.
const L3Col = struct {
    label: []const u8 = "",
    key: []const u8,
    /// Built from this variable where `key` is unset or `gate` is off.
    fallback: ?[]const u8 = null,
    /// The `Makefile.L3` flag that builds `key` (`USE_TRMM`); off, the
    /// column is `fallback`, or not built (`-`).
    gate: ?[]const u8 = null,
    /// What rstd has for it: the gemm tile's Zig port, or nothing yet.
    rstd: enum { tile, none } = .tile,
    /// Where set, the upstream sources rstd ports for it: any other tuned
    /// source is a gap, whatever `rstd` says.
    ported: ?[]const []const u8 = null,
};

const L3Table = struct { name: []const u8, cols: []const L3Col };
const L3Family = struct { dir: []const u8, tables: []const L3Table, note: []const u8 = "" };

/// `t` with `#` replaced by `p` and `@` by `f`.
fn fill(comptime t: []const u8, comptime p: []const u8, comptime f: []const u8) []const u8 {
    return comptime blk: {
        @setEvalBranchQuota(100_000);
        var out: []const u8 = "";
        for (t) |ch| {
            const piece: []const u8 = if (ch == '#') p else if (ch == '@') f else &[_]u8{ch};
            out = out ++ piece;
        }
        break :blk out;
    };
}

/// `base` for every precision `precs` and form `forms`, labeled `label`.
fn cols(comptime precs: []const []const u8, comptime forms: []const []const u8, comptime label: []const u8, comptime base: L3Col) []const L3Col {
    return comptime blk: {
        @setEvalBranchQuota(100_000);
        var out: []const L3Col = &.{};
        for (precs) |p| for (forms) |f| {
            var c = base;
            c.label = fill(label, p, f);
            c.key = fill(base.key, p, f);
            if (base.fallback) |fb| c.fallback = fill(fb, p, f);
            out = out ++ [_]L3Col{c};
        };
        const final = out[0..out.len].*;
        break :blk &final;
    };
}

const sdcz = &[_][]const u8{ "S", "D", "C", "Z" };
const cz = &[_][]const u8{ "C", "Z" };
const one = &[_][]const u8{""};
const gemm_kernels = cols(sdcz, one, "# gemm", .{ .key = "#GEMMKERNEL" });

fn direct(comptime routine: []const u8) []const L3Col {
    return cols(&.{"S"}, one, "S direct", .{ .key = "S" ++ routine ++ "DIRECTKERNEL_ALPHA_BETA", .gate = "USE_DIRECT_S" ++ routine, .rstd = .none });
}

/// The `direct` columns' caveat.
const sme = "arm64's `S direct` kernel (small matrices) runs only on a CPU with SME.";

/// `gemm`'s packing copies for precision `p`.
fn gemmCopies(comptime p: []const u8) []const L3Col {
    return &.{
        .{ .label = p ++ " incopy", .key = p ++ "GEMMINCOPY" },
        .{ .label = p ++ " itcopy", .key = p ++ "GEMMITCOPY" },
        .{ .label = p ++ " oncopy", .key = p ++ "GEMMONCOPY" },
        .{ .label = p ++ " otcopy", .key = p ++ "GEMMOTCOPY" },
    };
}

/// The small-matrix kernels (`SMALL_MATRIX_OPT`): the `B0` forms name the
/// same files on every core. rstd ports Skylake-X's (`gemm/x86_64/kernel.zig`)
/// and Cooper Lake's `sbgemm` ones; the others run `kernel/generic`'s, and
/// `sbgemm` nothing.
const small_forms = &[_][]const u8{ "NN", "NT", "TN", "TT" };
const small_kernels = cols(&.{ "S", "D" }, small_forms, "# small @", .{ .key = "#GEMM_SMALL_K_@", .ported = &.{
    "sgemm_small_kernel_nn_skylakex.c", "sgemm_small_kernel_nt_skylakex.c", "sgemm_small_kernel_tn_skylakex.c", "sgemm_small_kernel_tt_skylakex.c",
    "dgemm_small_kernel_nn_skylakex.c", "dgemm_small_kernel_nt_skylakex.c", "dgemm_small_kernel_tn_skylakex.c", "dgemm_small_kernel_tt_skylakex.c",
} }) ++ cols(&.{"SB"}, small_forms, "# small @", .{ .key = "#GEMM_SMALL_K_@", .rstd = .none, .ported = &.{
    "sbgemm_small_kernel_nn_cooperlake.c", "sbgemm_small_kernel_nt_cooperlake.c", "sbgemm_small_kernel_tn_cooperlake.c", "sbgemm_small_kernel_tt_cooperlake.c",
} });

const gemm3m = [_]L3Table{.{ .name = "kernels", .cols = cols(cz, one, "#", .{ .key = "#GEMM3MKERNEL", .gate = "USE_GEMM3M" }) }};

/// The 3M routines' note: their tile, and what upstream builds.
const three_m = "rstd's 3M kernel is the real gemm tile writing both parts of `C` (`gemm3m/kernel.zig`), and its copies `kernel/generic`'s 3M copies; " ++
    "upstream exports `gemm3m` on x86, x86_64, ia64 and MIPS only (`interface/Makefile`), and `symm3m` and `hemm3m` on none.";

/// The level-3 routines and the variables their kernels come from.
const l3_families = [_]L3Family{
    .{ .dir = "gemm", .tables = &.{
        .{ .name = "kernels", .cols = gemm_kernels },
        .{ .name = "beta", .cols = cols(sdcz, one, "# beta", .{ .key = "#GEMM_BETA" }) },
        .{ .name = "copies", .cols = gemmCopies("S") ++ gemmCopies("Z") },
        .{ .name = "small-matrix kernels", .cols = small_kernels },
    }, .note = "Only the `S` and `Z` copies are shown; `D` and `C` name the same layouts in nearly every source. " ++
        "The small-matrix kernels' `B0` forms name the same files; rstd ports Skylake-X's and Cooper Lake's (`gemm/x86_64/kernel.zig`), " ++
        "the latter in the portable half gemm's bits, and runs `kernel/generic`'s elsewhere (`gemm/small_kernel.zig`), where `sbgemm` has none." },
    .{ .dir = "symm", .note = sme, .tables = &.{
        .{ .name = "kernels", .cols = gemm_kernels ++ direct("SYMM") },
        .{ .name = "copies", .cols = cols(sdcz, &.{ "U", "L" }, "# @", .{ .key = "#SYMM@COPY_M" }) },
    } },
    .{ .dir = "hemm", .tables = &.{
        .{ .name = "kernels", .cols = cols(cz, one, "# gemm", .{ .key = "#GEMMKERNEL" }) },
        .{ .name = "copies", .cols = cols(cz, &.{ "UT", "LT" }, "# @", .{ .key = "#HEMM@COPY_M" }) },
    } },
    .{ .dir = "syrk", .tables = &.{.{ .name = "kernels", .cols = gemm_kernels ++ direct("SYRK") }}, .note = "`herk` runs the same kernels. " ++ sme },
    .{ .dir = "syr2k", .tables = &.{.{ .name = "kernels", .cols = gemm_kernels ++ direct("SYR2K") }}, .note = "`her2k` runs the same kernels. " ++ sme },
    .{ .dir = "gemmt", .tables = &.{.{ .name = "kernels", .cols = gemm_kernels }} },
    .{ .dir = "trmm", .tables = &.{
        .{ .name = "kernels", .cols = cols(sdcz, one, "#", .{ .key = "#TRMMKERNEL", .gate = "USE_TRMM", .fallback = "#GEMMKERNEL" }) ++
            cols(&.{"S"}, one, "S direct", .{ .key = "STRMMDIRECTKERNEL", .gate = "USE_DIRECT_STRMM", .rstd = .none }) },
        .{ .name = "copies", .cols = cols(sdcz, &.{ "UN", "LN", "UT", "LT" }, "# @", .{ .key = "#TRMM@COPY_M" }) },
    }, .note = "Without `USE_TRMM` upstream builds the `trmm` kernel from the gemm kernel's source (`-DTRMMKERNEL`). " ++ sme },
    .{ .dir = "trsm", .tables = &.{
        .{ .name = "kernels", .cols = cols(sdcz, &.{ "LN", "LT", "RN", "RT" }, "# @", .{ .key = "#TRSMKERNEL_@" }) },
        .{ .name = "copies", .cols = cols(&.{"real"}, &.{ "UN", "LN", "UT", "LT" }, "# @", .{ .key = "TRSMCOPY@_M" }) ++
            cols(&.{"complex"}, &.{ "UN", "LN", "UT", "LT" }, "# @", .{ .key = "ZTRSMCOPY@_M" }) },
    } },
    .{ .dir = "gemm3m", .tables = &gemm3m, .note = three_m },
    .{ .dir = "symm3m", .tables = &gemm3m, .note = "`symm3m` runs the `gemm3m` kernel. " ++ three_m },
    .{ .dir = "hemm3m", .tables = &gemm3m, .note = "`hemm3m` runs the `gemm3m` kernel. " ++ three_m },
};

/// Whether `Makefile.L3` sets `flag` to 1 for `core` of `arch`. Its `ifeq`
/// and `ifneq` tests on `ARCH`, `CORE` and `TARGET` are decided; any other
/// counts as taken, both sides.
fn flagSet(flag: []const u8, arch: []const u8, core: []const u8) bool {
    const text = map.read(map.path(&.{ map.root, "kernel", "Makefile.L3" }) catch return false) orelse return false;
    const Cond = struct { known: bool, taken: bool };
    var stack: [32]Cond = undefined;
    var depth: usize = 0;
    var value = false;
    var lines = std.mem.splitScalar(u8, text, '\n');
    while (lines.next()) |raw| {
        const line = std.mem.trim(u8, raw, " \t\r");
        const neg = std.mem.startsWith(u8, line, "ifneq");
        if (neg or std.mem.startsWith(u8, line, "ifeq")) {
            var c: Cond = .{ .known = false, .taken = true };
            const open = std.mem.find(u8, line, "$(");
            const close = if (open) |o| std.mem.findScalarPos(u8, line, o, ')') else null;
            const comma = if (close) |e| std.mem.findScalarPos(u8, line, e, ',') else null;
            if (comma) |cm| {
                const name = line[open.? + 2 .. close.?];
                const val = std.mem.trim(u8, line[cm + 1 ..], " )");
                const eq: ?bool = if (std.mem.eql(u8, name, "ARCH"))
                    std.mem.eql(u8, val, arch)
                else if (std.mem.eql(u8, name, "CORE") or std.mem.eql(u8, name, "TARGET") or std.mem.eql(u8, name, "TARGET_CORE"))
                    std.ascii.eqlIgnoreCase(val, core)
                else
                    null;
                if (eq) |e| c = .{ .known = true, .taken = e != neg };
            }
            stack[depth] = c;
            depth += 1;
        } else if (std.mem.startsWith(u8, line, "ifdef") or std.mem.startsWith(u8, line, "ifndef")) {
            stack[depth] = .{ .known = false, .taken = true };
            depth += 1;
        } else if (std.mem.startsWith(u8, line, "else")) {
            if (depth > 0 and stack[depth - 1].known) stack[depth - 1].taken = !stack[depth - 1].taken;
        } else if (std.mem.startsWith(u8, line, "endif")) {
            depth -|= 1;
        } else if (std.mem.startsWith(u8, line, flag)) {
            const rest = std.mem.trimStart(u8, line[flag.len..], " \t");
            if (rest.len == 0 or rest[0] != '=') continue;
            const live = for (stack[0..depth]) |c| {
                if (!c.taken) break false;
            } else true;
            if (live) value = std.mem.eql(u8, std.mem.trim(u8, rest[1..], " \t"), "1");
        }
    }
    return value;
}

/// The sources `core` builds column `c` from; null where it is not built.
fn l3Sources(ar: Arch, core: []const u8, c: L3Col) !?[]const []const u8 {
    var key: []const u8 = c.key;
    if (c.gate) |g| if (!flagSet(g, ar.name, core)) {
        key = c.fallback orelse return null;
    };
    const v = try map.resolveKey(ar.name, core, key, "Makefile.L3");
    if (v.len > 0) return v;
    if (c.fallback) |f| return try map.resolveKey(ar.name, core, f, "Makefile.L3");
    return v;
}

fn l3(init: std.process.Init, arena: std.mem.Allocator) !void {
    for (l3_families) |fam| try family3(init, arena, fam);
}

fn family3(init: std.process.Init, arena: std.mem.Allocator, fam: L3Family) !void {
    var w: std.ArrayList(u8) = .empty;
    try w.print(arena,
        \\# `{s}` kernel coverage
        \\
        \\Generated by `tools/blas_asm/survey.zig l3` from OpenBLAS `kernel/` (`{s}`);
        \\re-derive it rather than editing it. Each cell is the source that core's
        \\`KERNEL.<core>` selects for the column's make variable, falling back to the
        \\architecture's `KERNEL`, then to the `Makefile.L3` default; `generic` where
        \\that default is `kernel/generic`'s file for the core's unroll, `-` where
        \\upstream does not build it. `|` separates alternatives chosen by a
        \\conditional. `+asm` and `+intrinsics` mark C with inline assembly or
        \\intrinsics in the file or a `.c` it includes (a text search);
        \\`(unsupported)` marks assembly LLVM cannot assemble, `(no <feature>)`
        \\intrinsics the core cannot run, neither counted as tuned.
        \\
        \\The `rstd` row says what rstd serves for every core of the architecture:
        \\`tile` is the Zig port of `kernel/generic` over the `@Vector` gemm tile
        \\(`gemm/kernel.zig`), built into each core's kernel object (a core with none
        \\runs the portable build of the same source); `none` is nothing yet; `n/a`,
        \\no rstd target. A cell is **bold** where upstream is tuned and rstd has
        \\nothing (a gap); cores rstd does not know are left out of the gaps.
        \\
        \\{s}
    , .{ fam.dir, map.pin, if (fam.note.len > 0) try arena.print("{s}\n\n", .{fam.note}) else "" });

    var gaps: std.ArrayList(u8) = .empty;
    var n_gaps: usize = 0;

    for (arches) |ar| {
        const cores = try map.cores(ar.name) orelse continue;
        const known = try coreNames(arena, ar.cores);
        try w.print(arena, "## {s}\n\n", .{ar.name});
        var gap_cols: std.array_hash_map.String(std.ArrayList([]const u8)) = .empty;
        for (fam.tables) |t| {
            // every row first, to fold a table of defaults into a line
            var rows: std.ArrayList(u8) = .empty;
            var all_generic = true;
            for (cores) |core| {
                try rows.print(arena, "| {s} |", .{core});
                for (t.cols) |c| {
                    const v = try l3Sources(ar, core, c) orelse {
                        try rows.appendSlice(arena, " - |");
                        continue;
                    };
                    // the small-matrix kernels' default is `kernel/generic`'s own file
                    const defaults = c.ported != null and for (v) |x| {
                        if (!std.mem.startsWith(u8, x, "../generic/")) break false;
                    } else true;
                    if (v.len == 0 or defaults) {
                        try rows.appendSlice(arena, " generic |");
                        continue;
                    }
                    if (for (v) |x| {
                        if (x.len > 0) break false;
                    } else true) {
                        try rows.appendSlice(arena, " - |");
                        continue;
                    }
                    all_generic = false;
                    // an emptied variable is the default under `ifndef`
                    // (`KERNEL.SAPPHIRERAPIDS`'s `SBGEMM_SMALL_K_*`)
                    var named: std.ArrayList([]const u8) = .empty;
                    for (v) |x| if (x.len > 0 or c.ported == null) try named.append(arena, x);
                    const up = try upstreamCell(arena, ar, core, named.items);
                    const text = if (named.items.len < v.len) try arena.print("{s} \\| generic", .{up.text}) else up.text;
                    const served = if (c.ported) |p| for (v) |x| {
                        if (x.len > 0 and !contains(p, std.Io.Dir.path.basename(x))) break false;
                    } else true else c.rstd != .none;
                    const gap = up.tuned and ar.dir.len > 0 and !served and
                        contains(known, try std.ascii.allocLowerString(arena, core));
                    try rows.print(arena, " {s}{s}{s} |", .{ if (gap) "**" else "", text, if (gap) "**" else "" });
                    if (gap) {
                        n_gaps += 1;
                        const e = try gap_cols.getOrPut(arena, c.label);
                        if (!e.found_existing) e.value_ptr.* = .empty;
                        try e.value_ptr.append(arena, core);
                    }
                }
                try rows.append(arena, '\n');
            }
            if (all_generic) {
                try w.print(arena, "The {s} are `generic` on every core.\n\n", .{t.name});
                continue;
            }
            try w.print(arena, "### {s} {s}\n\n| core |", .{ ar.name, t.name });
            for (t.cols) |c| try w.print(arena, " {s} |", .{c.label});
            try w.appendSlice(arena, "\n|---|");
            for (t.cols) |_| try w.appendSlice(arena, "---|");
            try w.appendSlice(arena, "\n| rstd |");
            for (t.cols) |c| try w.print(arena, " {s} |", .{if (ar.dir.len == 0) "n/a" else @tagName(c.rstd)});
            try w.append(arena, '\n');
            try w.appendSlice(arena, rows.items);
            try w.append(arena, '\n');
        }
        for (gap_cols.keys(), gap_cols.values()) |label, list| try gaps.print(arena, "| {s} | {s} | {s} |\n", .{ ar.name, label, try std.mem.join(arena, ", ", list.items) });
    }

    try w.appendSlice(arena, "## Gaps\n\n");
    if (gaps.items.len == 0) {
        try w.appendSlice(arena, "None: the tile serves every tuned upstream kernel here.\n");
    } else {
        try w.print(arena, "{d} (core, column) cells where upstream is tuned and rstd is not.\n\n| architecture | column | cores |\n|---|---|---|\n", .{n_gaps});
        try w.appendSlice(arena, gaps.items);
    }

    const dest = try arena.print(blas ++ "/blas3/{s}/coverage.md", .{fam.dir});
    try Io.Dir.cwd().writeFile(init.io, .{ .sub_path = dest, .data = try alignTables(arena, w.items) });
    std.debug.print("wrote {s}: {d} gaps\n", .{ dest, n_gaps });
}

/// `text` with every table's columns padded to one width, which markdownlint's
/// aligned table style (MD060) wants.
fn alignTables(arena: std.mem.Allocator, text: []const u8) ![]const u8 {
    var out: std.ArrayList(u8) = .empty;
    var lines = std.mem.splitScalar(u8, text, '\n');
    var block: std.ArrayList([]const u8) = .empty;
    while (lines.next()) |line| {
        if (line.len > 0 and line[0] == '|') {
            try block.append(arena, line);
            continue;
        }
        try flushTable(arena, &out, block.items);
        block.clearRetainingCapacity();
        try out.appendSlice(arena, line);
        if (lines.peek() != null) try out.append(arena, '\n');
    }
    try flushTable(arena, &out, block.items);
    return out.items;
}

fn splitCells(arena: std.mem.Allocator, line: []const u8) ![]const []const u8 {
    var cells: std.ArrayList([]const u8) = .empty;
    var start: usize = 1;
    var i: usize = 1;
    while (i < line.len) : (i += 1) {
        if (line[i] == '\\') {
            i += 1;
        } else if (line[i] == '|') {
            try cells.append(arena, std.mem.trim(u8, line[start..i], " "));
            start = i + 1;
        }
    }
    return cells.items;
}

/// `c`'s width in characters, as markdownlint counts them.
fn width(c: []const u8) usize {
    return std.unicode.utf8CountCodepoints(c) catch c.len;
}

fn flushTable(arena: std.mem.Allocator, out: *std.ArrayList(u8), rows: []const []const u8) !void {
    if (rows.len == 0) return;
    var split: std.ArrayList([]const []const u8) = .empty;
    var widths: std.ArrayList(usize) = .empty;
    for (rows, 0..) |row, r| {
        const cells = try splitCells(arena, row);
        try split.append(arena, cells);
        for (cells, 0..) |c, k| {
            if (k >= widths.items.len) try widths.append(arena, 3);
            if (r != 1) widths.items[k] = @max(widths.items[k], width(c));
        }
    }
    for (split.items, 0..) |cells, r| {
        try out.append(arena, '|');
        for (cells, 0..) |c, k| {
            try out.append(arena, ' ');
            if (r == 1) {
                try out.appendNTimes(arena, '-', widths.items[k]);
            } else {
                try out.appendSlice(arena, c);
                try out.appendNTimes(arena, ' ', widths.items[k] - width(c));
            }
            try out.appendSlice(arena, " |");
        }
        try out.append(arena, '\n');
    }
}

/// One level-2 table: its routine and the make variable of each column.
const L2Family = struct { name: []const u8, cols: []const [2][]const u8 };

const l2_families = [_]L2Family{
    .{ .name = "gemv", .cols = &.{
        .{ "S n", "SGEMVNKERNEL" }, .{ "S t", "SGEMVTKERNEL" },
        .{ "D n", "DGEMVNKERNEL" }, .{ "D t", "DGEMVTKERNEL" },
        .{ "C n", "CGEMVNKERNEL" }, .{ "C t", "CGEMVTKERNEL" },
        .{ "Z n", "ZGEMVNKERNEL" }, .{ "Z t", "ZGEMVTKERNEL" },
        .{ "Q n", "QGEMVNKERNEL" }, .{ "Q t", "QGEMVTKERNEL" },
        .{ "X n", "XGEMVNKERNEL" }, .{ "X t", "XGEMVTKERNEL" },
    } },
    .{ .name = "ger", .cols = &.{
        .{ "S", "SGERKERNEL" },    .{ "D", "DGERKERNEL" },    .{ "Q", "QGERKERNEL" },
        .{ "C u", "CGERUKERNEL" }, .{ "C c", "CGERCKERNEL" }, .{ "Z u", "ZGERUKERNEL" },
        .{ "Z c", "ZGERCKERNEL" }, .{ "X u", "XGERUKERNEL" }, .{ "X c", "XGERCKERNEL" },
    } },
    .{ .name = "symv", .cols = &.{
        .{ "S U", "SSYMV_U_KERNEL" }, .{ "S L", "SSYMV_L_KERNEL" },
        .{ "D U", "DSYMV_U_KERNEL" }, .{ "D L", "DSYMV_L_KERNEL" },
        .{ "Q U", "QSYMV_U_KERNEL" }, .{ "Q L", "QSYMV_L_KERNEL" },
        .{ "C U", "CSYMV_U_KERNEL" }, .{ "C L", "CSYMV_L_KERNEL" },
        .{ "Z U", "ZSYMV_U_KERNEL" }, .{ "Z L", "ZSYMV_L_KERNEL" },
        .{ "X U", "XSYMV_U_KERNEL" }, .{ "X L", "XSYMV_L_KERNEL" },
    } },
    .{ .name = "hemv", .cols = &.{
        .{ "C U", "CHEMV_U_KERNEL" }, .{ "C L", "CHEMV_L_KERNEL" }, .{ "C V", "CHEMV_V_KERNEL" }, .{ "C M", "CHEMV_M_KERNEL" },
        .{ "Z U", "ZHEMV_U_KERNEL" }, .{ "Z L", "ZHEMV_L_KERNEL" }, .{ "Z V", "ZHEMV_V_KERNEL" }, .{ "Z M", "ZHEMV_M_KERNEL" },
        .{ "X U", "XHEMV_U_KERNEL" }, .{ "X L", "XHEMV_L_KERNEL" }, .{ "X V", "XHEMV_V_KERNEL" }, .{ "X M", "XHEMV_M_KERNEL" },
    } },
    .{ .name = "half-precision gemv", .cols = &.{
        .{ "SB n", "SBGEMVNKERNEL" }, .{ "SB t", "SBGEMVTKERNEL" },
        .{ "B n", "BGEMVNKERNEL" },   .{ "B t", "BGEMVTKERNEL" },
        .{ "SH n", "SHGEMVNKERNEL" }, .{ "SH t", "SHGEMVTKERNEL" },
    } },
};

/// What a C source holds besides plain C: the `.c` files it `#include`s
/// (OpenBLAS's per-core microkernels), and whether it, or one of those, has
/// inline assembly or vector intrinsics. A text search, so a heuristic.
const CInfo = struct {
    inline_asm: bool = false,
    intrinsics: bool = false,
    /// `name (asm)` / `name (intrinsics)` / `name (C)` per included `.c`.
    includes: std.ArrayList([]const u8) = .empty,
};

fn hasAsm(t: []const u8) bool {
    for ([_][]const u8{ "__asm__", "asm volatile", "__asm ", "__asm(", " asm(" }) |n| if (std.mem.find(u8, t, n) != null) return true;
    return false;
}

fn hasIntrinsics(t: []const u8) bool {
    for ([_][]const u8{ "immintrin.h", "arm_neon.h", "arm_sve.h", "arm_sme.h", "altivec.h", "riscv_vector.h", "lsxintrin.h", "lasxintrin.h", "vecintrin.h", "msa.h", "wasm_simd128.h", "__vector", "_mm256_", "_mm512_", "_mm_", "vld1q_", "svld1", "__riscv_v", "RISCV_RVV(", "__lsx_", "__lasx_", "vec_xl", "vec_madd" }) |n| if (std.mem.find(u8, t, n) != null) return true;
    return false;
}

fn classifyC(arena: std.mem.Allocator, arch: []const u8, source: []const u8) !CInfo {
    const p = try map.path(&.{ map.root, "kernel", arch, source });
    var info: CInfo = .{};
    const text = map.read(p) orelse {
        try info.includes.append(arena, "(missing)");
        return info;
    };
    info.inline_asm = hasAsm(text);
    info.intrinsics = hasIntrinsics(text);
    const dir = std.Io.Dir.path.dirname(p) orelse ".";
    var lines = std.mem.splitScalar(u8, text, '\n');
    while (lines.next()) |raw| {
        const line = std.mem.trim(u8, raw, " \t\r");
        if (!std.mem.startsWith(u8, line, "#include")) continue;
        const q = std.mem.findScalar(u8, line, '"') orelse continue;
        const e = std.mem.findScalarPos(u8, line, q + 1, '"') orelse continue;
        const name = line[q + 1 .. e];
        if (!std.mem.endsWith(u8, name, ".c")) continue;
        const sub = map.read(try map.path(&.{ dir, name }));
        const kind = if (sub) |t| (if (hasAsm(t)) "asm" else if (hasIntrinsics(t)) "intrinsics" else "C") else "missing";
        const entry = try arena.print("{s} ({s})", .{ name, kind });
        for (info.includes.items) |o| {
            if (std.mem.eql(u8, o, entry)) break;
        } else try info.includes.append(arena, entry);
    }
    return info;
}

const Use = struct {
    roles: std.array_hash_map.String(void) = .empty,
    cores: std.array_hash_map.String(void) = .empty,
};

/// The slots of `Vtable.zig`: its `<slot>: ?*const fn` fields.
fn vtableSlots(arena: std.mem.Allocator) ![]const []const u8 {
    const text = map.read(blas ++ "/Vtable.zig") orelse map.die("cannot read Vtable.zig", .{});
    var out: std.ArrayList([]const u8) = .empty;
    var lines = std.mem.splitScalar(u8, text, '\n');
    while (lines.next()) |line| {
        const at = std.mem.find(u8, line, ": ?*const fn") orelse continue;
        try out.append(arena, line[0..at]);
    }
    return out.items;
}

fn l2(init: std.process.Init, arena: std.mem.Allocator, dest: []const u8) !void {
    var w: std.ArrayList(u8) = .empty;
    try w.print(arena,
        \\# Level-2 kernel coverage
        \\
        \\Generated by `tools/blas_asm/survey.zig l2` from OpenBLAS `kernel/` (`{s}`) and
        \\rstd's kernel indexes; re-derive it rather than editing it. Each cell is
        \\the source that core's `KERNEL.<core>` selects for the column's make
        \\variable, falling back to the architecture's `KERNEL`, then to the
        \\`Makefile.L2` default (`|` separates alternatives chosen by a
        \\conditional), and after the arrow what rstd serves for that core, for
        \\every role built from it: the complex `gemv` forms `r`/`o`/`s` build the
        \\`n` source and `c`/`u`/`d` the `t` source (`CONJ`/`XCONJ`). `ger`'s
        \\`v`/`d` forms, built from the `u`/`c` source with `XCONJ`, have no
        \\vtable slot (`Vtable.zig`) and are left out. The half-precision kernels
        \\are built only with `BUILD_BFLOAT16`/`BUILD_HFLOAT16`; rstd's are
        \\`ext/bf16/kernel.zig`'s, the same on every core (`portable` here).
        \\
        \\
    ++ legend ++ "\n", .{map.pin});

    const text = map.read(try map.path(&.{ map.root, "kernel", "Makefile.L2" })).?;
    const all_roles = try roles.parse(arena, text);
    const slots = try vtableSlots(arena);
    const fills = try vectorSlots(arena);
    const displacing = try readCores(arena, "displacing", true);
    const flushing = try readCores(arena, "flushing", false);
    var t: Tally = .{};
    // every distinct (architecture, source): its roles and cores
    var seen: std.array_hash_map.String(Use) = .empty;

    for (arches) |ar| {
        const arch = ar.name;
        const cores = try map.cores(arch) orelse continue;
        const known = try coreNames(arena, ar.cores);
        const routing: Routing = .{ .rule = try readRule(arena, ar), .displacing = displacing, .flushing = flushing };
        try w.print(arena, "## {s}\n\n", .{arch});
        // each core's gaps over every table
        var core_gaps: std.array_hash_map.String(std.ArrayList([]const u8)) = .empty;
        if (ar.dir.len == 0) try w.appendSlice(arena, "rstd has no target here.\n\n");
        for (l2_families) |fam| {
            const s = try readServed(arena, ar, "blas2", try roles.dirFor(arena, fam.cols[0][1]));
            // the roles of each column that have a slot
            var groups: std.ArrayList([]const roles.Role) = .empty;
            for (fam.cols) |c| {
                var g: std.ArrayList(roles.Role) = .empty;
                for (all_roles) |r| {
                    if (std.mem.eql(u8, r.key, c[1]) and r.slot != null and contains(slots, r.slot.?)) try g.append(arena, r);
                }
                try groups.append(arena, g.items);
            }
            try w.print(arena, "### {s} {s}\n\n", .{ arch, fam.name });
            if (ar.dir.len > 0) try note(arena, &w, ar, s);
            try w.appendSlice(arena, "| core |");
            for (fam.cols) |c| try w.print(arena, " {s} |", .{c[0]});
            try w.appendSlice(arena, "\n|---|");
            for (fam.cols) |_| try w.appendSlice(arena, "---|");
            try w.append(arena, '\n');
            for (cores) |core| {
                const lower = try std.ascii.allocLowerString(arena, core);
                const here = ar.dir.len > 0 and contains(known, lower);
                const cg = try core_gaps.getOrPut(arena, core);
                if (!cg.found_existing) cg.value_ptr.* = .empty;
                try w.print(arena, "| {s} |", .{core});
                for (fam.cols, groups.items) |c, group| {
                    const v = try map.resolveKey(arch, core, c[1], "Makefile.L2");
                    for (v) |src| {
                        const e = try seen.getOrPut(arena, try arena.print("{s}/{s}", .{ arch, src }));
                        if (!e.found_existing) e.value_ptr.* = .{};
                        try e.value_ptr.roles.put(arena, c[1], {});
                        try e.value_ptr.cores.put(arena, core, {});
                    }
                    const up = try upstreamCell(arena, ar, core, v);
                    const rstd = try t.cell(arena, ar, s, routing, lower, here, group, up.tuned, fills, cg.value_ptr);
                    try w.print(arena, " {s} → {s} |", .{ up.text, rstd });
                }
                try w.append(arena, '\n');
            }
            try w.append(arena, '\n');
        }
        for (core_gaps.keys(), core_gaps.values()) |core, list| try t.row(arena, arch, core, list.items);
    }
    try t.write(arena, &w);
    try w.append(arena, '\n');

    // the sources, by kind
    try w.appendSlice(arena,
        \\## Sources
        \\
        \\Every distinct (architecture, source). `kind` is `S` for assembly; for C,
        \\`generic` (plain C), `intrinsics` or `inline asm` in the file itself, or
        \\`C + includes`; `includes` lists the `.c` files it includes, each marked by
        \\what it holds (a text search: `asm` means inline assembly).
        \\
        \\| architecture | source | kind | includes | roles | cores |
        \\|---|---|---|---|---|---|
        \\
    );
    var n_asm: usize = 0;
    var n_c: usize = 0;
    var n_c_asm: usize = 0;
    var n_c_intr: usize = 0;
    var n_missing: usize = 0;
    for (seen.keys(), seen.values()) |k, u| {
        const slash = std.mem.findScalar(u8, k, '/').?;
        const arch = k[0..slash];
        const source = k[slash + 1 ..];
        var kind: []const u8 = "S";
        var inc: []const u8 = "";
        // a default no file backs: a precision this architecture never builds
        if (map.read(try map.path(&.{ map.root, "kernel", arch, source })) == null) {
            kind = "missing";
            n_missing += 1;
        } else if (isAsm(source)) n_asm += 1 else {
            n_c += 1;
            const info = try classifyC(arena, arch, source);
            var any_asm = info.inline_asm;
            var any_intr = info.intrinsics;
            for (info.includes.items) |i| {
                if (std.mem.endsWith(u8, i, "(asm)")) any_asm = true;
                if (std.mem.endsWith(u8, i, "(intrinsics)")) any_intr = true;
            }
            if (any_asm) n_c_asm += 1 else if (any_intr) n_c_intr += 1;
            kind = if (info.inline_asm) "inline asm" else if (info.intrinsics) "intrinsics" else if (info.includes.items.len > 0) "C + includes" else "generic";
            inc = try join(info.includes.items);
        }
        try w.print(arena, "| {s} | `{s}` | {s} | {s} | {s} | {s} |\n", .{ arch, source, kind, inc, try std.mem.join(arena, ", ", u.roles.keys()), try std.mem.join(arena, ", ", u.cores.keys()) });
    }
    try w.print(arena, "\n{d} distinct (architecture, source) pairs: {d} assembly, {d} C ({d} with inline assembly in the file or an included `.c`, {d} more with intrinsics only), {d} missing.\n", .{ seen.count(), n_asm, n_c, n_c_asm, n_c_intr, n_missing });

    try Io.Dir.cwd().writeFile(init.io, .{ .sub_path = dest, .data = try alignTables(arena, w.items) });
    std.debug.print("wrote {s}: {d} gaps; {d} sources ({d} asm, {d} C, {d} C with asm)\n", .{ dest, t.n_gaps, seen.count(), n_asm, n_c, n_c_asm });
}
