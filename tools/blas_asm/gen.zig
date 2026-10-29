//! Converts OpenBLAS's level-1 and level-2 `.S` kernels into Zig files whose
//! bodies are the same assembly, as naked functions.
//!
//!     OPENBLAS=<checkout> zig run tools/blas_asm/gen.zig [-- <arch|level|routine>...]
//!
//! With nothing named, converts every routine for every architecture that
//! `arch/<name>.zig` describes; otherwise only the architectures, levels
//! (`blas1`, `blas2`) and routines (`blas<n>/<routine>/` names, such as
//! `axpy` or `gemv`) named. Output goes under `$BLAS_OUT` (default
//! `src/linalg/blas`), as `blas<n>/<routine>/<arch>/asm/`.
//!
//! The source is preprocessed once per core and variant, since OpenBLAS
//! selects code with `#ifdef <CORE>` (and the ABI with `WINDOWS_ABI` and the
//! like). Which cores use a source, and for which precisions, comes from the
//! `KERNEL.<core>` files (`kernel_map.zig`). Cores whose expansions are
//! identical, or nearly (`similar`), share one file.
//!
//! Each function is the source's own, unchanged: it reads its arguments where
//! the calling convention puts them, saves what it uses, and returns. Only the
//! directives that declare the function or its section go (Zig places it),
//! tabs become spaces (Zig's tokenizer requires it inside `\`), named labels
//! become numeric local ones (one object holds many of these), and x86's
//! `.align` becomes `.balign` (Mach-O reads `.align` as a power of two).
//! Before that, `arch/<name>.zig`'s patches fix upstream bugs; its `retired`
//! ones, for sources no longer converted, are only checked against upstream,
//! and a full run stops on a live patch that applied to nothing.
//!
//! `asm/__root__.zig` lists which file fills which vtable slot; each file
//! says which cores use it. The kernels `roles.skip` names, which no core
//! would install, are not converted; the index lists the cores of those it
//! keeps them for (`skipped`).

const std = @import("std");
const Io = std.Io;
const map = @import("kernel_map.zig");
const arch = @import("arch/__root__.zig");
const roles = @import("roles.zig");

const Arch = arch.Arch;
const Kernel = arch.Kernel;
const Prec = arch.Prec;

const arches = [_]Arch{
    @import("arch/x86_64.zig").x86_64,
    @import("arch/arm64.zig").arm64,
    @import("arch/loongarch64.zig").loongarch64,
    @import("arch/mips64.zig").mips64,
    @import("arch/arm.zig").arm,
    @import("arch/power.zig").power,
    @import("arch/x86.zig").x86,
    @import("arch/sparc.zig").sparc,
};

/// OpenBLAS's names: `daxpy_atom` for the double kernel of `axpy_atom.S`,
/// `zaxpyc_sse2` for the conjugated complex one of `zaxpy_sse2.S`, then the
/// build variant and, where cores differ within one, the first core.
fn fileName(a: std.mem.Allocator, k: Kernel, width: []const u8, first: ?[]const u8) ![]const u8 {
    const stem = k.stem();
    // what follows the routine's name in the source's: `_sse2` of `dot_sse2`
    var rest = if (std.mem.findScalar(u8, stem, '_')) |at| stem[at..] else "";
    // and its form, which the role (`zgemv_r`, `zhemv_U`) already names:
    // `_atom` of `zgemv_n_atom`, `_sse2` of `zsymv_U_sse2`
    if (std.mem.findScalar(u8, k.role_name, '_') != null and rest.len >= 2 and
        std.mem.findScalar(u8, "ntUL", rest[1]) != null and (rest.len == 2 or rest[2] == '_'))
        rest = rest[2..];
    const core = if (first) |f| try a.print("_{s}", .{f}) else "";
    return a.print("{s}{s}{s}{s}.zig", .{ try k.role(a), rest, width, core });
}

/// Tabs become spaces to the next multiple of eight.
fn expandTabs(a: std.mem.Allocator, s: []const u8) ![]const u8 {
    var out: std.ArrayList(u8) = .empty;
    var col: usize = 0;
    for (s) |c| {
        if (c == '\t') {
            const n = 8 - col % 8;
            try out.appendNTimes(a, ' ', n);
            col += n;
        } else {
            try out.append(a, c);
            col += 1;
        }
    }
    return out.items;
}

/// Runs the preprocessor over the source the way OpenBLAS's build does.
///
/// `-DDYNAMIC_ARCH`: this project always builds that way (one object per
/// core, like OpenBLAS's own `DYNAMIC_ARCH=1`), so sources that branch on it
/// should convert as they would for that build, not a single-target one.
fn runCpp(a: std.mem.Allocator, ar: Arch, inc: []const u8, k: Kernel, core: []const u8, v: arch.Variant, dump_macros: bool) ![]const u8 {
    var argv: std.ArrayList([]const u8) = .empty;
    try argv.appendSlice(a, &.{ "zig", "cc", "-E", "-P", "-x", "assembler-with-cpp", "-w", "-DDYNAMIC_ARCH", "-I", inc, "-I", map.root, "-I", try map.path(&.{ map.root, "kernel", ar.name }), "-target", v.triple orelse ar.triple });
    try argv.appendSlice(a, ar.defines);
    if (k.prec.double()) try argv.append(a, "-DDOUBLE");
    if (k.prec.complex()) try argv.append(a, "-DCOMPLEX");
    if (k.conj) try argv.append(a, "-DCONJ");
    try argv.appendSlice(a, k.defines);
    try argv.appendSlice(a, try ar.core_defines(a, k, core));
    try argv.appendSlice(a, v.defines);
    if (dump_macros) try argv.append(a, "-dM");
    try argv.append(a, try patched(a, ar, inc, k.source, v));
    try argv.appendSlice(a, &.{ "-o", "-" });

    const r = std.process.run(a, map.io, .{ .argv = argv.items }) catch |e|
        map.die("cannot run the preprocessor: {t}", .{e});
    if (r.stderr.len > 0) std.debug.print("{s}", .{r.stderr});
    return r.stdout;
}

/// The path of `source` with our patches for variant `v` applied: a copy under
/// `inc` when any apply (its includes still resolve through `-I`), else the
/// upstream file.
fn patched(a: std.mem.Allocator, ar: Arch, inc: []const u8, source: []const u8, v: arch.Variant) ![]const u8 {
    const upstream = try map.path(&.{ map.root, "kernel", ar.name, source });
    var text: ?[]const u8 = null;
    for ([_][]const arch.Patch{ ar.patches, v.patches }) |list| for (list) |*p| {
        if (!std.mem.eql(u8, p.source, source)) continue;
        text = try apply(a, text orelse map.read(upstream) orelse map.die("cannot read {s}", .{upstream}), p.*);
        try applied.put(a, p, {});
    };
    const t = text orelse return upstream;
    const dir = try map.path(&.{ inc, "patched", if (v.name.len == 0) "_" else v.name });
    try Io.Dir.cwd().createDirPath(map.io, dir);
    const out = try map.path(&.{ dir, source });
    try Io.Dir.cwd().writeFile(map.io, .{ .sub_path = out, .data = t });
    return out;
}

/// `text` with `p` applied, which must match `p.count` times.
fn apply(a: std.mem.Allocator, text: []const u8, p: arch.Patch) ![]const u8 {
    const n = std.mem.count(u8, text, p.find);
    if (n != p.count) map.die("patch for {s} matches {d} times, not {d}:\n{s}", .{ p.source, n, p.count, p.find });
    return std.mem.replaceOwned(u8, a, text, p.find, p.replace);
}

/// The patches `patched` applied, by address: a full run stops on any other.
var applied: std.array_hash_map.Auto(*const arch.Patch, void) = .empty;

/// Stops on a retired patch that no longer matches upstream; each source's
/// apply in order, as `patched` applies live ones.
fn checkRetired(a: std.mem.Allocator, ar: Arch) !void {
    for (ar.retired, 0..) |p, i| {
        const first = for (ar.retired[0..i]) |q| {
            if (std.mem.eql(u8, q.source, p.source)) break false;
        } else true;
        if (!first) continue;
        const upstream = try map.path(&.{ map.root, "kernel", ar.name, p.source });
        var t = map.read(upstream) orelse map.die("cannot read {s}", .{upstream});
        for (ar.retired[i..]) |q| if (std.mem.eql(u8, q.source, p.source)) {
            t = try apply(a, t, q);
        };
    }
}

/// After a full run: stops on a live patch no converted source took, which
/// belongs in `retired`.
fn checkApplied(ar: Arch) void {
    for (ar.patches) |*p| if (!applied.contains(p)) unapplied(ar, p);
    for (ar.variants) |v| for (v.patches) |*p| if (!applied.contains(p)) unapplied(ar, p);
}

fn unapplied(ar: Arch, p: *const arch.Patch) noreturn {
    map.die("{s}: the patch for {s} applies to nothing converted; move it to `retired`", .{ ar.name, p.source });
}

/// The preprocessed source, tab-free and without empty lines.
fn preprocess(a: std.mem.Allocator, ar: Arch, inc: []const u8, k: Kernel, core: []const u8, v: arch.Variant) ![]const []const u8 {
    const text = try runCpp(a, ar, inc, k, core, v, false);
    var out: std.ArrayList([]const u8) = .empty;
    var it = std.mem.splitScalar(u8, text, '\n');
    while (it.next()) |raw| {
        const line = std.mem.trimEnd(u8, try expandTabs(a, raw), " \r");
        if (std.mem.trim(u8, line, " ").len == 0) continue;
        try out.append(a, line);
    }
    return out.items;
}

fn lowered(a: std.mem.Allocator, s: []const u8) ![]const u8 {
    return try std.ascii.allocLowerString(a, s);
}

fn put(b: *std.ArrayList(u8), a: std.mem.Allocator, comptime fmt: []const u8, args: anytype) !void {
    try b.print(a, fmt, args);
}

/// Whether `c` can be part of a label or macro name.
fn isName(c: u8) bool {
    return std.ascii.isAlphanumeric(c) or c == '_' or c == '.' or c == '$';
}

/// The end of the name token starting at `s[i]`, one past a trailing `\@`
/// (GAS's per-macro-invocation marker), which makes it part of the name.
fn nameEnd(s: []const u8, i: usize) usize {
    var e = i;
    while (true) {
        while (e < s.len and isName(s[e])) e += 1;
        // a macro's parameter spliced into the name (`.L_\XW\()_N_7`, in
        // loongarch64's gemv): one name, so it becomes a numeric label too
        if (std.mem.startsWith(u8, s[e..], "\\()")) {
            e += 3;
        } else if (e > i and e + 1 < s.len and s[e] == '\\' and (std.ascii.isAlphabetic(s[e + 1]) or s[e + 1] == '_')) {
            e += 1;
        } else break;
    }
    if (e + 1 < s.len and s[e] == '\\' and s[e + 1] == '@') e += 2;
    return e;
}

/// The named label a line starts by defining (`name:`), if any. Numeric
/// labels are already local and stay as the source has them.
fn definedLabel(line: []const u8) ?[]const u8 {
    const t = std.mem.trim(u8, line, " ");
    const e = nameEnd(t, 0);
    if (e == 0 or e >= t.len or t[e] != ':' or std.ascii.isDigit(t[0])) return null;
    return t[0..e];
}

/// The function's own label; Zig declares the function.
fn isOwnLabel(name: []const u8) bool {
    return std.mem.eql(u8, name, "k") or std.mem.eql(u8, name, "k_");
}

/// The first number a named label becomes, above any the sources use.
const label_base = 20000;

/// The source's lines as a function body. The function's own label goes, and
/// every named label becomes a numeric one, `names` giving its number (the
/// same across the cores of one file, so their shared lines stay equal): the
/// assembler keeps numeric labels local however many functions in one object
/// define them, so each core's function can be instantiated beside another's.
/// The macros the source defines are purged again at the end.
fn finish(a: std.mem.Allocator, lines: []const []const u8, names: *const std.array_hash_map.String(void)) ![]const []const u8 {
    // which `.macro` ... `.endm` each line is in (0 for none), and where each
    // label is defined
    const region = try a.alloc(usize, lines.len);
    var macros: std.ArrayList([]const u8) = .empty;
    var defs: std.StringHashMapUnmanaged(std.ArrayList(usize)) = .empty;
    var cur: usize = 0;
    for (lines, 0..) |line, i| {
        const t = std.mem.trim(u8, line, " ");
        if (std.mem.startsWith(u8, t, ".macro ")) {
            var it = std.mem.tokenizeAny(u8, t[".macro ".len..], " ,");
            if (it.next()) |name| try macros.append(a, name);
            cur = macros.items.len;
        }
        region[i] = cur;
        if (std.mem.eql(u8, t, ".endm")) cur = 0;
        if (definedLabel(line)) |n| {
            const e = try defs.getOrPut(a, n);
            if (!e.found_existing) e.value_ptr.* = .empty;
            try e.value_ptr.append(a, i);
        }
    }

    var out: std.ArrayList([]const u8) = .empty;
    for (lines, 0..) |line, li| {
        if (definedLabel(line)) |n| if (isOwnLabel(n)) continue;
        const start = line.len - std.mem.trimStart(u8, line, " ").len;
        var buf: std.ArrayList(u8) = .empty;
        var i: usize = 0;
        while (i < line.len) {
            if (isName(line[i]) and (i == 0 or !isName(line[i - 1]))) {
                const e = nameEnd(line, i);
                const tok = line[i..e];
                if (names.getIndex(tok)) |num| {
                    try buf.print(a, "{d}", .{label_base + num});
                    // a reference names the definition by direction: in a
                    // macro, the one in that macro's text, else the one outside
                    if (!(i == start and e < line.len and line[e] == ':')) {
                        var fwd = false;
                        var found = false;
                        for (defs.get(tok).?.items) |d| {
                            if (region[d] != region[li]) continue;
                            found = true;
                            if (d > li) fwd = true;
                        }
                        if (!found) fwd = outsideDirection(lines, region, macros.items, defs.get(tok).?.items, region[li]) orelse
                            std.debug.panic("label {s} is referenced across a macro's boundary: {s}", .{ tok, line });
                        try buf.append(a, if (fwd) 'f' else 'b');
                    }
                } else try buf.appendSlice(a, tok);
                i = e;
            } else {
                try buf.append(a, line[i]);
                i += 1;
            }
        }
        try out.append(a, buf.items);
    }
    for (macros.items) |m| try out.append(a, try a.print(" .purgem {s}", .{m}));
    return out.items;
}

/// For a macro (`region`) jumping to a label defined once outside any macro
/// (`.L_END` in loongarch64's `sgemv_n_8_lasx.S`): forward when every
/// top-level invocation precedes that label, back when every one follows it,
/// null when they straddle it or there is no such invocation.
fn outsideDirection(lines: []const []const u8, region: []const usize, macros: []const []const u8, defs: []const usize, macro: usize) ?bool {
    if (macro == 0) return null;
    var def: ?usize = null;
    for (defs) |d| if (region[d] == 0) {
        if (def != null) return null;
        def = d;
    };
    const at = def orelse return null;
    const name = macros[macro - 1];
    var before = false;
    var after = false;
    for (lines, 0..) |l, i| {
        if (region[i] != 0) continue;
        const t = std.mem.trim(u8, l, " ");
        if (!std.mem.startsWith(u8, t, name) or (t.len > name.len and t[name.len] != ' ')) continue;
        if (i < at) before = true else after = true;
    }
    if (before == after) return null;
    return before;
}

const Op = union(enum) { eq: [2]usize, del: usize, ins: usize };

/// The shortest edit script from `x` to `y` (Myers), or null when it takes
/// more than `max_d` insertions and deletions.
fn diff(a: std.mem.Allocator, x: []const []const u8, y: []const []const u8, max_d: usize) !?[]const Op {
    var scratch: std.heap.ArenaAllocator = .init(std.heap.page_allocator);
    defer scratch.deinit();
    const s = scratch.allocator();
    const n: isize = @intCast(x.len);
    const m: isize = @intCast(y.len);
    const max: isize = @intCast(@min(max_d, x.len + y.len));
    const off = max + 1;
    const v = try s.alloc(isize, @intCast(2 * max + 3));
    @memset(v, 0);
    // `trace[d]` is `v[-d - 1 .. d + 1]` before step `d`
    var trace: std.ArrayList([]const isize) = .empty;
    var d: isize = 0;
    found: while (d <= max) : (d += 1) {
        try trace.append(s, try s.dupe(isize, v[@intCast(off - d - 1)..@intCast(off + d + 2)]));
        var k = -d;
        while (k <= d) : (k += 2) {
            var px: isize = if (k == -d or (k != d and v[@intCast(off + k - 1)] < v[@intCast(off + k + 1)]))
                v[@intCast(off + k + 1)]
            else
                v[@intCast(off + k - 1)] + 1;
            var py = px - k;
            while (px < n and py < m and std.mem.eql(u8, x[@intCast(px)], y[@intCast(py)])) {
                px += 1;
                py += 1;
            }
            v[@intCast(off + k)] = px;
            if (px >= n and py >= m) break :found;
        }
    } else return null;

    var rev: std.ArrayList(Op) = .empty;
    var px = n;
    var py = m;
    while (d > 0) : (d -= 1) {
        const t = trace.items[@intCast(d)];
        const at = struct {
            fn f(tr: []const isize, dd: isize, kk: isize) isize {
                return tr[@intCast(kk + dd + 1)];
            }
        }.f;
        const k = px - py;
        const down = k == -d or (k != d and at(t, d, k - 1) < at(t, d, k + 1));
        const pk = if (down) k + 1 else k - 1;
        const prev_x = at(t, d, pk);
        const prev_y = prev_x - pk;
        const mid_x = if (down) prev_x else prev_x + 1;
        while (px > mid_x) {
            px -= 1;
            py -= 1;
            try rev.append(s, .{ .eq = .{ @intCast(px), @intCast(py) } });
        }
        if (down) try rev.append(s, .{ .ins = @intCast(prev_y) }) else try rev.append(s, .{ .del = @intCast(prev_x) });
        px = prev_x;
        py = prev_y;
    }
    while (px > 0) {
        px -= 1;
        py -= 1;
        try rev.append(s, .{ .eq = .{ @intCast(px), @intCast(py) } });
    }
    const ops = try a.alloc(Op, rev.items.len);
    for (ops, 0..) |*o, i| o.* = rev.items[rev.items.len - 1 - i];
    return ops;
}

/// Whether two expansions are close enough to share a file: at most a tenth
/// of their lines differ.
fn similar(a: std.mem.Allocator, x: []const []const u8, y: []const []const u8) !bool {
    return (try diff(a, x, y, (x.len + y.len) / 10)) != null;
}

/// A line of a merged body, and which of the merged texts have it.
const MLine = struct { text: []const u8, mask: u128 };

/// One sequence that holds every text in order: each text is exactly the
/// lines whose mask has its bit.
fn merge(a: std.mem.Allocator, texts: []const []const []const u8) ![]const MLine {
    var m: std.ArrayList(MLine) = .empty;
    for (texts[0]) |l| try m.append(a, .{ .text = l, .mask = 1 });
    for (texts[1..], 1..) |t, i| {
        const bit = @as(u128, 1) << @intCast(i);
        const old = try a.alloc([]const u8, m.items.len);
        for (m.items, old) |l, *o| o.* = l.text;
        const ops = (try diff(a, old, t, std.math.maxInt(usize))).?;
        var next: std.ArrayList(MLine) = .empty;
        for (ops) |op| switch (op) {
            .eq => |p| try next.append(a, .{ .text = m.items[p[0]].text, .mask = m.items[p[0]].mask | bit }),
            .del => |ai| try next.append(a, m.items[ai]),
            .ins => |bi| try next.append(a, .{ .text = t[bi], .mask = bit }),
        };
        m = next;
    }
    return m.items;
}

/// One expansion of a source: a core, built as one variant.
const Text = struct { core: []const u8, variant: []const u8, lines: []const []const u8 };

/// The text of one Zig file: the merged expansions of `texts`, one naked
/// function per core and variant, the lines only some have under `parts`.
fn emit(a: std.mem.Allocator, ar: Arch, k: Kernel, texts: []const Text, lines: []const MLine) ![]const u8 {
    var w: std.ArrayList(u8) = .empty;
    const role = try k.role(a);
    const all: u128 = if (texts.len == 128) std.math.maxInt(u128) else (@as(u128, 1) << @intCast(texts.len)) - 1;

    // the masks of the lines not every text has, each a part
    var conds: std.array_hash_map.Auto(u128, void) = .empty;
    for (lines) |l| if (l.mask != all) try conds.put(a, l.mask, {});
    const part_count = conds.count();

    // header: the cores per variant
    var who: std.ArrayList(u8) = .empty;
    for (ar.variants) |v| {
        var first = true;
        for (texts) |t| if (std.mem.eql(u8, t.variant, v.name)) {
            if (first and who.items.len > 0) try who.appendSlice(a, "; ");
            try who.appendSlice(a, if (first) "" else ", ");
            try who.appendSlice(a, t.core);
            first = false;
        };
        if (!first and v.name.len > 0) try put(&who, a, " ({s})", .{v.name});
    }
    try put(&w, a,
        \\//! `{s}`, from OpenBLAS `kernel/{s}/{s}`, as expanded for {s}.
        \\//!
        \\//! Generated by `tools/blas_asm/gen.zig`: the assembly of that source,
        \\//! unchanged, as a naked function per core; the lines only some cores
        \\//! have are under `parts`. Re-derive it rather than editing it.
        \\
        \\const std = @import("std");
        \\
        \\/// Each function's core and variant.
        \\const keys = [_][2][]const u8{{
        \\
    , .{ role, ar.name, k.source, who.items });
    for (texts) |t| try put(&w, a, "    .{{ \"{s}\", \"{s}\" }},\n", .{ t.core, t.variant });
    try put(&w, a, "}};\n\n", .{});
    if (part_count > 0) {
        try put(&w, a, "/// Which parts each key's function has: `__bk`'s arguments.\nconst parts = [_][]const u8{{\n", .{});
        for (texts, 0..) |_, i| {
            try put(&w, a, "    \"", .{});
            for (conds.keys(), 0..) |c, j| try put(&w, a, "{s}{s}", .{ if (j == 0) "" else ", ", if (c & (@as(u128, 1) << @intCast(i)) != 0) "1" else "0" });
            try put(&w, a, "\",\n", .{});
        }
        try put(&w, a, "}};\n\n", .{});
    }
    try put(&w, a,
        \\/// Whether any build of `core` uses this file.
        \\pub fn serves(comptime core: []const u8) bool {{
        \\    inline for (keys) |key| if (comptime std.mem.eql(u8, key[0], core)) return true;
        \\    return false;
        \\}}
        \\
        \\/// The function for `core` built as `variant`, if this file has it.
        \\pub fn kernel(comptime core: []const u8, comptime variant: []const u8) ?type {{
        \\
    , .{});
    if (part_count > 0) {
        try w.appendSlice(a,
            \\    inline for (keys, 0..) |key, i| {
            \\        if (comptime std.mem.eql(u8, key[0], core) and std.mem.eql(u8, key[1], variant)) return Body(parts[i]);
            \\
        );
    } else {
        try w.appendSlice(a,
            \\    inline for (keys) |key| {
            \\        if (comptime std.mem.eql(u8, key[0], core) and std.mem.eql(u8, key[1], variant)) return Body;
            \\
        );
    }
    try w.appendSlice(a,
        \\    }
        \\    return null;
        \\}
        \\
        \\
    );
    if (part_count == 0) {
        try w.appendSlice(a, "const Body = struct {\npub fn entry() callconv(.naked) void {\nasm volatile (\n");
        for (lines) |l| try put(&w, a, "\\\\{s}\n", .{l.text});
        try w.appendSlice(a, "\\\\\n);\n}\n};\n");
        return w.items;
    }

    // The parts are chosen by the assembler, not by joining strings at
    // comptime, which costs minutes per file: the text is one macro, each
    // part between `.if \bpN` and `.endif`, and a function expands it with
    // its part flags. A macro, not symbols, so nothing reaches the object.
    try w.appendSlice(a,
        \\/// The functions' one text, a macro whose parameters select the parts.
        \\const text =
        \\
    );
    try w.appendSlice(a, "\\\\ .macro __bk");
    for (0..part_count) |j| try put(&w, a, "{s} bp{d}", .{ if (j == 0) "" else ",", j });
    try w.appendSlice(a, "\n");
    var i: usize = 0;
    while (i < lines.len) {
        const mask = lines[i].mask;
        var e = i;
        while (e < lines.len and lines[e].mask == mask) e += 1;
        if (mask != all) try put(&w, a, "\\\\ .if \\bp{d}\n", .{conds.getIndex(mask).?});
        for (lines[i..e]) |l| try put(&w, a, "\\\\{s}\n", .{l.text});
        if (mask != all) try w.appendSlice(a, "\\\\ .endif\n");
        i = e;
    }
    try w.appendSlice(a, "\\\\ .endm\n\\\\\n;\n\n");
    try w.appendSlice(a,
        \\fn Body(comptime args: []const u8) type {
        \\    return struct {
        \\        pub fn entry() callconv(.naked) void {
        \\            asm volatile (text ++ (" __bk " ++ args ++ "\n .purgem __bk\n"));
        \\        }
        \\    };
        \\}
        \\
    );
    return w.items;
}

/// A source and precision, and the cores that use it.
const Job = struct {
    kernel: Kernel,
    cores: std.ArrayList([]const u8) = .empty,
};

fn findJob(jobs: *std.ArrayList(Job), k: Kernel) ?*Job {
    for (jobs.items) |*j| {
        if (std.mem.eql(u8, j.kernel.role_name, k.role_name) and std.mem.eql(u8, j.kernel.source, k.source)) return j;
    }
    return null;
}

const Entry = struct { slot: []const u8, file: []const u8 };

fn convert(a: std.mem.Allocator, ar: Arch, role_list: []const roles.Role, out_root: []const u8, makefile: []const u8) !void {
    const io = map.io;
    const inc = try a.print(".zig-cache/blas_asm_inc_{s}", .{ar.name});
    try Io.Dir.cwd().createDirPath(io, inc);
    try Io.Dir.cwd().writeFile(io, .{ .sub_path = try map.path(&.{ inc, "config.h" }), .data = ar.config_h });

    const dir = try a.print("{s}/{s}/asm", .{ out_root, ar.name });

    // which sources and precisions the cores use
    var jobs: std.ArrayList(Job) = .empty;
    const all_cores = try map.cores(ar.name) orelse return;
    for (all_cores) |core| {
        for (role_list) |r| {
            for (try map.resolveKey(ar.name, core, r.key, makefile)) |source| {
                if (!ar.is_asm(source)) continue;
                const k: Kernel = .{ .source = source, .prec = r.prec, .conj = r.conj, .role_name = r.name, .slot = r.slot.?, .defines = r.defines };
                const job = findJob(&jobs, k) orelse blk: {
                    try jobs.append(a, .{ .kernel = k });
                    break :blk &jobs.items[jobs.items.len - 1];
                };
                try job.cores.append(a, core);
            }
        }
    }
    if (jobs.items.len == 0) return;

    // a core with more than one alternative for a kernel (chosen upstream by
    // compiler version or byte order) is not converted for it
    var ambiguous: std.StringHashMapUnmanaged(void) = .empty;
    for (all_cores) |core| {
        for (role_list) |r| {
            if ((try map.resolveKey(ar.name, core, r.key, makefile)).len < 2) continue;
            try ambiguous.put(a, try a.print("{s}:{s}", .{ try lowered(a, core), r.name }), {});
        }
    }

    var entries: std.ArrayList(Entry) = .empty;
    // the cores of each skipped slot that keeps them (`roles.Skip.cores`)
    var skips: std.array_hash_map.String(std.array_hash_map.String(void)) = .empty;
    for (jobs.items) |job| {
        const k = job.kernel;
        // sources the toolchain cannot assemble are not converted at all
        const unsupported = for (ar.unsupported) |u| {
            if (std.mem.eql(u8, u, k.source)) break true;
        } else std.mem.find(u8, k.source, "_hummer") != null;
        if (unsupported) continue;

        // nor are kernels no core would install (`roles.skip`)
        if (roles.skipped(ar.name, k.slot)) |s| {
            if (!s.cores) continue;
            const e = try skips.getOrPut(a, k.slot);
            if (!e.found_existing) e.value_ptr.* = .empty;
            for (job.cores.items) |core| {
                const lower = try lowered(a, core);
                if (ambiguous.contains(try a.print("{s}:{s}", .{ lower, try k.role(a) }))) continue;
                try e.value_ptr.put(a, lower, {});
            }
            continue;
        }

        // every core's expansion, per variant
        var texts: std.ArrayList(Text) = .empty;
        for (ar.variants) |v| for (job.cores.items) |core| {
            if (!v.builds(core)) continue;
            const lower = try lowered(a, core);
            if (ambiguous.contains(try a.print("{s}:{s}", .{ lower, try k.role(a) }))) continue;
            var lines = try arch.nakedLines(a, ar, try preprocess(a, ar, inc, k, core, v));
            if (ar.check) |check| for (lines) |l| check(k.source, l);
            if (v.rewrite) |rewrite| {
                const out = try a.alloc([]const u8, lines.len);
                for (lines, out) |l, *o| o.* = try rewrite(a, l);
                lines = out;
            }
            try texts.append(a, .{ .core = lower, .variant = v.name, .lines = lines });
        };
        if (texts.items.len == 0) continue;

        // close expansions share a file
        var clusters: std.ArrayList(std.ArrayList(usize)) = .empty;
        for (texts.items, 0..) |t, i| {
            for (clusters.items) |*c| {
                if (try similar(a, texts.items[c.items[0]].lines, t.lines)) {
                    try c.append(a, i);
                    break;
                }
            } else {
                var c: std.ArrayList(usize) = .empty;
                try c.append(a, i);
                try clusters.append(a, c);
            }
        }

        // a file is named for its variants unless it has them all, and for its
        // first core where that alone leaves two files the same name
        const variant_tags = try a.alloc([]const u8, clusters.items.len);
        for (clusters.items, variant_tags) |c, *vt| {
            var names: std.ArrayList(u8) = .empty;
            var count: usize = 0;
            for (ar.variants) |v| {
                const has = for (c.items) |i| {
                    if (std.mem.eql(u8, texts.items[i].variant, v.name)) break true;
                } else false;
                if (!has) continue;
                count += 1;
                if (names.items.len > 0) try names.append(a, '_');
                try names.appendSlice(a, v.name);
            }
            vt.* = if (count == ar.variants.len) "" else names.items;
        }

        for (clusters.items, variant_tags, 0..) |c, variant_tag, ci| {
            if (c.items.len > 128) map.die("{s}: more than 128 expansions in one file", .{k.source});
            const members = try a.alloc(Text, c.items.len);
            for (c.items, members) |i, *mbr| mbr.* = texts.items[i];
            var names: std.array_hash_map.String(void) = .empty;
            for (members) |t| for (t.lines) |l| if (definedLabel(l)) |n| if (!isOwnLabel(n)) try names.put(a, n, {});
            const bodies = try a.alloc([]const []const u8, members.len);
            for (members, bodies) |t, *b| b.* = try finish(a, t.lines, &names);
            const merged = try merge(a, bodies);

            const shared = for (variant_tags, 0..) |other, oi| {
                if (oi != ci and std.mem.eql(u8, other, variant_tag)) break true;
            } else false;
            var tag: std.ArrayList(u8) = .empty;
            if (shared) try tag.appendSlice(a, members[0].core);
            if (shared and variant_tag.len > 0) try tag.append(a, '_');
            try tag.appendSlice(a, variant_tag);
            const file = try fileName(a, k, "", if (tag.items.len > 0) tag.items else null);
            const dest = try map.path(&.{ dir, file });
            try Io.Dir.cwd().createDirPath(io, dir);
            try Io.Dir.cwd().writeFile(io, .{ .sub_path = dest, .data = try emit(a, ar, k, members, merged) });
            _ = try std.process.run(a, io, .{ .argv = &.{ "zig", "fmt", dest } });
            try entries.append(a, .{ .slot = k.slot, .file = file });
            std.debug.print("wrote {s} ({d} lines, {d} expansions)\n", .{ dest, merged.len, members.len });
        }
    }

    // which file serves which kernel; none where every kernel was skipped
    // and none keeps its cores
    if (entries.items.len == 0 and skips.count() == 0) return;
    try Io.Dir.cwd().createDirPath(io, dir);
    var idx: std.ArrayList(u8) = .empty;
    try put(&idx, a,
        \\//! Which converted file has which kernel; each says which cores it serves.
        \\//! Generated by `tools/blas_asm/gen.zig`; re-derive it rather than editing it.
        \\
        \\const registry = @import("../../../../asm_registry.zig");
        \\
        \\pub const Entry = registry.Entry;
        \\
        \\pub const entries = [_]Entry{{
        \\
    , .{});
    for (entries.items) |e| try put(&idx, a, "    .{{ .slot = .{s}, .module = @import(\"{s}\") }},\n", .{ e.slot, e.file });
    try put(&idx, a, "}};\n", .{});
    if (skips.count() > 0) {
        try put(&idx, a,
            \\
            \\/// The kernels upstream builds here that are not converted
            \\/// (`tools/blas_asm/roles.zig`'s `skip`), with the cores it builds them for.
            \\pub const skipped = [_]registry.Skipped{{
            \\
        , .{});
        for (skips.keys(), skips.values()) |slot, cores| {
            const names = try a.dupe([]const u8, cores.keys());
            std.mem.sort([]const u8, names, {}, struct {
                fn lt(_: void, x: []const u8, y: []const u8) bool {
                    return std.mem.lessThan(u8, x, y);
                }
            }.lt);
            try put(&idx, a, "    .{{ .slot = .{s}, .cores = &.{{", .{slot});
            for (names) |c| try put(&idx, a, " \"{s}\",", .{c});
            try put(&idx, a, " }} }},\n", .{});
        }
        try put(&idx, a, "}};\n", .{});
    }
    const idx_path = try map.path(&.{ dir, "__root__.zig" });
    try Io.Dir.cwd().writeFile(io, .{ .sub_path = idx_path, .data = idx.items });
    _ = try std.process.run(a, io, .{ .argv = &.{ "zig", "fmt", idx_path } });
}

pub fn main(init: std.process.Init) !void {
    const a = init.arena.allocator();
    const checkout = init.environ_map.get("OPENBLAS") orelse
        map.die("set OPENBLAS to the OpenBLAS checkout", .{});
    try map.init(a, init.io, checkout);
    const base = init.environ_map.get("BLAS_OUT") orelse "src/linalg/blas";

    var want: std.ArrayList([]const u8) = .empty;
    var args = try init.minimal.args.iterateAllocator(a);
    _ = args.next();
    while (args.next()) |arg| try want.append(a, arg);
    const wants = struct {
        fn any(list: []const []const u8, names: []const []const u8) bool {
            for (list) |w| for (names) |n| if (std.mem.eql(u8, w, n)) return true;
            return false;
        }
    }.any;

    // every rule's role, grouped into its routine's `blas<n>/<dir>/`
    const Group = struct { level: Level, dir: []const u8, roles: std.ArrayList(roles.Role) = .empty };
    var groups: std.array_hash_map.String(Group) = .empty;
    for (levels) |lv| {
        const text = map.read(try map.path(&.{ map.root, "kernel", lv.makefile })).?;
        for (try roles.parse(a, text)) |r| {
            if (r.extended or r.slot == null) continue;
            const d = try roles.dirFor(a, r.key);
            const e = try groups.getOrPut(a, try a.print("{s}/{s}", .{ lv.dir, d }));
            if (!e.found_existing) e.value_ptr.* = .{ .level = lv, .dir = d };
            try e.value_ptr.roles.append(a, r);
        }
    }
    var arch_named = false;
    var dir_named = false;
    var level_named = false;
    for (arches) |ar| arch_named = arch_named or wants(want.items, &.{ar.name});
    for (groups.values()) |g| dir_named = dir_named or wants(want.items, &.{g.dir});
    for (levels) |lv| level_named = level_named or wants(want.items, &.{lv.dir});

    for (groups.keys(), groups.values()) |path, g| {
        if (dir_named and !wants(want.items, &.{g.dir})) continue;
        if (level_named and !wants(want.items, &.{g.level.dir})) continue;
        const out_root = try a.print("{s}/{s}", .{ base, path });
        for (arches) |ar| {
            if (arch_named and !wants(want.items, &.{ar.name})) continue;
            try convert(a, ar, g.roles.items, out_root, g.level.makefile);
        }
    }
    for (arches) |ar| if (!arch_named or wants(want.items, &.{ar.name})) try checkRetired(a, ar);
    if (want.items.len == 0) for (arches) |ar| checkApplied(ar);
}

/// A BLAS level: the rules that build its kernels, and its directory.
const Level = struct { makefile: []const u8, dir: []const u8 };

const levels = [_]Level{
    .{ .makefile = "Makefile.L1", .dir = "blas1" },
    .{ .makefile = "Makefile.L2", .dir = "blas2" },
};
