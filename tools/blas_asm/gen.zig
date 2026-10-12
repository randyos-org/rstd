//! Converts OpenBLAS's level-1 `.S` kernels into Zig files whose bodies are
//! the same assembly, as naked functions.
//!
//!     OPENBLAS=<checkout> zig run tools/blas_asm/gen.zig [-- <arch|routine>...]
//!
//! With nothing named, converts every routine for every architecture that
//! `arch/<name>.zig` describes; otherwise only the architectures and routines
//! (`blas1/<routine>/` names, such as `axpy`) named.
//!
//! The source is preprocessed once per core and variant, since OpenBLAS
//! selects code with `#ifdef <CORE>` (and the ABI with `WINDOWS_ABI` and the
//! like). Which cores use a source, and for which precisions, comes from the
//! `KERNEL.<core>` files (`kernel_map.zig`). Cores whose expansions are
//! identical share one file.
//!
//! Each function is the source's own, unchanged: it reads its arguments where
//! the calling convention puts them, saves what it uses, and returns. Only the
//! directives that declare the function go (Zig declares it), tabs become
//! spaces (Zig's tokenizer requires it inside `\`), and labels get the
//! file's name appended (one object holds many of these).
//!
//! `asm/__root__.zig` lists which cores use which file.

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
    const rest = if (std.mem.indexOfScalar(u8, stem, '_')) |at| stem[at..] else "";
    const core = if (first) |f| try std.fmt.allocPrint(a, "_{s}", .{f}) else "";
    return std.fmt.allocPrint(a, "{s}{s}{s}{s}.zig", .{ try k.role(a), rest, width, core });
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
    try argv.appendSlice(a, &.{ "zig", "cc", "-E", "-P", "-x", "assembler-with-cpp", "-w", "-DDYNAMIC_ARCH", "-I", inc, "-I", map.root, "-I", try map.path(&.{ map.root, "kernel", ar.name }), "-target", ar.triple });
    try argv.appendSlice(a, ar.defines);
    if (k.prec.double()) try argv.append(a, "-DDOUBLE");
    if (k.prec.complex()) try argv.append(a, "-DCOMPLEX");
    if (k.conj) try argv.append(a, "-DCONJ");
    try argv.appendSlice(a, k.defines);
    try argv.appendSlice(a, try ar.core_defines(a, k, core));
    try argv.appendSlice(a, v.defines);
    if (dump_macros) try argv.append(a, "-dM");
    try argv.append(a, try map.path(&.{ map.root, "kernel", ar.name, k.source }));
    try argv.appendSlice(a, &.{ "-o", "-" });

    const r = std.process.run(a, map.io, .{ .argv = argv.items }) catch |e|
        map.die("cannot run the preprocessor: {t}", .{e});
    if (r.stderr.len > 0) std.debug.print("{s}", .{r.stderr});
    return r.stdout;
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
    while (e < s.len and isName(s[e])) e += 1;
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
fn finish(a: std.mem.Allocator, lines: []const []const u8, names: *const std.StringArrayHashMapUnmanaged(void)) ![]const []const u8 {
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
                        if (!found) std.debug.panic("label {s} is referenced across a macro's boundary: {s}", .{ tok, line });
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
    for (macros.items) |m| try out.append(a, try std.fmt.allocPrint(a, " .purgem {s}", .{m}));
    return out.items;
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
    var conds: std.AutoArrayHashMapUnmanaged(u128, void) = .empty;
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
        try put(&w, a, "/// Which parts each key's function has.\nconst parts = [_][{d}]bool{{\n", .{part_count});
        for (texts, 0..) |_, i| {
            try put(&w, a, "    .{{", .{});
            for (conds.keys(), 0..) |c, j| try put(&w, a, "{s}{s}", .{ if (j == 0) " " else ", ", if (c & (@as(u128, 1) << @intCast(i)) != 0) "true" else "false" });
            try put(&w, a, " }},\n", .{});
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
    if (part_count > 0) {
        try put(&w, a, "fn Body(comptime on: [{d}]bool) type {{\n    return struct {{\n", .{part_count});
    } else {
        try put(&w, a, "const Body = struct {{\n", .{});
    }
    try put(&w, a, "pub fn {s}() callconv(.naked) void {{\nasm volatile (\n", .{role});
    var i: usize = 0;
    var first = true;
    while (i < lines.len) {
        const mask = lines[i].mask;
        var e = i;
        while (e < lines.len and lines[e].mask == mask) e += 1;
        if (!first) try put(&w, a, "++ ", .{});
        if (mask != all) try put(&w, a, "(if (on[{d}])\n", .{conds.getIndex(mask).?});
        for (lines[i..e]) |l| try put(&w, a, "\\\\{s}\n", .{l.text});
        try put(&w, a, "\\\\\n", .{});
        if (mask != all) try put(&w, a, "else \"\")\n", .{});
        first = false;
        i = e;
    }
    try w.appendSlice(a, ");\n}\n");
    try w.appendSlice(a, if (part_count > 0) "};\n}\n" else "};\n");
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

const Entry = struct { role: []const u8, file: []const u8 };

fn convert(a: std.mem.Allocator, ar: Arch, role_list: []const roles.Role, out_root: []const u8) !void {
    const io = map.io;
    const inc = try std.fmt.allocPrint(a, ".zig-cache/blas_asm_inc_{s}", .{ar.name});
    try Io.Dir.cwd().createDirPath(io, inc);
    try Io.Dir.cwd().writeFile(io, .{ .sub_path = try map.path(&.{ inc, "config.h" }), .data = ar.config_h });

    const dir = try std.fmt.allocPrint(a, "{s}/{s}/asm", .{ out_root, ar.name });

    // which sources and precisions the cores use
    var jobs: std.ArrayList(Job) = .empty;
    const all_cores = try map.cores(ar.name) orelse return;
    for (all_cores) |core| {
        for (role_list) |r| {
            for (try map.resolveKey(ar.name, core, r.key, "Makefile.L1")) |source| {
                if (!ar.is_asm(source)) continue;
                const k: Kernel = .{ .source = source, .prec = r.prec, .conj = r.conj, .role_name = r.name, .defines = r.defines };
                const job = findJob(&jobs, k) orelse blk: {
                    try jobs.append(a, .{ .kernel = k });
                    break :blk &jobs.items[jobs.items.len - 1];
                };
                try job.cores.append(a, core);
            }
        }
    }
    if (jobs.items.len == 0) return;
    try Io.Dir.cwd().createDirPath(io, dir);

    // a core with more than one alternative for a kernel (chosen upstream by
    // compiler version or byte order) is not converted for it
    var ambiguous: std.StringHashMapUnmanaged(void) = .empty;
    for (all_cores) |core| {
        for (role_list) |r| {
            if ((try map.resolveKey(ar.name, core, r.key, "Makefile.L1")).len < 2) continue;
            try ambiguous.put(a, try std.fmt.allocPrint(a, "{s}:{s}", .{ try lowered(a, core), r.name }), {});
        }
    }

    var entries: std.ArrayList(Entry) = .empty;
    for (jobs.items) |job| {
        const k = job.kernel;
        // sources the toolchain cannot assemble are not converted at all
        const unsupported = for (ar.unsupported) |u| {
            if (std.mem.eql(u8, u, k.source)) break true;
        } else std.mem.indexOf(u8, k.source, "_hummer") != null;
        if (unsupported) continue;

        // every core's expansion, per variant
        var texts: std.ArrayList(Text) = .empty;
        for (ar.variants) |v| for (job.cores.items) |core| {
            const lower = try lowered(a, core);
            if (ambiguous.contains(try std.fmt.allocPrint(a, "{s}:{s}", .{ lower, try k.role(a) }))) continue;
            try texts.append(a, .{ .core = lower, .variant = v.name, .lines = try arch.nakedLines(a, try preprocess(a, ar, inc, k, core, v)) });
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

        // a file is named for its variant where there are several, and for its
        // first core where that alone leaves two files the same name
        const variant_tags = try a.alloc([]const u8, clusters.items.len);
        for (clusters.items, variant_tags) |c, *vt| {
            const v = texts.items[c.items[0]].variant;
            const one = for (c.items) |i| {
                if (!std.mem.eql(u8, texts.items[i].variant, v)) break false;
            } else true;
            vt.* = if (one and ar.variants.len > 1) v else "";
        }

        for (clusters.items, variant_tags, 0..) |c, variant_tag, ci| {
            if (c.items.len > 128) map.die("{s}: more than 128 expansions in one file", .{k.source});
            const members = try a.alloc(Text, c.items.len);
            for (c.items, members) |i, *mbr| mbr.* = texts.items[i];
            var names: std.StringArrayHashMapUnmanaged(void) = .empty;
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
            try Io.Dir.cwd().writeFile(io, .{ .sub_path = dest, .data = try emit(a, ar, k, members, merged) });
            _ = try std.process.run(a, io, .{ .argv = &.{ "zig", "fmt", dest } });
            try entries.append(a, .{ .role = try k.role(a), .file = file });
            std.debug.print("wrote {s} ({d} lines, {d} expansions)\n", .{ dest, merged.len, members.len });
        }
    }

    // which file serves which kernel
    var idx: std.ArrayList(u8) = .empty;
    try put(&idx, a,
        \\//! Which converted file has which kernel; each says which cores it serves.
        \\//! Generated by `tools/blas_asm/gen.zig`; re-derive it rather than editing it.
        \\
        \\pub const Entry = struct {{ role: []const u8, module: type }};
        \\
        \\pub const entries = [_]Entry{{
        \\
    , .{});
    for (entries.items) |e| try put(&idx, a, "    .{{ .role = \"{s}\", .module = @import(\"{s}\") }},\n", .{ e.role, e.file });
    try put(&idx, a, "}};\n", .{});
    const idx_path = try map.path(&.{ dir, "__root__.zig" });
    try Io.Dir.cwd().writeFile(io, .{ .sub_path = idx_path, .data = idx.items });
    _ = try std.process.run(a, io, .{ .argv = &.{ "zig", "fmt", idx_path } });
}

pub fn main(init: std.process.Init) !void {
    const a = init.arena.allocator();
    const checkout = init.environ_map.get("OPENBLAS") orelse
        map.die("set OPENBLAS to the OpenBLAS checkout", .{});
    try map.init(a, init.io, checkout);
    const base = init.environ_map.get("BLAS_OUT") orelse "src/linalg/blas/blas1";

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

    // every level-1 rule's role, grouped into its routine's `blas1/<dir>/`
    const text = map.read(try map.path(&.{ map.root, "kernel", "Makefile.L1" })).?;
    var by_dir: std.StringArrayHashMapUnmanaged(std.ArrayList(roles.Role)) = .empty;
    for (try roles.parse(a, text)) |r| {
        const e = try by_dir.getOrPut(a, try roles.dirFor(a, r.key));
        if (!e.found_existing) e.value_ptr.* = .empty;
        try e.value_ptr.append(a, r);
    }
    var arch_named = false;
    var dir_named = false;
    for (arches) |ar| arch_named = arch_named or wants(want.items, &.{ar.name});
    for (by_dir.keys()) |d| dir_named = dir_named or wants(want.items, &.{d});

    var dirs = by_dir.iterator();
    while (dirs.next()) |kv| {
        if (dir_named and !wants(want.items, &.{kv.key_ptr.*})) continue;
        const out_root = try std.fmt.allocPrint(a, "{s}/{s}", .{ base, kv.key_ptr.* });
        for (arches) |ar| {
            if (arch_named and !wants(want.items, &.{ar.name})) continue;
            try convert(a, ar, kv.value_ptr.items, out_root);
        }
    }
}
