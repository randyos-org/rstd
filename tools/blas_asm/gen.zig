//! Converts OpenBLAS `axpy` `.S` kernels into Zig files whose bodies are the
//! same assembly, written inline.
//!
//!     OPENBLAS=<checkout> zig run tools/blas_asm/gen.zig [-- <arch>...]
//!
//! With no architecture named, converts all that `arch/<name>.zig` describes.
//!
//! The source is preprocessed once per core, since OpenBLAS selects code with
//! `#ifdef <CORE>`. Which cores use a source, and for which precisions, comes
//! from the `KERNEL.<core>` files (`kernel_map.zig`). Cores whose expansions
//! are identical share one file, named for the first of them, with the rest
//! listed in its header; when that is every core that uses it, there is no
//! core in the name. Merging the files into fewer, idiomatic ones is done by
//! hand afterward; the header of each says where its text came from.
//!
//! The body is left as the C preprocessor emits it. Only these change, none
//! of them instructions:
//!
//!   1. the directives that declare the function are dropped; Zig declares it,
//!   2. the argument shuffling the source does by hand is dropped; Zig's
//!      calling convention delivers the arguments, into the registers the
//!      source picked, pinned by the operand list,
//!   3. each `ret` becomes a jump to the end of the block, since the block is
//!      inside a Zig function that has its own epilogue,
//!   4. tabs become spaces, which Zig's tokenizer requires inside `\\`.
//!
//! The registers stay exactly as OpenBLAS chose them. Nothing is reallocated.
//! The clobber list is read off the body, so a register the source starts
//! using is declared without anyone remembering to.
//!
//! Each function has the reference's kernel signature, so it is exported as
//! the kernel directly. `asm/__root__.zig` lists which cores use which file.

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

/// The axpy kernels keep pinned operands; every other routine is a naked
/// function of the source's own.
fn isAxpy(k: Kernel) bool {
    return std.mem.indexOf(u8, k.role_name, "axpy") != null;
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

/// The final macro definitions of the source, for `arch.macros`.
fn readMacros(a: std.mem.Allocator, ar: Arch, inc: []const u8, k: Kernel, core: []const u8, v: arch.Variant) !void {
    arch.macros.clearRetainingCapacity();
    var it = std.mem.splitScalar(u8, try runCpp(a, ar, inc, k, core, v, true), '\n');
    while (it.next()) |raw| {
        const line = std.mem.trim(u8, raw, " \r");
        const prefix = "#define ";
        if (!std.mem.startsWith(u8, line, prefix)) continue;
        const rest = line[prefix.len..];
        const sp = std.mem.indexOfScalar(u8, rest, ' ') orelse continue;
        try arch.macros.put(a, rest[0..sp], std.mem.trim(u8, rest[sp..], " "));
    }
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

fn same(x: []const []const u8, y: []const []const u8) bool {
    if (x.len != y.len) return false;
    for (x, y) |p, q| if (!std.mem.eql(u8, p, q)) return false;
    return true;
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
/// (GAS's per-macro-invocation unique-label marker) if there is one: the
/// assembler resolves that itself, per invocation, so it is left alone and
/// just carried along as part of the name for suffixing purposes.
fn nameEnd(s: []const u8, i: usize) usize {
    var e = i;
    while (e < s.len and isName(s[e])) e += 1;
    if (e + 1 < s.len and s[e] == '\\' and s[e + 1] == '@') e += 2;
    return e;
}

/// Every label and macro the assembler sees is global to the object, and one
/// object holds several of these kernels, so each name gets the file's name
/// appended and each macro is removed again at the end. Numeric labels are
/// local to their block and stay as they are.
fn finish(a: std.mem.Allocator, lines: []const []const u8, suffix: []const u8) ![]const []const u8 {
    var labels: std.StringArrayHashMapUnmanaged(void) = .empty;
    var macros: std.ArrayList([]const u8) = .empty;
    for (lines) |line| {
        const t = std.mem.trim(u8, line, " ");
        if (std.mem.startsWith(u8, t, ".macro ")) {
            var it = std.mem.tokenizeAny(u8, t[".macro ".len..], " ,");
            if (it.next()) |name| try macros.append(a, name);
            continue;
        }
        const e = nameEnd(t, 0);
        if (e > 0 and e < t.len and t[e] == ':' and !std.ascii.isDigit(t[0])) try labels.put(a, t[0..e], {});
    }

    var out: std.ArrayList([]const u8) = .empty;
    for (lines) |line| {
        // the function's own label; Zig declares the function
        const t = std.mem.trim(u8, line, " ");
        if (std.mem.eql(u8, t, "k:") or std.mem.eql(u8, t, "k_:")) continue;
        var buf: std.ArrayList(u8) = .empty;
        var i: usize = 0;
        while (i < line.len) {
            if (isName(line[i]) and (i == 0 or !isName(line[i - 1]))) {
                const e = nameEnd(line, i);
                const tok = line[i..e];
                try buf.appendSlice(a, tok);
                if (labels.contains(tok)) try buf.appendSlice(a, suffix);
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

/// The text of one Zig file.
fn emit(a: std.mem.Allocator, ar: Arch, k: Kernel, cores: []const u8, lines: []const []const u8) ![]const u8 {
    var w: std.ArrayList(u8) = .empty;
    const F = k.prec.float();
    const role = try k.role(a);
    if (ar.naked or !isAxpy(k)) {
        // the whole function is the source's: it reads its arguments from the
        // ABI's places itself, saves what it uses, and returns
        try put(&w, a,
            \\//! `{s}`, from OpenBLAS `kernel/{s}/{s}`, as expanded for {s}.
            \\//!
            \\//! Generated by `tools/blas_asm/gen.zig`; the function is the assembly of
            \\//! that source, unchanged, as a naked function, so it reads its own
            \\//! arguments as the reference's does. Re-derive it rather than editing it.
            \\
            \\pub fn {s}() callconv(.naked) void {{
            \\    asm volatile (
            \\
        , .{ role, ar.name, k.source, cores, role });
        for (lines) |l| try put(&w, a, "        \\\\{s}\n", .{l});
        try put(&w, a, "    );\n}}\n", .{});
        return w.items;
    }
    const pins = try ar.pins(a, k);

    try put(&w, a,
        \\//! `{s}`, from OpenBLAS `kernel/{s}/{s}`, as expanded for {s}.
        \\//!
        \\//! Generated by `tools/blas_asm/gen.zig`; the body is the assembly of
        \\//! that source, unchanged, with OpenBLAS's own registers pinned by the
        \\//! operand list. Re-derive it rather than editing it.
        \\
        \\pub noinline fn {s}(
        \\    n: isize,
        \\    _: isize,
        \\    _: isize,
        \\
    , .{ role, ar.name, k.source, cores, role });
    if (k.prec.complex()) {
        try put(&w, a, "    alpha_r: {s},\n    alpha_i: {s},\n", .{ F, F });
    } else {
        try put(&w, a, "    alpha: {s},\n", .{F});
    }
    try put(&w, a, "    x: [*]{s},\n    inc_x: isize,\n    y: [*]{s},\n    inc_y: isize,\n    _: [*]{s},\n    _: isize,\n) callconv(.c) i32 {{\n", .{ F, F, F });

    for (pins) |pin| if (pin.init) |i| try put(&w, a, "    {s} {s} = {s};\n", .{ if (pin.inout) "var" else "const", pin.variable, i });
    for (try ar.decls(a, k)) |d| try put(&w, a, "    {s}\n", .{d});
    try put(&w, a, "    asm volatile (\n", .{});
    for (lines) |l| try put(&w, a, "        \\\\{s}\n", .{l});

    // outputs, then inputs
    var first = true;
    for (pins) |pin| {
        if (!pin.inout) continue;
        try put(&w, a, "        {s} [{s}] \"+{{{s}}}\" ({s}),\n", .{ if (first) ":" else " ", pin.name, pin.reg, pin.variable });
        first = false;
    }
    var any_in = false;
    for (pins) |pin| {
        if (pin.inout) continue;
        try put(&w, a, "        {s} [{s}] \"{{{s}}}\" ({s}),\n", .{ if (!any_in) ":" else " ", pin.name, pin.reg, pin.variable });
        any_in = true;
    }
    if (!any_in) try put(&w, a, "        :\n", .{});
    try put(&w, a, "        : .{{", .{});
    var sep: []const u8 = " ";
    for (ar.always_clobbers) |c| {
        try put(&w, a, "{s}.{s} = true", .{ sep, c });
        sep = ", ";
    }
    for (try ar.clobbers(a, lines, pins)) |c| try put(&w, a, ", .{s} = true", .{c});
    try put(&w, a, " }});\n    return 0;\n}}\n", .{});
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

const Entry = struct { role: []const u8, variant: []const u8, cores: []const u8, file: []const u8 };

fn convert(a: std.mem.Allocator, ar: Arch, role_list: []const roles.Role, out_root: []const u8) !void {
    const io = map.io;
    const inc = try std.fmt.allocPrint(a, ".zig-cache/blas_asm_inc_{s}", .{ar.name});
    try Io.Dir.cwd().createDirPath(io, inc);
    try Io.Dir.cwd().writeFile(io, .{ .sub_path = try map.path(&.{ inc, "config.h" }), .data = ar.config_h });

    const dir = try std.fmt.allocPrint(a, "{s}/{s}/asm", .{ out_root, ar.name });
    try Io.Dir.cwd().createDirPath(io, dir);

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

    // a core with more than one alternative for a kernel (chosen upstream by
    // compiler version or byte order) is not indexed for it
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
        for (ar.variants) |v| {
            const Group = struct { text: []const []const u8, names: std.ArrayList(u8) = .empty, cores: std.ArrayList(u8) = .empty, first: []const u8 };
            var groups: std.ArrayList(Group) = .empty;
            try readMacros(a, ar, inc, k, job.cores.items[0], v);
            outer: for (job.cores.items) |core| {
                const lines = if (isAxpy(k)) try ar.transform(a, k, try preprocess(a, ar, inc, k, core, v)) else try arch.nakedLines(a, try preprocess(a, ar, inc, k, core, v));
                const lower = try lowered(a, core);
                const key = try std.fmt.allocPrint(a, "{s}:{s}", .{ lower, try k.role(a) });
                const indexed = !ambiguous.contains(key);
                for (groups.items) |*g| {
                    if (same(g.text, lines)) {
                        try g.names.print(a, ", {s}", .{lower});
                        if (indexed) try g.cores.print(a, ",{s}", .{lower});
                        continue :outer;
                    }
                }
                var g: Group = .{ .text = lines, .first = lower };
                try g.names.print(a, "{s}", .{lower});
                if (indexed) try g.cores.print(a, "{s}", .{lower});
                try groups.append(a, g);
            }

            for (groups.items) |g| {
                const cores = g.names.items;
                const width = if (v.name.len > 0) try std.fmt.allocPrint(a, "_{s}", .{v.name}) else "";
                const file = try fileName(a, k, width, if (groups.items.len == 1) null else g.first);
                // sources the toolchain cannot assemble are not converted at all
                const unsupported = for (ar.unsupported) |u| {
                    if (std.mem.eql(u8, u, k.source)) break true;
                } else std.mem.indexOf(u8, k.source, "_hummer") != null;
                if (unsupported) continue;
                const dest = try map.path(&.{ dir, file });
                const suffix = try std.fmt.allocPrint(a, "_{s}", .{file[0 .. file.len - ".zig".len]});
                try Io.Dir.cwd().writeFile(io, .{ .sub_path = dest, .data = try emit(a, ar, k, cores, try finish(a, g.text, suffix)) });
                _ = try std.process.run(a, io, .{ .argv = &.{ "zig", "fmt", dest } });
                if (g.cores.items.len > 0) try entries.append(a, .{ .role = try k.role(a), .variant = v.name, .cores = g.cores.items, .file = file });
                std.debug.print("wrote {s} ({d} body lines): {s}\n", .{ dest, g.text.len, cores });
            }
        }
    }

    // which file serves which core
    var idx: std.ArrayList(u8) = .empty;
    try put(&idx, a,
        \\//! Which converted file serves which cores, per kernel. Generated by
        \\//! `tools/blas_asm/gen.zig`; re-derive it rather than editing it.
        \\
        \\const std = @import("std");
        \\
        \\pub const Entry = struct {{
        \\    role: []const u8,
        \\    /// Empty when the file is the same however the source is built.
        \\    variant: []const u8,
        \\    cores: []const []const u8,
        \\    module: type,
        \\}};
        \\
        \\pub const entries = [_]Entry{{
        \\
    , .{});
    for (entries.items) |e| {
        try put(&idx, a, "    .{{ .role = \"{s}\", .variant = \"{s}\", .cores = &.{{", .{ e.role, e.variant });
        var it = std.mem.splitScalar(u8, e.cores, ',');
        while (it.next()) |c| try put(&idx, a, " \"{s}\",", .{c});
        try put(&idx, a, " }}, .module = @import(\"{s}\") }},\n", .{e.file});
    }
    try put(&idx, a, "}};\n\n", .{});
    try idx.appendSlice(a,
        \\/// The file whose function `role` is, for `core` built as `variant`, if it
        \\/// is assembly there.
        \\pub fn source(comptime role: []const u8, comptime core: []const u8, comptime variant: []const u8) ?type {
        \\    @setEvalBranchQuota(1_000_000);
        \\    inline for (entries) |e| {
        \\        if (comptime std.mem.eql(u8, e.role, role) and (e.variant.len == 0 or std.mem.eql(u8, e.variant, variant))) {
        \\            inline for (e.cores) |c| {
        \\                if (comptime std.mem.eql(u8, c, core)) return e.module;
        \\            }
        \\        }
        \\    }
        \\    return null;
        \\}
        \\
    );
    const idx_path = try map.path(&.{ dir, "__root__.zig" });
    try Io.Dir.cwd().writeFile(io, .{ .sub_path = idx_path, .data = idx.items });
    _ = try std.process.run(a, io, .{ .argv = &.{ "zig", "fmt", idx_path } });
}

pub fn main(init: std.process.Init) !void {
    const a = init.arena.allocator();
    const checkout = init.environ_map.get("OPENBLAS") orelse
        map.die("set OPENBLAS to the OpenBLAS checkout", .{});
    try map.init(a, init.io, checkout);
    if (init.environ_map.get("BLAS_OUT")) |o| map.out_root = o;

    var want: std.ArrayList([]const u8) = .empty;
    var others = false;
    var args = try init.minimal.args.iterateAllocator(a);
    _ = args.next();
    while (args.next()) |arg| {
        if (std.mem.eql(u8, arg, "--others")) others = true else try want.append(a, arg);
    }

    // the roles: axpy, kept as it was, or every other routine of the level-1
    // rules, each going to its own `blas1/<dir>/`, as axpy already does
    const text = map.read(try map.path(&.{ map.root, "kernel", "Makefile.L1" })).?;
    var by_dir: std.StringArrayHashMapUnmanaged(std.ArrayList(roles.Role)) = .empty;
    for (try roles.parse(a, text)) |r| {
        const is_axpy = std.mem.indexOf(u8, r.name, "axpy") != null and std.mem.indexOf(u8, r.name, "axpby") == null;
        if (is_axpy == others) continue;
        const dir = if (others) try roles.dirFor(a, r.key) else "axpy";
        const e = try by_dir.getOrPut(a, dir);
        if (!e.found_existing) e.value_ptr.* = .empty;
        try e.value_ptr.append(a, r);
    }
    const base = init.environ_map.get("BLAS_OUT_KERNELS") orelse "src/linalg/blas/blas1";

    var dirs = by_dir.iterator();
    while (dirs.next()) |kv| {
        // axpy keeps `BLAS_OUT` (its whole out_root, as before); every other
        // routine nests under `base/<dir>`
        const out_root = if (others) try std.fmt.allocPrint(a, "{s}/{s}", .{ base, kv.key_ptr.* }) else map.out_root;
        for (arches) |ar| {
            if (want.items.len > 0) {
                const listed = for (want.items) |w| {
                    if (std.mem.eql(u8, w, ar.name)) break true;
                } else false;
                if (!listed) continue;
            }
            try convert(a, ar, kv.value_ptr.items, out_root);
        }
    }
}
