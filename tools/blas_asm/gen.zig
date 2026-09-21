//! Converts OpenBLAS x86_64 `axpy` `.S` kernels into Zig files whose bodies
//! are the same assembly, written inline.
//!
//!     OPENBLAS=<checkout> zig run tools/blas_asm/gen.zig
//!
//! The source is preprocessed once per core, since OpenBLAS selects code with
//! `#ifdef <CORE>`. Which cores use a source comes from the `KERNEL.<core>`
//! files (`kernel_map.zig`). Cores whose expansions are identical share one
//! file, named for the first of them, with the rest listed in its header.
//! Merging the files into fewer, idiomatic ones is done by hand afterward;
//! the header of each says where its text came from.
//!
//! The body is left as the C preprocessor emits it. Only five things change,
//! none of them instructions:
//!
//!   1. the assembler directives that declare the function are dropped; Zig
//!      declares it,
//!   2. the argument shuffling the source does by hand is dropped; Zig's
//!      calling convention delivers the arguments into the same registers
//!      the source picked, pinned by the operand list,
//!   3. `alpha` is loaded from a pointer, as the source's own Windows path
//!      loads it from the stack,
//!   4. each `ret` becomes a jump to the end of the block, since the block
//!      is inside a Zig function that has its own epilogue,
//!   5. tabs become spaces, which Zig's tokenizer requires inside `\\`.
//!
//! The registers stay exactly as OpenBLAS chose them. Nothing is reallocated.
//! The clobber list is read off the body, so a register the source starts
//! using is declared without anyone remembering to.

const std = @import("std");
const Io = std.Io;
const map = @import("kernel_map.zig");

const arch = "x86_64";

/// A core, the macro OpenBLAS's headers test for it, and whether its build
/// defines `HAVE_SSE3`.
const Core = struct {
    name: []const u8,
    macro: []const u8,
    sse3: bool,
};

/// Real or complex, single or double.
const Prec = enum {
    s,
    d,
    c,
    z,

    fn complex(p: Prec) bool {
        return p == .c or p == .z;
    }
    fn double(p: Prec) bool {
        return p == .d or p == .z;
    }
    fn float(p: Prec) []const u8 {
        return if (p.double()) "f64" else "f32";
    }
};

/// One source, and the precision it is built for.
const Kernel = struct {
    /// Under `kernel/x86_64/` in the OpenBLAS checkout.
    source: []const u8,
    prec: Prec,
    /// The conjugated complex form.
    conj: bool = false,

    fn stem(k: Kernel) []const u8 {
        return k.source[0 .. k.source.len - 2];
    }
    /// `saxpy`, `zaxpyc`, ...
    fn name(k: Kernel) []const u8 {
        return std.fmt.allocPrint(map.arena, "{t}axpy{s}", .{ k.prec, if (k.conj) "c" else "" }) catch @panic("OOM");
    }
    /// The precision letter `Makefile.L1` keys the source by.
    fn letter(k: Kernel) []const u8 {
        return switch (k.prec) {
            .s => "S",
            .d => "D",
            .c => "C",
            .z => "Z",
        };
    }
};

const kernels = [_]Kernel{
    .{ .source = "axpy_sse.S", .prec = .s },
    .{ .source = "axpy_sse2.S", .prec = .d },
    .{ .source = "axpy_atom.S", .prec = .d },
    .{ .source = "zaxpy_sse.S", .prec = .c },
    .{ .source = "zaxpy_sse.S", .prec = .c, .conj = true },
    .{ .source = "zaxpy_sse2.S", .prec = .z },
    .{ .source = "zaxpy_sse2.S", .prec = .z, .conj = true },
    .{ .source = "zaxpy_atom.S", .prec = .z },
    .{ .source = "zaxpy_atom.S", .prec = .z, .conj = true },
};

/// The registers the body expects its arguments in, and the Zig variables
/// that carry them. `rdx` holds the address of `alpha`.
const operands = [_][3][]const u8{
    .{ "m", "rdi", "m" },
    .{ "x", "rcx", "xp" },
    .{ "incx", "r8", "ix" },
    .{ "y", "r9", "yp" },
    .{ "incy", "rsi", "iy" },
    .{ "ap", "rdx", "ap" },
};

/// General registers by every name they go by, to the 64-bit one.
fn canonical(reg: []const u8) ?[]const u8 {
    const legacy = [_][]const u8{ "ax", "bx", "cx", "dx", "si", "di", "bp" };
    for (legacy) |l| {
        const wide = std.fmt.allocPrint(map.arena, "r{s}", .{l}) catch @panic("OOM");
        if (std.mem.eql(u8, reg, wide) or std.mem.eql(u8, reg[1..], l) and reg[0] == 'e') return wide;
        if (std.mem.eql(u8, reg, l)) return wide;
    }
    // al, bl, cl, dl and sil, dil, bpl
    const low = [_][2][]const u8{
        .{ "al", "rax" },  .{ "ah", "rax" },  .{ "bl", "rbx" },  .{ "bh", "rbx" },
        .{ "cl", "rcx" },  .{ "ch", "rcx" },  .{ "dl", "rdx" },  .{ "dh", "rdx" },
        .{ "sil", "rsi" }, .{ "dil", "rdi" }, .{ "bpl", "rbp" },
    };
    for (low) |l| if (std.mem.eql(u8, reg, l[0])) return l[1];
    if (reg.len >= 2 and reg[0] == 'r') {
        var end = reg.len;
        if (std.mem.indexOfScalar(u8, "dwb", reg[end - 1]) != null and end > 2 and std.ascii.isDigit(reg[end - 2])) end -= 1;
        const n = std.fmt.parseInt(u8, reg[1..end], 10) catch return null;
        if (n >= 8 and n <= 15) return std.fmt.allocPrint(map.arena, "r{d}", .{n}) catch @panic("OOM");
    }
    return null;
}

/// The registers `text` mentions that are not operands, as Zig clobber names.
fn clobbers(text: []const []const u8) []const []const u8 {
    var seen: std.StringArrayHashMapUnmanaged(void) = .empty;
    for (text) |line| {
        var i: usize = 0;
        while (std.mem.indexOfScalarPos(u8, line, i, '%')) |p| {
            var e = p + 1;
            while (e < line.len and std.ascii.isAlphanumeric(line[e])) e += 1;
            const reg = line[p + 1 .. e];
            i = e;
            if (reg.len == 0) continue;
            const is_vec = std.mem.startsWith(u8, reg, "xmm") or std.mem.startsWith(u8, reg, "ymm") or std.mem.startsWith(u8, reg, "zmm");
            const name = if (is_vec) reg else canonical(reg) orelse continue;
            if (std.mem.eql(u8, name, "rsp")) continue;
            const pinned = for (operands) |o| {
                if (std.mem.eql(u8, o[1], name)) break true;
            } else false;
            if (!pinned) seen.put(map.arena, name, {}) catch @panic("OOM");
        }
    }
    const names = seen.keys();
    std.mem.sort([]const u8, names, {}, struct {
        fn lt(_: void, a: []const u8, b: []const u8) bool {
            return std.mem.lessThan(u8, a, b);
        }
    }.lt);
    return names;
}

/// Runs the preprocessor over the source the way OpenBLAS's build does.
fn preprocess(inc: []const u8, k: Kernel, c: Core) []const u8 {
    const a = map.arena;
    var argv: std.ArrayList([]const u8) = .empty;
    argv.appendSlice(a, &.{
        "zig",                                    "cc",
        "-E",                                     "-P",
        "-x",                                     "assembler-with-cpp",
        "-w",                                     "-I",
        inc,                                      "-I",
        map.root,                                 "-I",
        map.path(&.{ map.root, "kernel", arch }), "-DASSEMBLER",
        "-DARCH_X86_64",                          "-DOS_LINUX",
        "-DASMNAME=k",                            "-DASMFNAME=k_",
        "-DNAME=k_",                              "-DCNAME=k",
        "-target",                                "x86_64-linux-gnu",
        "-DHAVE_SSE",                             "-DHAVE_SSE2",
    }) catch @panic("OOM");
    argv.append(a, std.fmt.allocPrint(a, "-D{s}", .{c.macro}) catch @panic("OOM")) catch @panic("OOM");
    argv.append(a, std.fmt.allocPrint(a, "-DCORE_{s}", .{c.macro}) catch @panic("OOM")) catch @panic("OOM");
    if (c.sse3) argv.append(a, "-DHAVE_SSE3") catch @panic("OOM");
    if (k.prec.double()) argv.append(a, "-DDOUBLE") catch @panic("OOM");
    if (k.prec.complex()) argv.append(a, "-DCOMPLEX") catch @panic("OOM");
    if (k.conj) argv.append(a, "-DCONJ") catch @panic("OOM");
    argv.append(a, map.path(&.{ map.root, "kernel", arch, k.source })) catch @panic("OOM");
    argv.appendSlice(a, &.{ "-o", "-" }) catch @panic("OOM");

    const r = std.process.run(a, map.io, .{ .argv = argv.items }) catch |e|
        map.die("cannot run the preprocessor: {t}", .{e});
    if (r.stderr.len > 0) std.debug.print("{s}", .{r.stderr});
    return r.stdout;
}

fn has(line: []const u8, needle: []const u8) bool {
    return std.mem.indexOf(u8, line, needle) != null;
}

/// Tabs become spaces to the next multiple of eight.
fn expandTabs(s: []const u8) []const u8 {
    var out: std.ArrayList(u8) = .empty;
    var col: usize = 0;
    for (s) |c| {
        if (c == '\t') {
            const n = 8 - col % 8;
            out.appendNTimes(map.arena, ' ', n) catch @panic("OOM");
            col += n;
        } else {
            out.append(map.arena, c) catch @panic("OOM");
            col += 1;
        }
    }
    return out.items;
}

/// The body, as the lines that go inside the `asm` block.
fn body(k: Kernel, expanded: []const u8) []const []const u8 {
    var out: std.ArrayList([]const u8) = .empty;
    // `alpha` arrives by pointer, as the source's Windows path also does; the
    // real part goes where the source's Linux path expects it, `xmm0`, and
    // the imaginary part in `xmm1`
    const load = if (k.prec.double()) " movsd" else " movss";
    const step: usize = if (k.prec.double()) 8 else 4;
    out.append(map.arena, std.fmt.allocPrint(map.arena, "{s} (%rdx), %xmm0", .{load}) catch @panic("OOM")) catch @panic("OOM");
    if (k.prec.complex()) out.append(map.arena, std.fmt.allocPrint(map.arena, "{s} {d}(%rdx), %xmm1", .{ load, step }) catch @panic("OOM")) catch @panic("OOM");
    var it = std.mem.splitScalar(u8, expanded, '\n');
    while (it.next()) |raw| {
        const line = std.mem.trimEnd(u8, expandTabs(raw), " \r");
        if (std.mem.trim(u8, line, " ").len == 0) continue;
        // the directives that declare the function, and the ABI shuffling
        if (has(line, ".globl") or has(line, ".size k") or has(line, ".section .note")) continue;
        if (has(line, "movq 8(%rsp)")) continue;
        if (has(line, "(%rsp)")) map.die("the body reads the stack: {s}", .{line});
        // `ret` leaves the Zig function without its epilogue; jump to the end
        if (std.mem.eql(u8, std.mem.trim(u8, line, " "), "ret")) {
            out.append(map.arena, " jmp 9f") catch @panic("OOM");
            continue;
        }
        out.append(map.arena, line) catch @panic("OOM");
    }
    out.append(map.arena, "9:") catch @panic("OOM");
    return out.items;
}

fn emit(w: *std.ArrayList(u8), k: Kernel, cores: []const u8, lines: []const []const u8) void {
    const a = map.arena;
    const F = k.prec.float();
    w.print(a,
        \\//! `{s}`, from OpenBLAS `kernel/x86_64/{s}`, as expanded for {s}.
        \\//!
        \\//! Generated by `tools/blas_asm/gen.zig`; the body is the assembly of
        \\//! that source, unchanged, with OpenBLAS's own registers pinned by the
        \\//! operand list. Re-derive it rather than editing it.
        \\
        \\pub noinline fn {s}(
    , .{ k.name(), k.source, cores, k.name() }) catch @panic("OOM");
    if (k.prec.complex()) {
        w.print(a, "n: isize, alpha_r: {s}, alpha_i: {s}, x: [*]{s}, inc_x: isize, y: [*]{s}, inc_y: isize", .{ F, F, F, F }) catch @panic("OOM");
    } else {
        w.print(a, "n: isize, alpha: {s}, x: [*]{s}, inc_x: isize, y: [*]{s}, inc_y: isize", .{ F, F, F }) catch @panic("OOM");
    }
    w.print(a, ") void {{\n    var m = n;\n    var xp = x;\n    var ix = inc_x;\n    var yp = y;\n    var iy = inc_y;\n", .{}) catch @panic("OOM");
    if (k.prec.complex()) {
        w.print(a, "    const al = [2]{s}{{ alpha_r, alpha_i }};\n", .{F}) catch @panic("OOM");
    } else {
        w.print(a, "    const al = [1]{s}{{alpha}};\n", .{F}) catch @panic("OOM");
    }
    w.print(a, "    var ap: [*]const {s} = &al;\n    asm volatile (\n", .{F}) catch @panic("OOM");
    for (lines) |l| w.print(a, "        \\\\{s}\n", .{l}) catch @panic("OOM");
    for (operands, 0..) |o, i| {
        w.print(a, "        {s} [{s}] \"+{{{s}}}\" ({s}),\n", .{ if (i == 0) ":" else " ", o[0], o[1], o[2] }) catch @panic("OOM");
    }
    w.print(a, "        :\n        : .{{ .memory = true, .cc = true", .{}) catch @panic("OOM");
    for (clobbers(lines)) |c| w.print(a, ", .{s} = true", .{c}) catch @panic("OOM");
    w.print(a, " }});\n}}\n", .{}) catch @panic("OOM");
}

fn same(a: []const []const u8, b: []const []const u8) bool {
    if (a.len != b.len) return false;
    for (a, b) |x, y| if (!std.mem.eql(u8, x, y)) return false;
    return true;
}

/// Whether `core` uses `k.source` for `k`'s precision.
fn uses(core: []const u8, k: Kernel) bool {
    for (map.resolve(arch, core, k.letter())) |s| if (std.mem.eql(u8, s, k.source)) return true;
    return false;
}

fn coreOf(name: []const u8) Core {
    // `OPTERON_SSE3` is `OPTERON` with SSE3; the headers know one macro
    const sse3 = !std.mem.eql(u8, name, "OPTERON");
    const macro = if (std.mem.eql(u8, name, "OPTERON_SSE3")) "OPTERON" else name;
    return .{ .name = name, .macro = macro, .sse3 = sse3 };
}

pub fn main(init: std.process.Init) !void {
    const arena = init.arena.allocator();
    const checkout = init.environ_map.get("OPENBLAS") orelse
        map.die("set OPENBLAS to the OpenBLAS checkout", .{});
    map.init(arena, init.io, checkout);
    const io = init.io;

    // the two headers the sources include; nothing in them matters here
    const inc = ".zig-cache/blas_asm_inc";
    try Io.Dir.cwd().createDirPath(io, inc);
    try Io.Dir.cwd().writeFile(io, .{
        .sub_path = inc ++ "/config.h",
        .data = "#define ARCH_X86_64 1\n#define BINARY64 1\n",
    });

    const dir = "src/linalg/blas/blas1/axpy/x86_64/asm";
    try Io.Dir.cwd().createDirPath(io, dir);
    const all_cores = map.cores(arch).?;

    for (kernels) |k| {
        // group the cores by the text they expand to
        const Group = struct { text: []const []const u8, names: std.ArrayList(u8) = .empty, first: []const u8 };
        var groups: std.ArrayList(Group) = .empty;
        outer: for (all_cores) |name| {
            if (!uses(name, k)) continue;
            const c = coreOf(name);
            const lines = body(k, preprocess(inc, k, c));
            for (groups.items) |*g| {
                if (same(g.text, lines)) {
                    try g.names.print(arena, ", {s}", .{lowered(c.name)});
                    continue :outer;
                }
            }
            var g: Group = .{ .text = lines, .first = lowered(c.name) };
            try g.names.print(arena, "{s}", .{lowered(c.name)});
            try groups.append(arena, g);
        }

        for (groups.items) |g| {
            var out: std.ArrayList(u8) = .empty;
            emit(&out, k, g.names.items, g.text);
            const file = try std.fmt.allocPrint(arena, "{s}{s}_{s}.zig", .{ k.stem(), if (k.conj) "_conj" else "", g.first });
            const dest = map.path(&.{ dir, file });
            try Io.Dir.cwd().writeFile(io, .{ .sub_path = dest, .data = out.items });
            _ = try std.process.run(arena, io, .{ .argv = &.{ "zig", "fmt", dest } });
            std.debug.print("wrote {s} ({d} body lines): {s}\n", .{ dest, g.text.len, g.names.items });
        }
    }
}

fn lowered(s: []const u8) []const u8 {
    return std.ascii.allocLowerString(map.arena, s) catch @panic("OOM");
}
