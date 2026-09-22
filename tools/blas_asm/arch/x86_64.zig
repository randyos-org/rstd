//! x86_64: the SSE, SSE2 and Atom kernels.
//!
//! The sources read their arguments from the registers the SysV ABI puts them
//! in, except the stride of `y`, which they load from the stack; the operand
//! list pins the same registers and the Zig function delivers into them on any
//! OS. `alpha` is passed by address and loaded into `xmm0` (and `xmm1`), which
//! is where the source's own Linux path finds it.

const std = @import("std");
const arch = @import("__root__.zig");
const Error = arch.Error;

const Kernel = arch.Kernel;
const Pin = arch.Pin;

const pins_real = [_]Pin{
    .{ .name = "m", .reg = "rdi", .variable = "m", .init = "n" },
    .{ .name = "x", .reg = "rcx", .variable = "xp", .init = "x" },
    .{ .name = "incx", .reg = "r8", .variable = "ix", .init = "inc_x" },
    .{ .name = "y", .reg = "r9", .variable = "yp", .init = "y" },
    .{ .name = "incy", .reg = "rsi", .variable = "iy", .init = "inc_y" },
    .{ .name = "ap", .reg = "rdx", .variable = "ap", .init = null },
};

fn pins(_: std.mem.Allocator, _: Kernel) Error![]const Pin {
    return &pins_real;
}

fn decls(a: std.mem.Allocator, k: Kernel) Error![]const []const u8 {
    const F = k.prec.float();
    const arr = if (k.prec.complex())
        try std.fmt.allocPrint(a, "const al = [2]{s}{{ alpha_r, alpha_i }};", .{F})
    else
        try std.fmt.allocPrint(a, "const al = [1]{s}{{alpha}};", .{F});
    const ptr = try std.fmt.allocPrint(a, "var ap: [*]const {s} = &al;", .{F});
    const out = try a.alloc([]const u8, 2);
    out[0] = arr;
    out[1] = ptr;
    return out;
}

fn coreDefines(a: std.mem.Allocator, _: Kernel, core: []const u8) Error![]const []const u8 {
    // `OPTERON_SSE3` is `OPTERON` with SSE3; the headers know one macro
    const macro = if (std.mem.eql(u8, core, "OPTERON_SSE3")) "OPTERON" else core;
    var out: std.ArrayList([]const u8) = .empty;
    try out.append(a, try std.fmt.allocPrint(a, "-D{s}", .{macro}));
    try out.append(a, try std.fmt.allocPrint(a, "-DCORE_{s}", .{macro}));
    try out.appendSlice(a, &.{ "-DHAVE_SSE", "-DHAVE_SSE2" });
    if (!std.mem.eql(u8, core, "OPTERON")) try out.append(a, "-DHAVE_SSE3");
    return out.items;
}

fn has(line: []const u8, needle: []const u8) bool {
    return std.mem.indexOf(u8, line, needle) != null;
}

fn transform(a: std.mem.Allocator, k: Kernel, lines: []const []const u8) Error![]const []const u8 {
    var out: std.ArrayList([]const u8) = .empty;
    // `alpha` arrives by pointer; the real part goes where the source's
    // Linux path expects it, `xmm0`, and the imaginary part in `xmm1`
    const load = if (k.prec.double()) " movsd" else " movss";
    const step: usize = if (k.prec.double()) 8 else 4;
    try out.append(a, try std.fmt.allocPrint(a, "{s} (%rdx), %xmm0", .{load}));
    if (k.prec.complex()) try out.append(a, try std.fmt.allocPrint(a, "{s} {d}(%rdx), %xmm1", .{ load, step }));
    for (lines) |line| {
        // the directives that declare the function, and the ABI shuffling
        if (has(line, ".globl") or has(line, ".size k") or has(line, ".section .note")) continue;
        if (has(line, "movq 8(%rsp)")) continue;
        if (has(line, "(%rsp)")) std.debug.panic("the body reads the stack: {s}", .{line});
        // `ret` leaves the Zig function without its epilogue; jump to the end
        if (std.mem.eql(u8, std.mem.trim(u8, line, " "), "ret")) {
            try out.append(a, " jmp " ++ arch.end_label ++ "f");
            continue;
        }
        try out.append(a, line);
    }
    try out.append(a, arch.end_label ++ ":");
    return out.items;
}

/// General registers by every name they go by, to the 64-bit one.
fn canonical(a: std.mem.Allocator, reg: []const u8) !?[]const u8 {
    const legacy = [_][]const u8{ "ax", "bx", "cx", "dx", "si", "di", "bp" };
    for (legacy) |l| {
        const wide = try std.fmt.allocPrint(a, "r{s}", .{l});
        if (std.mem.eql(u8, reg, wide) or (std.mem.eql(u8, reg[1..], l) and reg[0] == 'e')) return wide;
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
        if (n >= 8 and n <= 15) return try std.fmt.allocPrint(a, "r{d}", .{n});
    }
    return null;
}

fn clobbers(a: std.mem.Allocator, text: []const []const u8, pinned: []const Pin) Error![]const []const u8 {
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
            const name = if (is_vec) reg else try canonical(a, reg) orelse continue;
            if (std.mem.eql(u8, name, "rsp")) continue;
            const is_pin = for (pinned) |o| {
                if (std.mem.eql(u8, o.reg, name)) break true;
            } else false;
            if (!is_pin) try seen.put(a, name, {});
        }
    }
    return seen.keys();
}

pub const x86_64: arch.Arch = .{
    .name = "x86_64",
    .triple = "x86_64-linux-gnu",
    .config_h = "#define ARCH_X86_64 1\n#define BINARY64 1\n",
    .defines = &.{ "-DASSEMBLER", "-DARCH_X86_64", "-DOS_LINUX", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k" },
    .core_defines = coreDefines,
    .pins = pins,
    .decls = decls,
    .transform = transform,
    .clobbers = clobbers,
    .always_clobbers = &.{ "memory", "cc" },
};
