//! PowerPC and POWER. Each source names its argument registers in `#define`s
//! (`N`, `X`, `INCX`, `Y`, `INCY`), one set per pointer width, which are read
//! from the preprocessor rather than copied here; the two builds are two
//! variants. `alpha` is in `f1` (and `f2`) on entry.
//!
//! An argument that does not fit in registers, the stride of `y`, is read from
//! the caller's frame by the source. Zig delivers it in `r12` instead, and the
//! load becomes a register move. The registers the source saves in its own
//! frame, and the frame, are dropped: the compiler saves what is declared
//! clobbered.
//!
//! Registers are bare numbers in this assembly, `fmr 24, 1` being `f24` and
//! `f1`, so a register's class cannot be read off the text; every number an
//! instruction mentions is taken to be both a general and a floating-point
//! register. That over-declares, which costs the compiler a few saves and is
//! otherwise harmless.

const std = @import("std");
const arch = @import("__root__.zig");
const Error = arch.Error;

const Kernel = arch.Kernel;
const Pin = arch.Pin;

/// Where the stride of `y` arrives; the sources do not use it.
const spare = "r12";

fn macro(name: []const u8) ?[]const u8 {
    return arch.macros.get(name);
}

const power8_ints = [_]Pin{
    .{ .name = "n", .reg = "r3", .variable = "m", .init = "n" },
    .{ .name = "x", .reg = "r8", .variable = "xp", .init = "x" },
    .{ .name = "incx", .reg = "r9", .variable = "ix", .init = "inc_x" },
    .{ .name = "y", .reg = "r10", .variable = "yp", .init = "y" },
    .{ .name = "incy", .reg = "r12", .variable = "iy", .init = "inc_y", .inout = false },
};

fn pins(a: std.mem.Allocator, k: Kernel) Error![]const Pin {
    var out: std.ArrayList(Pin) = .empty;
    if (macro("N")) |n| {
        try out.append(a, .{ .name = "n", .reg = n, .variable = "m", .init = "n" });
        try out.append(a, .{ .name = "x", .reg = macro("X").?, .variable = "xp", .init = "x" });
        try out.append(a, .{ .name = "incx", .reg = macro("INCX").?, .variable = "ix", .init = "inc_x" });
        try out.append(a, .{ .name = "y", .reg = macro("Y").?, .variable = "yp", .init = "y" });
        try out.append(a, .{ .name = "incy", .reg = macro("INCY").?, .variable = "iy", .init = "inc_y" });
        try out.append(a, .{ .name = "incy_arg", .reg = spare, .variable = "iy2", .init = "inc_y", .inout = false });
    } else {
        // compiler-generated sources with no macros: the ELFv2 argument registers
        try out.appendSlice(a, &power8_ints);
    }
    if (k.prec.complex()) {
        try out.append(a, .{ .name = "alpha_r", .reg = "f1", .variable = "dar", .init = "alpha_r" });
        try out.append(a, .{ .name = "alpha_i", .reg = "f2", .variable = "dai", .init = "alpha_i" });
    } else {
        try out.append(a, .{ .name = "alpha", .reg = "f1", .variable = "da", .init = "alpha" });
    }
    return out.items;
}

fn decls(_: std.mem.Allocator, _: Kernel) Error![]const []const u8 {
    return &.{};
}

fn coreDefines(a: std.mem.Allocator, _: Kernel, core: []const u8) Error![]const []const u8 {
    const out = try a.alloc([]const u8, 1);
    out[0] = try std.fmt.allocPrint(a, "-D{s}", .{core});
    return out;
}

/// The first integer operand of `t` after its mnemonic, if the whole operand
/// is a base register in `(1)` form, and the displacement before it.
fn stackOperand(t: []const u8) ?struct { rt: []const u8, off: i64 } {
    const close = std.mem.lastIndexOf(u8, t, "(1)") orelse return null;
    if (close + 3 != t.len) return null;
    const sp = std.mem.indexOfScalar(u8, t, ' ') orelse return null;
    var it = std.mem.splitScalar(u8, t[sp + 1 .. close], ',');
    const rt = std.mem.trim(u8, it.next() orelse return null, " ");
    const disp = std.mem.trim(u8, it.next() orelse return null, " ");
    const off = std.fmt.parseInt(i64, disp, 10) catch return null;
    return .{ .rt = rt, .off = off };
}

fn mnemonic(t: []const u8) []const u8 {
    const end = std.mem.indexOfAny(u8, t, " \t") orelse t.len;
    return t[0..end];
}

fn has(line: []const u8, needle: []const u8) bool {
    return std.mem.indexOf(u8, line, needle) != null;
}

fn transform(a: std.mem.Allocator, _: Kernel, lines: []const []const u8) Error![]const []const u8 {
    var out: std.ArrayList([]const u8) = .empty;
    for (lines) |line| {
        const t = std.mem.trim(u8, line, " ");
        // what declares the function
        const declares = [_][]const u8{ ".globl", ".hidden", ".type", ".size", ".abiversion", ".file", ".localentry", ".previous", ".TOC." };
        var skip = false;
        for (declares) |d| skip = skip or has(t, d);
        if (std.mem.startsWith(u8, t, ".section") or std.mem.startsWith(u8, t, ".machine")) skip = true;
        if (std.mem.startsWith(u8, t, "caxpy") and t[t.len - 1] == ':') skip = true;
        if (std.mem.eql(u8, t, ".LCF0:")) skip = true;
        if (skip) continue;
        // the frame: arguments read from it become a move, the rest is dropped
        if (has(t, "(1)")) {
            const m = mnemonic(t);
            if ((std.mem.eql(u8, m, "ld") or std.mem.eql(u8, m, "lwz"))) {
                if (stackOperand(t)) |op| if (op.off >= 96) {
                    try out.append(a, try std.fmt.allocPrint(a, " mr {s}, 12", .{op.rt}));
                };
            }
            continue;
        }
        if (std.mem.startsWith(u8, t, "addi 1, 1,")) continue;
        // the return
        if (std.mem.eql(u8, t, "blr")) {
            try out.append(a, " b " ++ arch.end_label ++ "f");
            continue;
        }
        try out.append(a, line);
    }
    try out.append(a, arch.end_label ++ ":");
    return out.items;
}

fn isVector(m: []const u8) bool {
    const prefixes = [_][]const u8{ "lxv", "stxv", "xx", "xv", "xs", "lxs", "stxs", "v" };
    for (prefixes) |p| if (std.mem.startsWith(u8, m, p)) return true;
    return false;
}

fn clobbers(a: std.mem.Allocator, text: []const []const u8, pinned: []const Pin) Error![]const []const u8 {
    var seen: std.StringArrayHashMapUnmanaged(void) = .empty;
    // the volatile registers and the condition and count registers are
    // always the body's to use
    var n: u8 = 0;
    while (n < 14) : (n += 1) {
        if (n != 1 and n != 2 and n != 13) try seen.put(a, try std.fmt.allocPrint(a, "r{d}", .{n}), {});
        try seen.put(a, try std.fmt.allocPrint(a, "f{d}", .{n}), {});
    }
    for ([_][]const u8{ "ctr", "xer", "cr0", "cr1", "cr5", "cr6", "cr7" }) |r| try seen.put(a, r, {});
    var vector = false;
    for (text) |line| {
        const t = std.mem.trim(u8, line, " ");
        if (t.len == 0 or t[0] == '.' or t[t.len - 1] == ':') continue;
        const m = mnemonic(t);
        if (isVector(m)) vector = true;
        // every number in the operands
        var i: usize = m.len;
        while (i < t.len) {
            if (std.ascii.isDigit(t[i]) and (i == 0 or !(std.ascii.isAlphanumeric(t[i - 1]) or t[i - 1] == '.' or t[i - 1] == '_'))) {
                var e = i;
                while (e < t.len and std.ascii.isDigit(t[e])) e += 1;
                if (e < t.len and std.ascii.isAlphabetic(t[e])) {
                    i = e;
                    continue;
                }
                const v = std.fmt.parseInt(u8, t[i..e], 10) catch 255;
                if (v <= 31 and v != 1 and v != 2 and v != 13) {
                    try seen.put(a, try std.fmt.allocPrint(a, "r{d}", .{v}), {});
                    try seen.put(a, try std.fmt.allocPrint(a, "f{d}", .{v}), {});
                }
                i = e;
            } else i += 1;
        }
    }
    if (vector) {
        var v: u8 = 0;
        while (v < 64) : (v += 1) try seen.put(a, try std.fmt.allocPrint(a, "vs{d}", .{v}), {});
        v = 0;
        while (v < 32) : (v += 1) try seen.put(a, try std.fmt.allocPrint(a, "v{d}", .{v}), {});
    }
    for (pinned) |p| _ = seen.swapRemove(p.reg);
    return seen.keys();
}

pub const power: arch.Arch = .{
    .name = "power",
    .triple = "powerpc64le-linux-gnu",
    .config_h = "#define ARCH_POWER 1\n#define BINARY64 1\n",
    .defines = &.{ "-DASSEMBLER", "-DARCH_POWER", "-DOS_LINUX", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k", "-D__BYTE_ORDER__=__ORDER_LITTLE_ENDIAN__" },
    .core_defines = coreDefines,
    .pins = pins,
    .decls = decls,
    .transform = transform,
    .clobbers = clobbers,
    .always_clobbers = &.{"memory"},
    // the Blue Gene double-precision FPU instructions are not in LLVM
    .unsupported = &.{ "axpy_hummer.S", "zaxpy_hummer.S" },
    .variants = &.{
        .{ .name = "32", .defines = &.{} },
        .{ .name = "64", .defines = &.{"-D__64BIT__"} },
    },
};
