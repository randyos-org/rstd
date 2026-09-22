//! MIPS64 (n64). The sources read their arguments from the n64 argument
//! registers, which the function's own C ABI delivers them in: `n` in `$4`,
//! `x`, `incx`, `y`, `incy` in `$8` to `$11`, and `alpha`, the fourth
//! argument, in `$f15`. As on LoongArch, LLVM cannot constrain a double to a
//! register here, so `alpha` goes in as bits through `$12`, which the sources
//! leave alone, and the preamble moves it.
//!
//! The sources run with `.set noreorder`, scheduling their own delay slots, so
//! the block does too, and restores the assembler's mode after. Register
//! numbers are also written into `.word`s for the Loongson extensions, where
//! no scan of the text finds them; every floating-point register is declared.

const std = @import("std");
const arch = @import("__root__.zig");
const Error = arch.Error;

const Kernel = arch.Kernel;
const Pin = arch.Pin;

const ints = [_]Pin{
    .{ .name = "n", .reg = "$4", .variable = "m", .init = "n" },
    .{ .name = "x", .reg = "$8", .variable = "xp", .init = "x" },
    .{ .name = "incx", .reg = "$9", .variable = "ix", .init = "inc_x" },
    .{ .name = "y", .reg = "$10", .variable = "yp", .init = "y" },
    .{ .name = "incy", .reg = "$11", .variable = "iy", .init = "inc_y" },
};
const real_s = ints ++ [_]Pin{
    .{ .name = "alpha", .reg = "$12", .variable = "dab", .init = "@as(u32, @bitCast(alpha))", .inout = false },
};
const real_d = ints ++ [_]Pin{
    .{ .name = "alpha", .reg = "$12", .variable = "dab", .init = "@as(u64, @bitCast(alpha))", .inout = false },
};
const cplx_s = ints ++ [_]Pin{
    .{ .name = "alpha_r", .reg = "$12", .variable = "dar", .init = "@as(u32, @bitCast(alpha_r))", .inout = false },
    .{ .name = "alpha_i", .reg = "$13", .variable = "dai", .init = "@as(u32, @bitCast(alpha_i))", .inout = false },
};
const cplx_d = ints ++ [_]Pin{
    .{ .name = "alpha_r", .reg = "$12", .variable = "dar", .init = "@as(u64, @bitCast(alpha_r))", .inout = false },
    .{ .name = "alpha_i", .reg = "$13", .variable = "dai", .init = "@as(u64, @bitCast(alpha_i))", .inout = false },
};

fn pins(_: std.mem.Allocator, k: Kernel) Error![]const Pin {
    return switch (k.prec) {
        .s => &real_s,
        .d => &real_d,
        .c => &cplx_s,
        .z => &cplx_d,
    };
}

fn decls(_: std.mem.Allocator, _: Kernel) Error![]const []const u8 {
    return &.{};
}

fn coreDefines(a: std.mem.Allocator, _: Kernel, core: []const u8) Error![]const []const u8 {
    const out = try a.alloc([]const u8, 1);
    out[0] = try std.fmt.allocPrint(a, "-D{s}", .{core});
    return out;
}

fn has(line: []const u8, needle: []const u8) bool {
    return std.mem.indexOf(u8, line, needle) != null;
}

fn isBranch(line: []const u8) bool {
    const t = std.mem.trim(u8, line, " ");
    var it = std.mem.tokenizeScalar(u8, t, ' ');
    const op = it.next() orelse return false;
    return std.mem.startsWith(u8, op, "b") or std.mem.startsWith(u8, op, "j");
}

fn transform(a: std.mem.Allocator, k: Kernel, lines: []const []const u8) Error![]const []const u8 {
    var out: std.ArrayList([]const u8) = .empty;
    try out.appendSlice(a, &.{ " .set push", " .set mips64", " .set noreorder", " .set nomacro" });
    const move = if (k.prec.double()) "dmtc1" else "mtc1";
    try out.append(a, try std.fmt.allocPrint(a, " {s} $12, $f15", .{move}));
    if (k.prec.complex()) try out.append(a, try std.fmt.allocPrint(a, " {s} $13, $f14", .{move}));
    var prev_branch = false;
    for (lines) |line| {
        const t = std.mem.trim(u8, line, " ");
        // what declares the function, and the assembler modes it sets
        if (has(t, ".globl") or has(t, ".ent ") or has(t, ".end ") or has(t, ".size")) continue;
        // the source saves and restores what the compiler now saves
        if (has(t, "$sp")) {
            // in a delay slot, the slot still needs an instruction
            if (prev_branch) try out.append(a, " nop");
            prev_branch = false;
            continue;
        }
        // the return
        if (std.mem.eql(u8, t, "j $31") or std.mem.eql(u8, t, "jr $31")) {
            try out.append(a, " j " ++ arch.end_label ++ "f");
            prev_branch = true;
            continue;
        }
        try out.append(a, line);
        prev_branch = isBranch(line);
    }
    try out.append(a, arch.end_label ++ ":");
    try out.append(a, " .set pop");
    return out.items;
}

fn clobbers(a: std.mem.Allocator, text: []const []const u8, pinned: []const Pin) Error![]const []const u8 {
    var seen: std.StringArrayHashMapUnmanaged(void) = .empty;
    // every floating-point register, since `.word`s name them by number
    var f: u8 = 0;
    while (f < 32) : (f += 1) try seen.put(a, try std.fmt.allocPrint(a, "f{d}", .{f}), {});
    for (text) |line| {
        var i: usize = 0;
        while (std.mem.indexOfScalarPos(u8, line, i, '$')) |p| {
            var e = p + 1;
            while (e < line.len and std.ascii.isAlphanumeric(line[e])) e += 1;
            const reg = line[p + 1 .. e];
            i = e;
            const n = std.fmt.parseInt(u8, reg, 10) catch continue;
            // the zero register, `$sp`, `$gp` and `$ra` are not the body's to use
            if (n == 0 or n == 28 or n == 29 or n == 31) continue;
            const name = try std.fmt.allocPrint(a, "r{d}", .{n});
            try seen.put(a, name, {});
        }
    }
    // the pinned registers are operands, not clobbers
    for (pinned) |pin| {
        const name = if (std.mem.startsWith(u8, pin.reg, "$")) try std.fmt.allocPrint(a, "r{s}", .{pin.reg[1..]}) else pin.reg;
        _ = seen.swapRemove(name);
    }
    // `$f14` and `$f15` are set by the preamble but are still clobbered
    return seen.keys();
}

pub const mips64: arch.Arch = .{
    .name = "mips64",
    .triple = "mips64el-linux-gnuabi64",
    .config_h = "#define ARCH_MIPS64 1\n#define BINARY64 1\n",
    .defines = &.{ "-DASSEMBLER", "-DARCH_MIPS64", "-DOS_LINUX", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k" },
    .core_defines = coreDefines,
    .pins = pins,
    .decls = decls,
    .transform = transform,
    .clobbers = clobbers,
    .always_clobbers = &.{"memory"},
};
