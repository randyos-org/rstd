//! LoongArch64. The sources read their arguments from the LP64D argument
//! registers (`$r4` and `$r7` to `$r10`), which the function's own C ABI
//! delivers them in, so those pins cost no move.

const std = @import("std");
const arch = @import("__root__.zig");
const Error = arch.Error;

const Kernel = arch.Kernel;
const Pin = arch.Pin;

const ints = [_]Pin{
    .{ .name = "n", .reg = "r4", .variable = "m", .init = "n" },
    .{ .name = "x", .reg = "r7", .variable = "xp", .init = "x" },
    .{ .name = "incx", .reg = "r8", .variable = "ix", .init = "inc_x" },
    .{ .name = "y", .reg = "r9", .variable = "yp", .init = "y" },
    .{ .name = "incy", .reg = "r10", .variable = "iy", .init = "inc_y" },
};
// LLVM cannot yet constrain a double to a register here, so `alpha` goes in
// as bits through two integer registers the sources leave alone, and the
// preamble moves them to `$f0` and `$f1`, where the sources expect it.
const real_s = ints ++ [_]Pin{
    .{ .name = "alpha", .reg = "r11", .variable = "dab", .init = "@as(u32, @bitCast(alpha))", .inout = false },
};
const real_d = ints ++ [_]Pin{
    .{ .name = "alpha", .reg = "r11", .variable = "dab", .init = "@as(u64, @bitCast(alpha))", .inout = false },
};
const cplx_s = ints ++ [_]Pin{
    .{ .name = "alpha_r", .reg = "r11", .variable = "dar", .init = "@as(u32, @bitCast(alpha_r))", .inout = false },
    .{ .name = "alpha_i", .reg = "r18", .variable = "dai", .init = "@as(u32, @bitCast(alpha_i))", .inout = false },
};
const cplx_d = ints ++ [_]Pin{
    .{ .name = "alpha_r", .reg = "r11", .variable = "dar", .init = "@as(u64, @bitCast(alpha_r))", .inout = false },
    .{ .name = "alpha_i", .reg = "r18", .variable = "dai", .init = "@as(u64, @bitCast(alpha_i))", .inout = false },
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

fn transform(a: std.mem.Allocator, k: Kernel, lines: []const []const u8) Error![]const []const u8 {
    var out: std.ArrayList([]const u8) = .empty;
    const w = if (k.prec.double()) "d" else "w";
    try out.append(a, try std.fmt.allocPrint(a, " movgr2fr.{s} $f0, $r11", .{w}));
    if (k.prec.complex()) try out.append(a, try std.fmt.allocPrint(a, " movgr2fr.{s} $f1, $r18", .{w}));
    for (lines) |line| {
        const t = std.mem.trim(u8, line, " ");
        // what declares the function
        if (has(t, ".globl") or has(t, ".type k") or has(t, ".hidden") or std.mem.startsWith(u8, t, ".end") or std.mem.startsWith(u8, t, ".size")) continue;
        // the return
        if (std.mem.startsWith(u8, t, "jirl $r0, $r1")) {
            try out.append(a, " b " ++ arch.end_label ++ "f");
            continue;
        }
        try out.append(a, line);
    }
    try out.append(a, arch.end_label ++ ":");
    return out.items;
}

/// `$r12` to `r12`; `$f3`, `$vr3` and `$xr3`, one register in three widths,
/// to `v3`; `$fcc0` to `fcc0`. Null for the zero register, `$r3` (the stack
/// pointer) and anything else.
fn canonical(a: std.mem.Allocator, reg: []const u8) !?[]const u8 {
    const prefixes = [_]struct { p: []const u8, key: []const u8, max: u8 }{
        .{ .p = "fcc", .key = "fcc", .max = 7 },
        .{ .p = "vr", .key = "v", .max = 31 },
        .{ .p = "xr", .key = "v", .max = 31 },
        .{ .p = "r", .key = "r", .max = 31 },
        .{ .p = "f", .key = "v", .max = 31 },
    };
    for (prefixes) |pf| {
        if (!std.mem.startsWith(u8, reg, pf.p)) continue;
        const n = std.fmt.parseInt(u8, reg[pf.p.len..], 10) catch continue;
        if (n > pf.max) return null;
        if (std.mem.eql(u8, pf.key, "r") and (n == 0 or n == 3)) return null;
        return try std.fmt.allocPrint(a, "{s}{d}", .{ pf.key, n });
    }
    return null;
}

fn clobbers(a: std.mem.Allocator, text: []const []const u8, pinned: []const Pin) Error![]const []const u8 {
    var seen: std.StringArrayHashMapUnmanaged(void) = .empty;
    for (text) |line| {
        var i: usize = 0;
        while (std.mem.indexOfScalarPos(u8, line, i, '$')) |p| {
            var e = p + 1;
            while (e < line.len and std.ascii.isAlphanumeric(line[e])) e += 1;
            const reg = line[p + 1 .. e];
            i = e;
            const key = try canonical(a, reg) orelse continue;
            const is_pin = for (pinned) |pin| {
                if (try canonical(a, pin.reg)) |c| if (std.mem.eql(u8, c, key)) break true;
            } else false;
            if (!is_pin) try seen.put(a, reg, {});
        }
    }
    return seen.keys();
}

pub const loongarch64: arch.Arch = .{
    .name = "loongarch64",
    .triple = "loongarch64-linux-gnu",
    .config_h = "#define ARCH_LOONGARCH64 1\n#define BINARY64 1\n",
    .defines = &.{ "-DASSEMBLER", "-DARCH_LOONGARCH64", "-DOS_LINUX", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k" },
    .core_defines = coreDefines,
    .pins = pins,
    .decls = decls,
    .transform = transform,
    .clobbers = clobbers,
    .always_clobbers = &.{"memory"},
};
