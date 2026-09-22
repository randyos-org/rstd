//! AArch64. The sources read their arguments straight from the AAPCS64
//! registers, which the function's own C ABI delivers them in: the pins name
//! those registers and cost no move. Every source here is assembly with no
//! stack use, so only the registers it touches are declared.

const std = @import("std");
const arch = @import("__root__.zig");
const Error = arch.Error;

const Kernel = arch.Kernel;
const Pin = arch.Pin;

const real_s = [_]Pin{
    .{ .name = "n", .reg = "x0", .variable = "m", .init = "n" },
    .{ .name = "x", .reg = "x3", .variable = "xp", .init = "x" },
    .{ .name = "incx", .reg = "x4", .variable = "ix", .init = "inc_x" },
    .{ .name = "y", .reg = "x5", .variable = "yp", .init = "y" },
    .{ .name = "incy", .reg = "x6", .variable = "iy", .init = "inc_y" },
    .{ .name = "alpha", .reg = "s0", .variable = "da", .init = "alpha" },
};
const real_d = real_s[0..5].* ++ [_]Pin{
    .{ .name = "alpha", .reg = "d0", .variable = "da", .init = "alpha" },
};
const cplx_s = real_s[0..5].* ++ [_]Pin{
    .{ .name = "alpha_r", .reg = "s0", .variable = "dar", .init = "alpha_r" },
    .{ .name = "alpha_i", .reg = "s1", .variable = "dai", .init = "alpha_i" },
};
const cplx_d = real_s[0..5].* ++ [_]Pin{
    .{ .name = "alpha_r", .reg = "d0", .variable = "dar", .init = "alpha_r" },
    .{ .name = "alpha_i", .reg = "d1", .variable = "dai", .init = "alpha_i" },
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

fn transform(a: std.mem.Allocator, _: Kernel, lines: []const []const u8) Error![]const []const u8 {
    var out: std.ArrayList([]const u8) = .empty;
    var in_prologue = false;
    for (lines) |line| {
        const t = std.mem.trim(u8, line, " ");
        // the macro that declares the function, and its use
        if (std.mem.startsWith(u8, t, ".macro PROLOGUE")) {
            in_prologue = true;
            continue;
        }
        if (in_prologue) {
            if (std.mem.eql(u8, t, ".endm")) in_prologue = false;
            continue;
        }
        if (std.mem.eql(u8, t, "PROLOGUE")) continue;
        if (std.mem.startsWith(u8, t, ".size ") or std.mem.startsWith(u8, t, ".globl ") or std.mem.startsWith(u8, t, ".global ")) continue;
        // `ret` leaves the Zig function without its epilogue; branch to the end
        if (std.mem.eql(u8, t, "ret")) {
            try out.append(a, " b " ++ arch.end_label ++ "f");
            continue;
        }
        try out.append(a, line);
    }
    try out.append(a, arch.end_label ++ ":");
    return out.items;
}

/// `x5`, `w5` to `x5`; `v3`, `q3`, `d3`, `s3`, `h3`, `b3` to `v3`.
fn canonical(a: std.mem.Allocator, reg: []const u8) !?[]const u8 {
    if (reg.len < 2 or !std.ascii.isDigit(reg[1])) return null;
    const n = std.fmt.parseInt(u8, reg[1..], 10) catch return null;
    return switch (reg[0]) {
        'x', 'w' => if (n <= 30) try std.fmt.allocPrint(a, "x{d}", .{n}) else null,
        'v', 'q', 'd', 's', 'h', 'b' => if (n <= 31) try std.fmt.allocPrint(a, "v{d}", .{n}) else null,
        else => null,
    };
}

fn clobbers(a: std.mem.Allocator, text: []const []const u8, pinned: []const Pin) Error![]const []const u8 {
    var seen: std.StringArrayHashMapUnmanaged(void) = .empty;
    for (text) |line| {
        var i: usize = 0;
        while (i < line.len) {
            if (!std.ascii.isAlphabetic(line[i]) or (i > 0 and (std.ascii.isAlphanumeric(line[i - 1]) or line[i - 1] == '_' or line[i - 1] == '.'))) {
                i += 1;
                continue;
            }
            var e = i;
            while (e < line.len and std.ascii.isAlphanumeric(line[e])) e += 1;
            const name = try canonical(a, line[i..e]);
            i = e;
            if (name == null) continue;
            const is_pin = for (pinned) |p| {
                if (try canonical(a, p.reg)) |c| if (std.mem.eql(u8, c, name.?)) break true;
            } else false;
            if (!is_pin) try seen.put(a, name.?, {});
        }
    }
    return seen.keys();
}

pub const arm64: arch.Arch = .{
    .name = "arm64",
    .triple = "aarch64-linux-gnu",
    .config_h = "#define ARCH_ARM64 1\n#define BINARY64 1\n",
    .defines = &.{ "-DASSEMBLER", "-DARCH_ARM64", "-DOS_LINUX", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k" },
    .core_defines = coreDefines,
    .pins = pins,
    .decls = decls,
    .transform = transform,
    .clobbers = clobbers,
    .always_clobbers = &.{ "memory", "nzcv" },
};
