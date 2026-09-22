//! 32-bit ARM with VFP. The source has two argument paths, one for the
//! hard-float ABI and one for the soft; it is expanded on the hard-float path
//! and the arguments are delivered to the registers that path uses, whichever
//! ABI the Zig function is built for: `n`, `y`, `incx`, `x`, `incy` in `r0` to
//! `r4`, where the source loads three of them from the stack itself, and
//! `alpha` in `s0`/`d0` (and `s1`/`d1`), moved from integer registers the
//! source leaves alone.
//!
//! The source saves and restores what the compiler now saves for it, and
//! builds a frame to do so; all of that is dropped. It assumes ARM state.

const std = @import("std");
const arch = @import("__root__.zig");
const Error = arch.Error;

const Kernel = arch.Kernel;
const Pin = arch.Pin;

const ints = [_]Pin{
    .{ .name = "n", .reg = "r0", .variable = "m", .init = "n" },
    .{ .name = "y", .reg = "r1", .variable = "yp", .init = "y" },
    .{ .name = "incx", .reg = "r2", .variable = "ix", .init = "inc_x" },
    .{ .name = "x", .reg = "r3", .variable = "xp", .init = "x" },
    .{ .name = "incy", .reg = "r4", .variable = "iy", .init = "inc_y" },
};
const real_s = ints ++ [_]Pin{
    .{ .name = "alpha", .reg = "r5", .variable = "da", .init = "@as(u32, @bitCast(alpha))", .inout = false },
};
const real_d = ints ++ [_]Pin{
    .{ .name = "alpha_lo", .reg = "r5", .variable = "da_lo", .init = "@as(u32, @truncate(@as(u64, @bitCast(alpha))))", .inout = false },
    .{ .name = "alpha_hi", .reg = "r6", .variable = "da_hi", .init = "@as(u32, @truncate(@as(u64, @bitCast(alpha)) >> 32))", .inout = false },
};
const cplx_s = ints ++ [_]Pin{
    .{ .name = "alpha_r", .reg = "r5", .variable = "dar", .init = "@as(u32, @bitCast(alpha_r))", .inout = false },
    .{ .name = "alpha_i", .reg = "r7", .variable = "dai", .init = "@as(u32, @bitCast(alpha_i))", .inout = false },
};
const cplx_d = ints ++ [_]Pin{
    .{ .name = "alpha_r_lo", .reg = "r5", .variable = "dar_lo", .init = "@as(u32, @truncate(@as(u64, @bitCast(alpha_r))))", .inout = false },
    .{ .name = "alpha_r_hi", .reg = "r6", .variable = "dar_hi", .init = "@as(u32, @truncate(@as(u64, @bitCast(alpha_r)) >> 32))", .inout = false },
    .{ .name = "alpha_i_lo", .reg = "r7", .variable = "dai_lo", .init = "@as(u32, @truncate(@as(u64, @bitCast(alpha_i))))", .inout = false },
    .{ .name = "alpha_i_hi", .reg = "r8", .variable = "dai_hi", .init = "@as(u32, @truncate(@as(u64, @bitCast(alpha_i)) >> 32))", .inout = false },
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
    var out: std.ArrayList([]const u8) = .empty;
    try out.append(a, try std.fmt.allocPrint(a, "-D{s}", .{core}));
    try out.appendSlice(a, &.{ "-DARMV7", "-DHAVE_VFP", "-DHAVE_VFPV3", "-D__ARM_PCS_VFP" });
    return out.items;
}

fn has(line: []const u8, needle: []const u8) bool {
    return std.mem.indexOf(u8, line, needle) != null;
}

/// Whether `line` mentions register `name` as a whole word.
fn mentions(line: []const u8, name: []const u8) bool {
    var i: usize = 0;
    while (std.mem.indexOfPos(u8, line, i, name)) |p| {
        const before = p == 0 or !std.ascii.isAlphanumeric(line[p - 1]);
        const after = p + name.len >= line.len or !std.ascii.isAlphanumeric(line[p + name.len]);
        if (before and after) return true;
        i = p + name.len;
    }
    return false;
}

fn transform(a: std.mem.Allocator, k: Kernel, lines: []const []const u8) Error![]const []const u8 {
    var out: std.ArrayList([]const u8) = .empty;
    if (k.prec.double()) {
        try out.append(a, " vmov d0, r5, r6");
        if (k.prec.complex()) try out.append(a, " vmov d1, r7, r8");
    } else {
        try out.append(a, " vmov s0, r5");
        if (k.prec.complex()) try out.append(a, " vmov s1, r7");
    }
    for (lines) |line| {
        const t = std.mem.trim(u8, line, " ");
        // what declares the function
        if (has(t, ".global") or has(t, ".type k") or std.mem.startsWith(u8, t, ".section .note") or std.mem.startsWith(u8, t, ".size")) continue;
        if (std.mem.eql(u8, t, ".arm")) continue;
        // the frame, the arguments the source loads from the stack, and the
        // registers it saves; the compiler saves what is declared clobbered
        if (mentions(t, "fp") or mentions(t, "sp") or std.mem.startsWith(u8, t, "push ") or std.mem.startsWith(u8, t, "pop ")) continue;
        if (std.mem.startsWith(u8, t, "vstm ") or std.mem.startsWith(u8, t, "vldm ")) continue;
        // the return
        if (std.mem.eql(u8, t, "bx lr")) {
            try out.append(a, " b " ++ arch.end_label ++ "f");
            continue;
        }
        try out.append(a, line);
    }
    try out.append(a, arch.end_label ++ ":");
    return out.items;
}

/// The double register a name belongs to: `s5` is half of `d2`, `q1` is
/// `d2` and `d3`. Returns the first; `q` is handled by the caller.
fn canonical(a: std.mem.Allocator, reg: []const u8) !?[]const u8 {
    if (reg.len < 2 or !std.ascii.isDigit(reg[1])) return null;
    const n = std.fmt.parseInt(u8, reg[1..], 10) catch return null;
    return switch (reg[0]) {
        'r' => if (n <= 12) try std.fmt.allocPrint(a, "r{d}", .{n}) else null,
        'd' => if (n <= 31) try std.fmt.allocPrint(a, "d{d}", .{n}) else null,
        's' => if (n <= 31) try std.fmt.allocPrint(a, "d{d}", .{n / 2}) else null,
        else => null,
    };
}

/// The registers named in `line`, with the `{ d4 - d7 }` ranges of
/// load and store multiples filled in.
fn scan(a: std.mem.Allocator, line: []const u8, seen: *std.StringArrayHashMapUnmanaged(void)) !void {
    var i: usize = 0;
    var prev: ?[]const u8 = null; // the last register, to close a range
    var range = false;
    while (i < line.len) {
        if (!std.ascii.isAlphabetic(line[i]) or (i > 0 and (std.ascii.isAlphanumeric(line[i - 1]) or line[i - 1] == '_' or line[i - 1] == '.'))) {
            if (line[i] == '-') range = true;
            i += 1;
            continue;
        }
        var e = i;
        while (e < line.len and std.ascii.isAlphanumeric(line[e])) e += 1;
        const tok = line[i..e];
        i = e;
        const c = try canonical(a, tok) orelse continue;
        if (range and prev != null and c[0] == prev.?[0]) {
            const lo = std.fmt.parseInt(u8, prev.?[1..], 10) catch 0;
            const hi = std.fmt.parseInt(u8, c[1..], 10) catch 0;
            var n = lo;
            while (n <= hi) : (n += 1) try seen.put(a, try std.fmt.allocPrint(a, "{c}{d}", .{ c[0], n }), {});
        }
        range = false;
        prev = c;
        try seen.put(a, c, {});
    }
}

fn clobbers(a: std.mem.Allocator, text: []const []const u8, pinned: []const Pin) Error![]const []const u8 {
    var seen: std.StringArrayHashMapUnmanaged(void) = .empty;
    for (text) |line| try scan(a, line, &seen);
    for (pinned) |p| _ = seen.swapRemove(p.reg);
    return seen.keys();
}

pub const arm: arch.Arch = .{
    .name = "arm",
    .triple = "arm-linux-gnueabihf",
    .config_h = "#define ARCH_ARM 1\n#define BINARY32 1\n",
    .defines = &.{ "-DASSEMBLER", "-DARCH_ARM", "-DOS_LINUX", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k" },
    .core_defines = coreDefines,
    .pins = pins,
    .decls = decls,
    .transform = transform,
    .clobbers = clobbers,
    .always_clobbers = &.{ "memory", "cpsr" },
};
