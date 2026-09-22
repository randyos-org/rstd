//! 32-bit x86: the x87 kernels. The source reads every argument from the
//! stack, saves the registers it uses and returns itself, so each function is
//! a naked one whose body is the source unchanged. (As a function around a
//! block, LLVM cannot assemble the frame's unwind directives across a body
//! this long.)

const std = @import("std");
const arch = @import("__root__.zig");
const Error = arch.Error;

const Kernel = arch.Kernel;
const Pin = arch.Pin;

fn macro(name: []const u8) []const u8 {
    const v = arch.macros.get(name) orelse std.debug.panic("the source does not define {s}", .{name});
    return std.mem.trimStart(u8, v, "%");
}

fn pins(a: std.mem.Allocator, _: Kernel) Error![]const Pin {
    const out = try a.alloc(Pin, 6);
    out[0] = .{ .name = "m", .reg = macro("M"), .variable = "m", .init = "n" };
    out[1] = .{ .name = "x", .reg = macro("X"), .variable = "xp", .init = "x" };
    out[2] = .{ .name = "incx", .reg = macro("INCX"), .variable = "ix", .init = "inc_x" };
    out[3] = .{ .name = "y", .reg = macro("Y"), .variable = "yp", .init = "y" };
    out[4] = .{ .name = "incy", .reg = macro("INCY"), .variable = "iy", .init = "inc_y" };
    out[5] = .{ .name = "ap", .reg = "eax", .variable = "ap", .init = null };
    return out;
}

fn decls(a: std.mem.Allocator, k: Kernel) Error![]const []const u8 {
    const F = k.prec.float();
    const arr = if (k.prec.complex())
        try std.fmt.allocPrint(a, "const al = [2]{s}{{ alpha_r, alpha_i }};", .{F})
    else
        try std.fmt.allocPrint(a, "const al = [1]{s}{{alpha}};", .{F});
    const out = try a.alloc([]const u8, 2);
    out[0] = arr;
    out[1] = try std.fmt.allocPrint(a, "var ap: [*]const {s} = &al;", .{F});
    return out;
}

fn coreDefines(a: std.mem.Allocator, _: Kernel, core: []const u8) Error![]const []const u8 {
    const out = try a.alloc([]const u8, 1);
    out[0] = try std.fmt.allocPrint(a, "-D{s}", .{core});
    return out;
}

fn has(line: []const u8, needle: []const u8) bool {
    return std.mem.indexOf(u8, line, needle) != null;
}

fn transform(a: std.mem.Allocator, _: Kernel, lines: []const []const u8) Error![]const []const u8 {
    var out: std.ArrayList([]const u8) = .empty;
    for (lines) |line| {
        const t = std.mem.trim(u8, line, " ");
        // only what declares the function goes; the function is the source's own
        if (has(t, ".globl") or has(t, ".type k") or has(t, ".hidden") or std.mem.startsWith(u8, t, ".section .note")) continue;
        try out.append(a, line);
    }
    return out.items;
}

/// The 32-bit register a name belongs to.
fn canonical(reg: []const u8) ?[]const u8 {
    const map = [_][2][]const u8{
        .{ "eax", "eax" }, .{ "ax", "eax" }, .{ "al", "eax" },  .{ "ah", "eax" },
        .{ "ebx", "ebx" }, .{ "bx", "ebx" }, .{ "bl", "ebx" },  .{ "bh", "ebx" },
        .{ "ecx", "ecx" }, .{ "cx", "ecx" }, .{ "cl", "ecx" },  .{ "ch", "ecx" },
        .{ "edx", "edx" }, .{ "dx", "edx" }, .{ "dl", "edx" },  .{ "dh", "edx" },
        .{ "esi", "esi" }, .{ "si", "esi" }, .{ "edi", "edi" }, .{ "di", "edi" },
        .{ "ebp", "ebp" }, .{ "bp", "ebp" },
    };
    for (map) |m| if (std.mem.eql(u8, reg, m[0])) return m[1];
    return null;
}

fn clobbers(a: std.mem.Allocator, text: []const []const u8, pinned: []const Pin) Error![]const []const u8 {
    var seen: std.StringArrayHashMapUnmanaged(void) = .empty;
    // the body works on the x87 stack
    for ([_][]const u8{ "st0", "st1", "st2", "st3", "st4", "st5", "st6", "st7", "fpsr" }) |r| try seen.put(a, r, {});
    for (text) |line| {
        var i: usize = 0;
        while (std.mem.indexOfScalarPos(u8, line, i, '%')) |p| {
            var e = p + 1;
            while (e < line.len and std.ascii.isAlphanumeric(line[e])) e += 1;
            const name = canonical(line[p + 1 .. e]) orelse {
                i = e;
                continue;
            };
            i = e;
            const is_pin = for (pinned) |o| {
                if (std.mem.eql(u8, o.reg, name)) break true;
            } else false;
            if (!is_pin) try seen.put(a, name, {});
        }
    }
    return seen.keys();
}

pub const x86: arch.Arch = .{
    .name = "x86",
    .triple = "x86-linux-gnu",
    .config_h = "#define ARCH_X86 1\n#define BINARY32 1\n",
    .defines = &.{ "-DASSEMBLER", "-DARCH_X86", "-DOS_LINUX", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k" },
    .core_defines = coreDefines,
    .pins = pins,
    .decls = decls,
    .transform = transform,
    .clobbers = clobbers,
    .always_clobbers = &.{ "memory", "cc" },
    .naked = true,
};
