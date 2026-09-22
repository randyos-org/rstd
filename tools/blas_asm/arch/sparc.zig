//! SPARC (V9, 32- and 64-bit). Each source `save`s its own register window
//! and reads the arguments beyond the sixth from the caller's frame; here the
//! Zig function's window is the source's, so its `%i` registers already hold
//! what the source expects, the register names being read from the source's
//! `#define`s (`N`, `X`, `INCX`, `Y`, `INCY`) as for POWER, and the frame
//! loads and the `save` go. `alpha` arrives in a floating-point register the
//! compiler picks; it is passed by address in `%o1` and loaded into the
//! registers the source reads it from (`%f6`, `%f7`, `%f8`, `%f9`), with the
//! source's own copies into its working registers left as they are.
//!
//! The return is `return %i7 + 8`, which restores the window; the block ends
//! with a branch instead, its delay slot being the source's own.

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
    out[0] = .{ .name = "n", .reg = macro("N"), .variable = "m", .init = "n" };
    out[1] = .{ .name = "x", .reg = macro("X"), .variable = "xp", .init = "x" };
    out[2] = .{ .name = "incx", .reg = macro("INCX"), .variable = "ix", .init = "inc_x" };
    out[3] = .{ .name = "y", .reg = macro("Y"), .variable = "yp", .init = "y" };
    out[4] = .{ .name = "incy", .reg = macro("INCY"), .variable = "iy", .init = "inc_y" };
    out[5] = .{ .name = "ap", .reg = "o1", .variable = "ap", .init = null, .inout = false };
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
    out[1] = try std.fmt.allocPrint(a, "const ap: [*]const {s} = &al;", .{F});
    return out;
}

fn coreDefines(a: std.mem.Allocator, _: Kernel, core: []const u8) Error![]const []const u8 {
    const out = try a.alloc([]const u8, 2);
    out[0] = try std.fmt.allocPrint(a, "-D{s}", .{core});
    out[1] = "-DSPARC";
    return out;
}

fn has(line: []const u8, needle: []const u8) bool {
    return std.mem.indexOf(u8, line, needle) != null;
}

fn transform(a: std.mem.Allocator, k: Kernel, lines: []const []const u8) Error![]const []const u8 {
    var out: std.ArrayList([]const u8) = .empty;
    const size: usize = if (k.prec.double()) 8 else 4;
    const ld = if (k.prec.double()) "ldd" else "ld";
    var from_frame: usize = 0;
    var seen_frame_alpha = false;
    for (lines) |line| {
        if (has(line, "[%sp")) {
            const t = std.mem.trim(u8, line, " ");
            // a floating-point load from the frame: `alpha`
            if (std.mem.startsWith(u8, t, "ld") and has(t, "%f")) {
                const comma = std.mem.lastIndexOf(u8, t, ", ").?;
                try out.append(a, try std.fmt.allocPrint(a, " {s} [%o1 + {d}], {s}", .{ ld, from_frame * size, t[comma + 2 ..] }));
                from_frame += 1;
                seen_frame_alpha = true;
            }
            continue;
        }
        const t = std.mem.trim(u8, line, " ");
        // what declares the function, and the window
        if (has(t, ".global") or has(t, ".type k") or std.mem.startsWith(u8, t, ".size") or std.mem.startsWith(u8, t, ".section .note")) continue;
        if (std.mem.startsWith(u8, t, "save ")) continue;
        if (std.mem.startsWith(u8, t, "return %i7")) {
            try out.append(a, " ba " ++ arch.end_label ++ "f");
            continue;
        }
        try out.append(a, line);
    }
    if (!seen_frame_alpha) {
        // 64-bit: `alpha` was in the argument registers; the copies the source
        // makes from them read what is loaded here
        var pre: std.ArrayList([]const u8) = .empty;
        if (k.prec.complex()) {
            if (k.prec.double()) {
                try pre.append(a, " ldd [%o1], %f6");
                try pre.append(a, " ldd [%o1 + 8], %f8");
            } else {
                try pre.append(a, " ld [%o1], %f7");
                try pre.append(a, " ld [%o1 + 4], %f9");
            }
        } else if (k.prec.double()) {
            try pre.append(a, " ldd [%o1], %f6");
        } else {
            try pre.append(a, " ld [%o1], %f7");
        }
        try pre.appendSlice(a, out.items);
        out = pre;
    }
    try out.append(a, arch.end_label ++ ":");
    return out.items;
}

/// The register a name is, in the form the Zig clobber list has it: `%f62`
/// is in `q15`, `%i3` is `i3`.
fn canonical(a: std.mem.Allocator, reg: []const u8) !?[]const u8 {
    if (reg.len < 2) return null;
    const n = std.fmt.parseInt(u8, reg[1..], 10) catch return null;
    switch (reg[0]) {
        'g' => return if (n >= 1 and n <= 7) try std.fmt.allocPrint(a, "g{d}", .{n}) else null,
        'o' => return if (n <= 5) try std.fmt.allocPrint(a, "o{d}", .{n}) else null,
        'l' => return if (n <= 7) try std.fmt.allocPrint(a, "l{d}", .{n}) else null,
        'i' => return if (n <= 5) try std.fmt.allocPrint(a, "i{d}", .{n}) else null,
        'f' => return if (n <= 63) try std.fmt.allocPrint(a, "q{d}", .{n / 4}) else null,
        else => return null,
    }
}

fn clobbers(a: std.mem.Allocator, text: []const []const u8, pinned: []const Pin) Error![]const []const u8 {
    var seen: std.StringArrayHashMapUnmanaged(void) = .empty;
    for (text) |line| {
        var i: usize = 0;
        while (std.mem.indexOfScalarPos(u8, line, i, '%')) |p| {
            var e = p + 1;
            while (e < line.len and std.ascii.isAlphanumeric(line[e])) e += 1;
            const name = try canonical(a, line[p + 1 .. e]);
            i = e;
            if (name == null) continue;
            const is_pin = for (pinned) |o| {
                if (try canonical(a, o.reg)) |c| if (std.mem.eql(u8, c, name.?)) break true;
            } else false;
            if (!is_pin) try seen.put(a, name.?, {});
        }
    }
    return seen.keys();
}

pub const sparc: arch.Arch = .{
    .name = "sparc",
    .triple = "sparc64-linux-gnu",
    .config_h = "#define ARCH_SPARC 1\n#define BINARY64 1\n",
    .defines = &.{ "-DASSEMBLER", "-DARCH_SPARC", "-DOS_LINUX", "-DASMNAME=k", "-DASMFNAME=k_", "-DNAME=k_", "-DCNAME=k" },
    .core_defines = coreDefines,
    .pins = pins,
    .decls = decls,
    .transform = transform,
    .clobbers = clobbers,
    .always_clobbers = &.{ "memory", "icc", "xcc", "fsr" },
    .variants = &.{
        .{ .name = "32", .defines = &.{} },
        .{ .name = "64", .defines = &.{"-D__64BIT__"} },
    },
};
