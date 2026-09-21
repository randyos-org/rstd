//! Reads OpenBLAS's kernel selection: which source each core uses for each
//! `axpy` precision, as `kernel/<arch>/KERNEL.<core>`, then `KERNEL`, then the
//! defaults in `kernel/Makefile.L1` decide it.

const std = @import("std");
const Io = std.Io;

pub var arena: std.mem.Allocator = undefined;
pub var io: Io = undefined;
/// The OpenBLAS checkout.
pub var root: []const u8 = undefined;

var defaults: []const u8 = undefined;

pub fn die(comptime f: []const u8, args: anytype) noreturn {
    std.debug.print(f ++ "\n", args);
    std.process.exit(1);
}

pub fn init(a: std.mem.Allocator, i: Io, checkout: []const u8) void {
    arena = a;
    io = i;
    root = checkout;
    defaults = read(path(&.{ root, "kernel", "Makefile.L1" })) orelse
        die("cannot read kernel/Makefile.L1", .{});
}

pub fn path(parts: []const []const u8) []const u8 {
    return std.mem.join(arena, "/", parts) catch @panic("OOM");
}

pub fn read(p: []const u8) ?[]const u8 {
    return Io.Dir.cwd().readFileAlloc(io, p, arena, .limited(1 << 22)) catch null;
}

/// Every value `<f>AXPYKERNEL` is assigned in `text`, in order, without repeats.
fn assigned(text: []const u8, f: []const u8) []const []const u8 {
    var out: std.ArrayList([]const u8) = .empty;
    const key = std.fmt.allocPrint(arena, "{s}AXPYKERNEL", .{f}) catch @panic("OOM");
    var lines = std.mem.splitScalar(u8, text, '\n');
    while (lines.next()) |raw| {
        const line = std.mem.trim(u8, raw, " \t\r");
        if (!std.mem.startsWith(u8, line, key)) continue;
        var rest = std.mem.trimStart(u8, line[key.len..], " \t");
        if (rest.len == 0 or rest[0] != '=') continue;
        rest = std.mem.trim(u8, rest[1..], " \t");
        const dup = for (out.items) |o| {
            if (std.mem.eql(u8, o, rest)) break true;
        } else false;
        if (!dup) out.append(arena, rest) catch @panic("OOM");
    }
    return out.items;
}

/// The cores of `arch`: the suffixes of its `KERNEL.<core>` files, sorted,
/// or `(all)` when it has only `KERNEL`.
pub fn cores(arch: []const u8) ?[]const []const u8 {
    var dir = Io.Dir.cwd().openDir(io, path(&.{ root, "kernel", arch }), .{ .iterate = true }) catch return null;
    defer dir.close(io);
    var out: std.ArrayList([]const u8) = .empty;
    var it = dir.iterate();
    while (it.next(io) catch null) |e| {
        if (e.kind != .file or !std.mem.startsWith(u8, e.name, "KERNEL.")) continue;
        out.append(arena, arena.dupe(u8, e.name["KERNEL.".len..]) catch @panic("OOM")) catch @panic("OOM");
    }
    std.mem.sort([]const u8, out.items, {}, struct {
        fn lt(_: void, a: []const u8, b: []const u8) bool {
            return std.ascii.lessThanIgnoreCase(a, b);
        }
    }.lt);
    if (out.items.len == 0) out.append(arena, "(all)") catch @panic("OOM");
    return out.items;
}

/// The sources `core` uses for precision `f` (`S`, `D`, `C`, `Z`, `Q`, `X`);
/// several when a conditional picks between them.
pub fn resolve(arch: []const u8, core: []const u8, f: []const u8) []const []const u8 {
    const dir = path(&.{ root, "kernel", arch });
    if (!std.mem.eql(u8, core, "(all)")) {
        const own = read(path(&.{ dir, std.fmt.allocPrint(arena, "KERNEL.{s}", .{core}) catch @panic("OOM") }));
        if (own) |t| {
            const v = assigned(t, f);
            if (v.len > 0) return v;
        }
    }
    if (read(path(&.{ dir, "KERNEL" }))) |t| {
        const v = assigned(t, f);
        if (v.len > 0) return v;
    }
    return assigned(defaults, f);
}
