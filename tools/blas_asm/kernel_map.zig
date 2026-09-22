//! Reads OpenBLAS's kernel selection: which source each core uses for each
//! `axpy` precision, as `kernel/<arch>/KERNEL.<core>`, then `KERNEL`, then the
//! defaults in `kernel/Makefile.L1` decide it.

const std = @import("std");
const Io = std.Io;

pub var arena: std.mem.Allocator = undefined;
pub var io: Io = undefined;
/// The OpenBLAS checkout.
pub var root: []const u8 = undefined;
/// Where the converted files go, under `<out_root>/<arch>/asm/`.
pub var out_root: []const u8 = "src/linalg/blas/blas1/axpy";

pub fn die(comptime f: []const u8, args: anytype) noreturn {
    std.debug.print(f ++ "\n", args);
    std.process.exit(1);
}

pub fn init(a: std.mem.Allocator, i: Io, checkout: []const u8) !void {
    arena = a;
    io = i;
    root = checkout;
    if (read(try path(&.{ root, "kernel", "Makefile.L1" })) == null) die("cannot read kernel/Makefile.L1", .{});
}

pub fn path(parts: []const []const u8) ![]const u8 {
    return try std.mem.join(arena, "/", parts);
}

pub fn read(p: []const u8) ?[]const u8 {
    return Io.Dir.cwd().readFileAlloc(io, p, arena, .limited(1 << 22)) catch null;
}

const Values = std.ArrayList([]const u8);

/// The alternatives one variable has been given so far.
const Assigns = struct {
    /// Set by a core's own file, before the defaults.
    values: Values = .empty,

    fn set(self: *Assigns, v: []const u8, conditional: bool) !void {
        if (!conditional) self.values.clearRetainingCapacity();
        for (self.values.items) |o| if (std.mem.eql(u8, o, v)) return;
        try self.values.append(arena, v);
    }
};

/// Applies the assignments to `<f>AXPYKERNEL` in `path`, in order, following
/// `include`s where they appear as make does. An unconditional assignment
/// replaces what was there; one under a conditional adds an alternative. With
/// `defaults` set, only a variable nothing has set yet is assigned, which is
/// what the `ifndef` blocks of the default files mean.
fn apply(a: *Assigns, file: []const u8, kernel_dir: []const u8, key: []const u8, defaults: bool) !void {
    const text = read(file) orelse return;
    var depth: usize = 0;
    var lines = std.mem.splitScalar(u8, text, '\n');
    while (lines.next()) |raw| {
        const line = std.mem.trim(u8, raw, " \t\r");
        if (std.mem.startsWith(u8, line, "ifeq") or std.mem.startsWith(u8, line, "ifneq") or
            std.mem.startsWith(u8, line, "ifdef") or std.mem.startsWith(u8, line, "ifndef"))
        {
            depth += 1;
            continue;
        }
        if (std.mem.startsWith(u8, line, "endif")) {
            depth -|= 1;
            continue;
        }
        if (std.mem.startsWith(u8, line, "include ") or std.mem.startsWith(u8, line, "-include ")) {
            const target = std.mem.trim(u8, line[std.mem.indexOfScalar(u8, line, ' ').?..], " ");
            const marker = "$(KERNELDIR)";
            if (std.mem.startsWith(u8, target, marker)) {
                const rest = target[marker.len..];
                try apply(a, try std.mem.concat(arena, u8, &.{ kernel_dir, rest }), kernel_dir, key, defaults);
            }
            continue;
        }
        if (!std.mem.startsWith(u8, line, key)) continue;
        var rest = std.mem.trimStart(u8, line[key.len..], " \t");
        if (rest.len == 0 or rest[0] != '=') continue;
        rest = std.mem.trim(u8, rest[1..], " \t");
        if (defaults) {
            if (a.values.items.len == 0) try a.set(rest, false);
        } else try a.set(rest, depth > 0);
    }
}

/// The cores of `arch`: the suffixes of its `KERNEL.<core>` files, sorted,
/// or `(all)` when it has only `KERNEL`.
pub fn cores(arch: []const u8) !?[]const []const u8 {
    var dir = Io.Dir.cwd().openDir(io, try path(&.{ root, "kernel", arch }), .{ .iterate = true }) catch return null;
    defer dir.close(io);
    var out: std.ArrayList([]const u8) = .empty;
    var it = dir.iterate();
    while (it.next(io) catch null) |e| {
        if (e.kind != .file or !std.mem.startsWith(u8, e.name, "KERNEL.")) continue;
        try out.append(arena, try arena.dupe(u8, e.name["KERNEL.".len..]));
    }
    std.mem.sort([]const u8, out.items, {}, struct {
        fn lt(_: void, a: []const u8, b: []const u8) bool {
            return std.ascii.lessThanIgnoreCase(a, b);
        }
    }.lt);
    if (out.items.len == 0) try out.append(arena, "(all)");
    return out.items;
}

/// The sources `core` uses for precision `f` (`S`, `D`, `C`, `Z`, `Q`, `X`);
/// several when a conditional picks between them.
pub fn resolve(arch: []const u8, core: []const u8, f: []const u8) ![]const []const u8 {
    const key = try std.fmt.allocPrint(arena, "{s}AXPYKERNEL", .{f});
    return try resolveKey(arch, core, key, "Makefile.L1");
}

/// The sources `core` assigns to the make variable `key`, with the defaults
/// of `kernel/<defaults>` when neither its `KERNEL.<core>` nor `KERNEL` does.
pub fn resolveKey(arch: []const u8, core: []const u8, key: []const u8, defaults: []const u8) ![]const []const u8 {
    const dir = try path(&.{ root, "kernel", arch });
    var a: Assigns = .{};
    if (!std.mem.eql(u8, core, "(all)")) {
        try apply(&a, try path(&.{ dir, try std.fmt.allocPrint(arena, "KERNEL.{s}", .{core}) }), dir, key, false);
    }
    try apply(&a, try path(&.{ dir, "KERNEL" }), dir, key, true);
    if (a.values.items.len == 0) try apply(&a, try path(&.{ root, "kernel", defaults }), dir, key, true);
    return a.values.items;
}
