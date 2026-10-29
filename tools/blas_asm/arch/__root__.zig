//! What `gen.zig` needs to know about one architecture. Each
//! `arch/<name>.zig` fills in an `Arch`; nothing here is specific to one.

const std = @import("std");

/// What the converters raise: the allocations they make are the only fallible calls.
pub const Error = std.mem.Allocator.Error;

/// Real or complex, single or double.
pub const Prec = enum {
    s,
    d,
    c,
    z,

    pub fn complex(p: Prec) bool {
        return p == .c or p == .z;
    }
    pub fn double(p: Prec) bool {
        return p == .d or p == .z;
    }
    /// The Zig type of one component.
    pub fn float(p: Prec) []const u8 {
        return if (p.double()) "f64" else "f32";
    }
    /// The letter `Makefile.L1` keys the source by.
    pub fn letter(p: Prec) []const u8 {
        return switch (p) {
            .s => "S",
            .d => "D",
            .c => "C",
            .z => "Z",
        };
    }
};

/// One source, built for one precision.
pub const Kernel = struct {
    /// Under `kernel/<arch>/` in the OpenBLAS checkout.
    source: []const u8,
    prec: Prec,
    /// The conjugated complex form.
    conj: bool = false,
    /// The vtable role (`cdotu`).
    role_name: []const u8 = "",
    /// Its vtable slot (`cdotu_k`).
    slot: []const u8 = "",
    /// `-D` symbols the role's rule adds.
    defines: []const []const u8 = &.{},

    /// `saxpy`, `zaxpyc`, ...
    pub fn role(k: Kernel, a: std.mem.Allocator) Error![]const u8 {
        if (k.role_name.len > 0) return k.role_name;
        return try a.print("{t}axpy{s}", .{ k.prec, if (k.conj) "c" else "" });
    }
    /// The source without its extension.
    pub fn stem(k: Kernel) []const u8 {
        return k.source[0..std.mem.findScalarLast(u8, k.source, '.').?];
    }
};

/// Whether `c` can be part of a plain assembly name (a coarser version of
/// the identically-named check in `gen.zig`, which this file cannot import).
fn isNameChar(c: u8) bool {
    return std.ascii.isAlphanumeric(c) or c == '_' or c == '.' or c == '$';
}

/// `-DDYNAMIC_ARCH` (`gen.zig` always sets it) turns some sources' calls into
/// another routine's kernel, upstream self-contained via that routine's own
/// `gotoblas`-style dispatch trampoline (e.g. `bl camax_k`), into a direct
/// call to one core's own kernel (`bl camax_k_LA264`): OpenBLAS's own name
/// for exactly the symbol this project's per-core kernel object already
/// builds and exports, under this project's own naming instead
/// (`util.kernelSymbol`: `blaskern_camax_k_la264`). Rewrites any such call to
/// that name, so it resolves within the same isolated per-core object.
fn rewriteDynamicArchCalls(a: std.mem.Allocator, line: []const u8) Error![]const u8 {
    if (std.mem.find(u8, line, "_k_") == null) return line;
    var out: std.ArrayList(u8) = .empty;
    var i: usize = 0;
    while (i < line.len) {
        if (isNameChar(line[i]) and (i == 0 or !isNameChar(line[i - 1]))) {
            var e = i;
            while (e < line.len and isNameChar(line[e])) e += 1;
            const tok = line[i..e];
            if (std.mem.find(u8, tok, "_k_")) |sep| rewrite: {
                const role = tok[0..sep];
                const core = tok[sep + "_k_".len ..];
                if (role.len == 0 or core.len == 0) break :rewrite;
                // digits are neither upper nor lower, so they pass either way
                for (role) |c| if (std.ascii.isUpper(c)) break :rewrite;
                for (core) |c| if (std.ascii.isLower(c)) break :rewrite;
                try out.print(a, "blaskern_{s}_k_{s}", .{ role, try std.ascii.allocLowerString(a, core) });
                i = e;
                continue;
            }
            try out.appendSlice(a, tok);
            i = e;
        } else {
            try out.append(a, line[i]);
            i += 1;
        }
    }
    return out.items;
}

/// The lines of a source as a naked function's body: the assembly of the source
/// as it is, less what declares the function, with `.align` as `.balign` where
/// `ar.align_bytes`.
pub fn nakedLines(a: std.mem.Allocator, ar: Arch, lines: []const []const u8) Error![]const []const u8 {
    var out: std.ArrayList([]const u8) = .empty;
    for (lines) |line| {
        // the prologue macro puts several statements on one line
        var parts = std.mem.splitScalar(u8, line, ';');
        while (parts.next()) |part| {
            const t = std.mem.trim(u8, part, " \t");
            if (t.len == 0) continue;
            if (std.mem.startsWith(u8, t, ".globl") or std.mem.startsWith(u8, t, ".global") or
                std.mem.startsWith(u8, t, ".hidden") or std.mem.startsWith(u8, t, ".type") or
                std.mem.startsWith(u8, t, ".size") or std.mem.startsWith(u8, t, ".section .note") or
                std.mem.startsWith(u8, t, ".file") or std.mem.startsWith(u8, t, ".ident") or
                // MIPS/Alpha/IA64-style PROLOGUE/EPILOGUE: `.ent NAME`/`.end NAME`
                // (not `.endm`/`.endif`/`.endr`, which need the space to miss).
                std.mem.startsWith(u8, t, ".ent ") or std.mem.eql(u8, t, ".end") or
                std.mem.startsWith(u8, t, ".end ")) continue;
            // Zig places the function in its own section; the source's `.text`
            // would move the code out of it where functions get a section each.
            // Any other section would hold data the code needs: not handled.
            if (std.mem.eql(u8, t, ".text")) continue;
            if (std.mem.startsWith(u8, t, ".section ")) {
                const name = std.mem.trim(u8, t[".section ".len..], " \t\"");
                if (std.mem.eql(u8, name, ".text")) continue;
                std.debug.panic("a naked function cannot switch sections: {s}", .{t});
            }
            if (ar.align_bytes and std.mem.startsWith(u8, t, ".align ")) {
                try out.append(a, try a.print(" .b{s}", .{t[1..]}));
                continue;
            }
            try out.append(a, try a.print(" {s}", .{try rewriteDynamicArchCalls(a, t)}));
        }
    }
    return out.items;
}

/// A way the same source is built more than once for one architecture, such
/// as 32-bit and 64-bit, each with its own file.
pub const Variant = struct {
    /// In the file name, and how the registry picks between them.
    name: []const u8,
    defines: []const []const u8,
    /// The target to preprocess for, where it is not `Arch.triple`.
    triple: ?[]const u8 = null,
    /// Fixes for sources built this way only, beyond `Arch.patches`.
    patches: []const Patch = &.{},
    /// A change to every statement of every source built this way, after
    /// preprocessing (`nakedLines`'s output, one statement each).
    rewrite: ?*const fn (a: std.mem.Allocator, statement: []const u8) Error![]const u8 = null,
    /// Cores (`KERNEL.<core>` names, any case) never built this way, such as
    /// 32-bit-only cores in a 64-bit variant (`cpu/_corespec.zig`'s `runsOn`).
    without: []const []const u8 = &.{},

    /// Whether `core` has a build of this variant.
    pub fn builds(v: Variant, core: []const u8) bool {
        for (v.without) |c| if (std.ascii.eqlIgnoreCase(c, core)) return false;
        return true;
    }
};

/// Our fix for a bug in an upstream source, applied to its text before it is
/// preprocessed. `find` must occur exactly `count` times, so a fix that no
/// longer matches the source stops the generator rather than vanishing.
pub const Patch = struct {
    /// The file under `kernel/<arch>/`.
    source: []const u8,
    find: []const u8,
    replace: []const u8,
    count: usize = 1,
};

pub const Arch = struct {
    /// The directory under OpenBLAS's `kernel/` and under each `blas1/<routine>/`.
    name: []const u8,
    /// The target the preprocessor is told, to pick the right headers.
    triple: []const u8,
    /// A stand-in for the `config.h` the OpenBLAS build generates.
    config_h: []const u8,
    /// `-D` symbols for every source.
    defines: []const []const u8,
    /// `-D` symbols for one kernel on one core, beyond `DOUBLE`, `COMPLEX`
    /// and `CONJ`. `core` is the upper-case `KERNEL.<core>` name.
    core_defines: *const fn (a: std.mem.Allocator, k: Kernel, core: []const u8) Error![]const []const u8,
    /// Sources the toolchain cannot assemble. They are not converted; the
    /// platforms and cores are listed in `src/linalg/blas/readme.md`.
    unsupported: []const []const u8 = &.{},
    /// Fixes for upstream bugs, for every variant.
    patches: []const Patch = &.{},
    /// Fixes to sources no longer converted (`roles.skip`), never applied:
    /// the record of the fixes `issues.md` words for upstream. Each run of
    /// the architecture checks they still match.
    retired: []const Patch = &.{},
    /// The ways to build each source; empty name for just one.
    variants: []const Variant = &.{.{ .name = "", .defines = &.{} }},
    /// Whether this source is one to convert (some are C, done elsewhere).
    is_asm: *const fn (source: []const u8) bool = isDotS,
    /// Whether the sources' `.align` takes a byte count, as x86's does for
    /// ELF and COFF. Mach-O reads it as a power of two instead, so these are
    /// written as `.balign`, which means bytes on every object format.
    align_bytes: bool = false,
    /// Stops the generator on a statement some target of this architecture
    /// cannot take (`nakedLines`'s output, one statement each).
    check: ?*const fn (source: []const u8, statement: []const u8) void = null,
};

pub fn isDotS(source: []const u8) bool {
    return std.mem.endsWith(u8, source, ".S") and std.mem.find(u8, source, "/") == null;
}
