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
    /// The vtable role (`cdotu`); empty for the axpy kernels, whose role follows
    /// from precision and conjugation.
    role_name: []const u8 = "",
    /// `-D` symbols the role's rule adds.
    defines: []const []const u8 = &.{},

    /// `saxpy`, `zaxpyc`, ...
    pub fn role(k: Kernel, a: std.mem.Allocator) Error![]const u8 {
        if (k.role_name.len > 0) return k.role_name;
        return try std.fmt.allocPrint(a, "{t}axpy{s}", .{ k.prec, if (k.conj) "c" else "" });
    }
    /// The source without its extension.
    pub fn stem(k: Kernel) []const u8 {
        return k.source[0..std.mem.lastIndexOfScalar(u8, k.source, '.').?];
    }
};

/// The lines of a source as a naked function's body: the assembly of the source
/// as it is, less what declares the function. Used for every routine but axpy.
pub fn nakedLines(a: std.mem.Allocator, lines: []const []const u8) Error![]const []const u8 {
    var out: std.ArrayList([]const u8) = .empty;
    for (lines) |line| {
        // the prologue macro puts several statements on one line
        var parts = std.mem.splitScalar(u8, line, ';');
        while (parts.next()) |part| {
            const t = std.mem.trim(u8, part, " ");
            if (t.len == 0) continue;
            if (std.mem.startsWith(u8, t, ".globl") or std.mem.startsWith(u8, t, ".global") or
                std.mem.startsWith(u8, t, ".hidden") or std.mem.startsWith(u8, t, ".type") or
                std.mem.startsWith(u8, t, ".size") or std.mem.startsWith(u8, t, ".section .note") or
                std.mem.startsWith(u8, t, ".file") or std.mem.startsWith(u8, t, ".ident")) continue;
            try out.append(a, try std.fmt.allocPrint(a, " {s}", .{t}));
        }
    }
    return out.items;
}

/// A register the body expects an argument in, and the Zig variable that
/// carries it.
pub const Pin = struct {
    /// The operand's name.
    name: []const u8,
    /// The register, as the compiler's constraint spells it.
    reg: []const u8,
    /// The Zig variable passed for it.
    variable: []const u8,
    /// Copied from the parameter the variable is initialised from; null
    /// when the architecture declares the variable itself.
    init: ?[]const u8,
    /// Whether the body writes it. Read-only pins go in the inputs.
    inout: bool = true,
};

/// The macros the preprocessor ends with for the source being converted, as
/// `-dM` lists them: where a source names its argument registers (`#define X
/// r7`), an architecture reads them here instead of keeping a copy. Set by
/// `gen.zig` before it asks for the pins.
pub var macros: std.StringHashMapUnmanaged([]const u8) = .empty;

/// A way the same source is built more than once for one architecture, such
/// as 32-bit and 64-bit, each with its own file.
pub const Variant = struct {
    /// In the file name, and how the registry picks between them.
    name: []const u8,
    defines: []const []const u8,
};

pub const Arch = struct {
    /// The directory under OpenBLAS's `kernel/` and rstd's `blas1/axpy/`.
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
    /// The registers the body reads its arguments from.
    pins: *const fn (a: std.mem.Allocator, k: Kernel) Error![]const Pin,
    /// Zig statements that declare what the pins refer to, if not by `init`.
    decls: *const fn (a: std.mem.Allocator, k: Kernel) Error![]const []const u8,
    /// Turns the preprocessed, tab-free, non-empty lines of a source into the
    /// lines of the `asm` block: drops what declares the function, does the
    /// argument shuffling Zig now does, and ends every return at `end_label`.
    transform: *const fn (a: std.mem.Allocator, k: Kernel, lines: []const []const u8) Error![]const []const u8,
    /// The Zig clobber names for a body, excluding the pinned registers.
    clobbers: *const fn (a: std.mem.Allocator, lines: []const []const u8, pins: []const Pin) Error![]const []const u8,
    /// Sources the toolchain cannot assemble. They are not converted; the
    /// platforms and cores are listed in `src/linalg/blas/readme.md`.
    unsupported: []const []const u8 = &.{},
    /// The ways to build each source; empty name for just one.
    variants: []const Variant = &.{.{ .name = "", .defines = &.{} }},
    /// Clobbers every body has: memory, and the flag register under the
    /// compiler's name for it.
    always_clobbers: []const []const u8,
    /// Emit each function as a naked function with the source unchanged,
    /// instead of a function around a block with pinned operands.
    naked: bool = false,
    /// Whether this source is one to convert (some are C, done elsewhere).
    is_asm: *const fn (source: []const u8) bool = isDotS,
};

pub fn isDotS(source: []const u8) bool {
    return std.mem.endsWith(u8, source, ".S") and std.mem.indexOf(u8, source, "/") == null;
}

/// The label every return in a body jumps to.
pub const end_label = "9999";
