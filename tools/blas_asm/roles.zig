//! The kernels a `kernel/Makefile.L*` builds, read from its rules:
//!
//!     $(KDIR)cdotu_k$(TSUFFIX).$(SUFFIX) ... : $(KERNELDIR)/$(CDOTKERNEL)
//!             $(CC) -c $(CFLAGS) -DCOMPLEX -UCONJ $< -o $@
//!
//! Each rule is one role (`cdotu`, the vtable's `cdotu_k`), the make variable
//! its source comes from, and the flags that build it: which of precision,
//! complex and conjugated it is, and any other `-D` symbol.

const std = @import("std");
const arch = @import("arch/__root__.zig");
const Error = arch.Error;

pub const Role = struct {
    /// `cdotu`.
    name: []const u8,
    /// `CDOTKERNEL`.
    key: []const u8,
    prec: arch.Prec,
    conj: bool,
    /// The `-D` symbols beyond `DOUBLE`, `COMPLEX` and `CONJ`.
    defines: []const []const u8,
};

/// The roles of `text`, a Makefile, skipping the extended-precision ones
/// (`XDOUBLE`), which no supported target builds.
///
/// The scan does not understand `ifdef`/`else`/`endif`: where a role's rule
/// is guarded by one (`dsdot`/`sdsdot`, built from `DSDOTKERNEL` if set, else
/// from `SDOTKERNEL`), both branches parse as separate roles of the same
/// name. Only the first (the `ifdef`-true branch, textually first) is kept,
/// matching Make's own precedence; a later same-named role is dropped.
pub fn parse(a: std.mem.Allocator, text: []const u8) Error![]const Role {
    var out: std.ArrayList(Role) = .empty;
    var lines = std.mem.splitScalar(u8, text, '\n');
    outer: while (lines.next()) |line| {
        if (!std.mem.startsWith(u8, line, "$(KDIR)")) continue;
        const colon = std.mem.indexOf(u8, line, " : $(KERNELDIR)/$(") orelse continue;
        const tail = line[colon + " : $(KERNELDIR)/$(".len ..];
        const close = std.mem.indexOfScalar(u8, tail, ')') orelse continue;
        const key = tail[0..close];

        const obj_start = "$(KDIR)".len;
        const obj_end = std.mem.indexOf(u8, line[obj_start..], "_k$(TSUFFIX)") orelse continue;
        const name = line[obj_start .. obj_start + obj_end];
        for (out.items) |r| {
            if (std.mem.eql(u8, r.name, name)) continue :outer;
        }

        // the recipe line
        const flags = lines.next() orelse break;
        if (std.mem.indexOf(u8, flags, "-DXDOUBLE") != null) continue;
        var double = false;
        var complex = false;
        var conj = false;
        var defines: std.ArrayList([]const u8) = .empty;
        var it = std.mem.tokenizeAny(u8, flags, " \t");
        while (it.next()) |tok| {
            if (std.mem.eql(u8, tok, "-DDOUBLE")) double = true else if (std.mem.eql(u8, tok, "-DCOMPLEX")) complex = true else if (std.mem.eql(u8, tok, "-DCONJ")) conj = true else if (std.mem.startsWith(u8, tok, "-D") and !std.mem.startsWith(u8, tok, "-DXDOUBLE")) try defines.append(a, try a.dupe(u8, tok));
        }
        const prec: arch.Prec = if (complex) (if (double) .z else .c) else (if (double) .d else .s);
        try out.append(a, .{ .name = try a.dupe(u8, name), .key = try a.dupe(u8, key), .prec = prec, .conj = conj, .defines = defines.items });
    }
    return out.items;
}

/// Keys that don't reduce to `<precision><base>KERNEL`, so `dirFor` cannot
/// strip them algorithmically: the bf16 dot forms (built from `sdot.c`/`ddot.c`
/// under different flags, so they keep the plain dot kernel's variable name
/// instead of one prefixed with their own precision) and the bf16 conversions.
const dir_overrides = [_]struct { key: []const u8, dir: []const u8 }{
    .{ .key = "DSDOTKERNEL", .dir = "dot" },
    .{ .key = "SBDOTKERNEL", .dir = "dot" },
    .{ .key = "BF16TOKERNEL", .dir = "bf16" },
    .{ .key = "TOBF16KERNEL", .dir = "bf16" },
};

/// The `blas1/<dir>` a role's kernel lives under, from its `key`: one BLAS
/// routine family per directory, precisions and the conjugated/index forms
/// together, the way `blas1/axpy/` already holds every `axpy` kernel. Regular
/// keys are `<precision><base>KERNEL`, or `I<precision><base>KERNEL` for an
/// index routine (kept as an `i` prefix on the directory, so `iamax` and
/// `amax` stay apart, as the routines themselves do).
pub fn dirFor(a: std.mem.Allocator, key: []const u8) ![]const u8 {
    for (dir_overrides) |o| {
        if (std.mem.eql(u8, o.key, key)) return o.dir;
    }
    var s = key;
    var prefix: []const u8 = "";
    if (s.len > 0 and s[0] == 'I') {
        prefix = "i";
        s = s[1..];
    }
    if (s.len > 0 and std.mem.indexOfScalar(u8, "SDCZQXB", s[0]) != null) s = s[1..];
    if (std.mem.endsWith(u8, s, "KERNEL")) s = s[0 .. s.len - "KERNEL".len];
    const lower = try std.ascii.allocLowerString(a, s);
    return std.fmt.allocPrint(a, "{s}{s}", .{ prefix, lower });
}
