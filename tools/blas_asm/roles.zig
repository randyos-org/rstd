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
pub fn parse(a: std.mem.Allocator, text: []const u8) Error![]const Role {
    var out: std.ArrayList(Role) = .empty;
    var lines = std.mem.splitScalar(u8, text, '\n');
    while (lines.next()) |line| {
        if (!std.mem.startsWith(u8, line, "$(KDIR)")) continue;
        const colon = std.mem.indexOf(u8, line, " : $(KERNELDIR)/$(") orelse continue;
        const tail = line[colon + " : $(KERNELDIR)/$(".len ..];
        const close = std.mem.indexOfScalar(u8, tail, ')') orelse continue;
        const key = tail[0..close];

        const obj_start = "$(KDIR)".len;
        const obj_end = std.mem.indexOf(u8, line[obj_start..], "_k$(TSUFFIX)") orelse continue;
        const name = line[obj_start .. obj_start + obj_end];

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
