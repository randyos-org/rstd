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
    /// The vtable slot, as the rule's object names it (`cdotu_k`, `sgemv_n`);
    /// null for one `gotoblas_t` has none for (`no_slot`).
    slot: ?[]const u8,
    /// `CDOTKERNEL`.
    key: []const u8,
    prec: arch.Prec,
    conj: bool,
    /// The `-D` symbols beyond `DOUBLE`, `COMPLEX` and `CONJ`.
    defines: []const []const u8,
    /// An x87 `q`/`x` rule (`XDOUBLE`), never converted: the portable kernels
    /// at `XDouble` serve those slots. `prec` is then that of the same flags
    /// without it.
    extended: bool = false,
};

/// Rules whose kernel has no vtable slot, so nothing would install it:
/// `sdsdot_k`, which OpenBLAS compiles per core as it does `dsdot_k`, though
/// its own `gotoblas_t` has no such slot either.
const no_slot = [_][]const u8{"sdsdot_k"};

/// A kernel left unconverted (`skip`): its slot, and the architecture (null
/// for every one).
pub const Skip = struct {
    slot: []const u8,
    arch: ?[]const u8 = null,
    /// Whether the slot is routed only on the cores that build this assembly
    /// (x86_64's `usesVector`): the index lists them (`skipped`), and the C
    /// kernels of the other cores stay (`uk.zig` converts them).
    cores: bool = false,
};

/// Kernels not converted, since no core would install them: the library's
/// `vector.routed` sends their slots to the portable or `@Vector` kernels
/// (keep the two in step). The extremes, `nrm2` and `sdot` everywhere; on
/// x86_64 the reductions whose assembly peels to an aligned operand
/// (`vector.peels`), routed only on the cores that build it.
pub const skip = blk: {
    var out: []const Skip = &.{};
    for ([_][]const u8{
        "samax_k",  "damax_k",  "camax_k",  "zamax_k",
        "samin_k",  "damin_k",  "camin_k",  "zamin_k",
        "smax_k",   "dmax_k",   "smin_k",   "dmin_k",
        "isamax_k", "idamax_k", "icamax_k", "izamax_k",
        "isamin_k", "idamin_k", "icamin_k", "izamin_k",
        "ismax_k",  "idmax_k",  "ismin_k",  "idmin_k",
        "snrm2_k",  "dnrm2_k",  "cnrm2_k",  "znrm2_k",
        "sdot_k",
    }) |s| out = out ++ [_]Skip{.{ .slot = s }};
    for ([_][]const u8{
        "sasum_k", "dasum_k", "casum_k", "zasum_k", "csum_k",  "zsum_k",
        "ddot_k",  "cdotu_k", "cdotc_k", "zdotu_k", "zdotc_k", "dgemv_n",
    }) |s| out = out ++ [_]Skip{.{ .slot = s, .arch = "x86_64", .cores = true }};
    break :blk out;
};

/// The `skip` entry for `slot` on `arch_name`, if any.
pub fn skipped(arch_name: []const u8, slot: []const u8) ?Skip {
    for (skip) |s| {
        if (!std.mem.eql(u8, s.slot, slot)) continue;
        if (s.arch) |sa| if (!std.mem.eql(u8, sa, arch_name)) continue;
        return s;
    }
    return null;
}

/// The roles of `text`, a Makefile; the extended-precision ones (`XDOUBLE`)
/// are marked `extended`.
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
        const colon = std.mem.find(u8, line, " : $(KERNELDIR)/$(") orelse continue;
        const tail = line[colon + " : $(KERNELDIR)/$(".len ..];
        const close = std.mem.findScalar(u8, tail, ')') orelse continue;
        const key = tail[0..close];

        // `sdot_k$(TSUFFIX)` is role `sdot`; level 2 has `sgemv_n$(TSUFFIX)`
        const obj_start = "$(KDIR)".len;
        const obj_end = std.mem.find(u8, line[obj_start..], "$(TSUFFIX)") orelse continue;
        const object = line[obj_start .. obj_start + obj_end];
        var name = object;
        if (std.mem.endsWith(u8, name, "_k")) name = name[0 .. name.len - 2];
        const slotted = for (no_slot) |n| {
            if (std.mem.eql(u8, n, object)) break false;
        } else true;
        for (out.items) |r| {
            if (std.mem.eql(u8, r.name, name)) continue :outer;
        }

        // the recipe line
        const flags = lines.next() orelse break;
        const ext = std.mem.find(u8, flags, "-DXDOUBLE") != null;
        var double = false;
        var complex = false;
        var conj = false;
        var defines: std.ArrayList([]const u8) = .empty;
        var it = std.mem.tokenizeAny(u8, flags, " \t");
        while (it.next()) |tok| {
            if (std.mem.eql(u8, tok, "-DDOUBLE")) double = true else if (std.mem.eql(u8, tok, "-DCOMPLEX")) complex = true else if (std.mem.eql(u8, tok, "-DCONJ")) conj = true else if (std.mem.startsWith(u8, tok, "-D") and !std.mem.startsWith(u8, tok, "-DXDOUBLE")) try defines.append(a, try a.dupe(u8, tok));
        }
        const prec: arch.Prec = if (complex) (if (double) .z else .c) else (if (double) .d else .s);
        try out.append(a, .{ .name = try a.dupe(u8, name), .slot = if (slotted) try a.dupe(u8, object) else null, .key = try a.dupe(u8, key), .prec = prec, .conj = conj, .defines = defines.items, .extended = ext });
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

/// The level-2 routines, as their make variables spell them.
const level2_dirs = [_][2][]const u8{
    .{ "GEMV", "gemv" },
    .{ "GER", "ger" },
    .{ "SYMV", "symv" },
    .{ "HEMV", "hemv" },
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
    // level 2 keys carry a form after the routine (`SGEMVNKERNEL`,
    // `CGERUKERNEL`, `ZHEMV_V_KERNEL`) and the half-precision gemvs a
    // two-letter prefix (`SBGEMVNKERNEL`): the routine is found by name
    for (level2_dirs) |d| {
        if (std.mem.find(u8, key, d[0]) != null) return d[1];
    }
    var s = key;
    var prefix: []const u8 = "";
    if (s.len > 0 and s[0] == 'I') {
        prefix = "i";
        s = s[1..];
    }
    if (s.len > 0 and std.mem.findScalar(u8, "SDCZQXB", s[0]) != null) s = s[1..];
    if (std.mem.endsWith(u8, s, "KERNEL")) s = s[0 .. s.len - "KERNEL".len];
    const lower = try std.ascii.allocLowerString(a, s);
    return a.print("{s}{s}", .{ prefix, lower });
}

test "extended-precision rules parse, marked extended" {
    var arena: std.heap.ArenaAllocator = .init(std.testing.allocator);
    defer arena.deinit();
    const text =
        "$(KDIR)qdot_k$(TSUFFIX).$(SUFFIX) $(KDIR)qdot_k$(TPSUFFIX).$(PSUFFIX) : $(KERNELDIR)/$(QDOTKERNEL)\n" ++
        "\t$(CC) -c $(CFLAGS) -UCOMPLEX -DXDOUBLE $< -o $@\n" ++
        "$(KDIR)xdotc_k$(TSUFFIX).$(SUFFIX) $(KDIR)xdotc_k$(TPSUFFIX).$(PSUFFIX) : $(KERNELDIR)/$(XDOTKERNEL)\n" ++
        "\t$(CC) -c $(CFLAGS) -DCOMPLEX -DXDOUBLE -DCONJ $< -o $@\n" ++
        "$(KDIR)ddot_k$(TSUFFIX).$(SUFFIX) $(KDIR)ddot_k$(TPSUFFIX).$(PSUFFIX) : $(KERNELDIR)/$(DDOTKERNEL)\n" ++
        "\t$(CC) -c $(CFLAGS) -UCOMPLEX -DDOUBLE $< -o $@\n";
    const rs = try parse(arena.allocator(), text);
    try std.testing.expectEqual(@as(usize, 3), rs.len);
    try std.testing.expect(rs[0].extended and rs[1].extended and !rs[2].extended);
    try std.testing.expect(rs[1].conj and rs[1].defines.len == 0 and rs[1].prec.complex());
    try std.testing.expectEqual(arch.Prec.d, rs[2].prec);
}
