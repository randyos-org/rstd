//! The library tables: columns, headers, rows and the contestants of a race.

const std = @import("std");
const rstd = @import("rstd");
const bench = @import("../bench.zig");
const merge = @import("../merge.zig");
const bb = @import("../blas_bench.zig");
const ext = @import("../ext.zig");
const lb = @import("../libs.zig");
const alone = @import("../alone.zig");

const Io = std.Io;
const blas = rstd.math.kernels.blas;
const arm = bb.arm;
const ob_asm = bb.ob_asm;
const ext_cols = ext.ext_cols;
const Src = lb.Src;
const loaded = lb.loaded;

/// A library table's column: a contestant, and its library if it is one.
pub const Col = struct {
    label: []const u8,
    short: []const u8,
    lib: ?Src = null,
    /// Whether this build has it, for a column that is no library.
    on: bool = true,
    /// Its rows of `n` up to this only (a slow baseline).
    max_n: usize = std.math.maxInt(usize),

    /// rstd's own columns (rstd, `asm`, `vec x4`, `ob .S`) are out in `alone`.
    pub fn shown(c: Col) bool {
        return if (c.lib) |l| loaded(l) else c.on and !bb.alone;
    }
};

/// `lib1` and `lib2`'s contestants. Ratios are over the first, OpenBLAS.
pub const lib_cols = [_]Col{
    .{ .label = "openblas", .short = "ob", .lib = .openblas },
    .{ .label = "rstd", .short = "rstd" },
    .{ .label = "blis", .short = "blis", .lib = .blis },
    .{ .label = "accelerate", .short = "accel", .lib = .accelerate },
    .{ .label = "asm", .short = "asm", .on = arm },
    .{ .label = "vec x4", .short = "x4", .on = arm },
};
/// `lib3`'s: `ob .S` is OpenBLAS's vortex `.S` through the bench's loops.
pub const lib3_cols = lib_cols[0..4].* ++ [_]Col{.{ .label = "ob .S", .short = "ob .S", .on = ob_asm }};
pub const rstd_col = 1;

/// Column `o` of a race in the order checked: rstd's first, whose output
/// the libraries' are compared with.
pub fn rstdFirst(o: usize) usize {
    return switch (o) {
        0 => rstd_col,
        rstd_col => 0,
        else => o,
    };
}

/// rstd's largest difference from a library's result, in eps of the type
/// relative to the size of the terms, and what accumulation order explains.
pub const Diff = struct { eps: f64, bound: f64 };

/// What accumulation order explains for a sum of `n` terms: a few ulps
/// times `n`.
pub fn sumBound(n: usize) f64 {
    return 2 * (@as(f64, @floatFromInt(n)) + 8);
}

/// A table's header: `lead` and `rule`, then each shown contestant's rate,
/// its ratio over OpenBLAS, rstd's difference from each library, rstd's
/// error, and a hash of rstd's output.
pub fn printHeader(comptime lead: []const u8, comptime rule: []const u8, cols: []const Col) void {
    var n: usize = 0;
    header_lead = lead;
    std.debug.print("| " ++ lead ++ " |", .{});
    for (cols) |c| if (c.shown()) {
        std.debug.print(" {s} |", .{c.label});
        n += 1;
    };
    for (cols[1..]) |c| if (c.shown()) {
        std.debug.print(" {s}/ob |", .{c.short});
        n += 1;
    };
    if (fair()) {
        std.debug.print(" ob alone | rstd/ob alone |", .{});
        n += 2;
    }
    if (bb.alone) {
        std.debug.print("\n| " ++ rule ++ " |", .{});
        for (0..n) |_| std.debug.print(" -: |", .{});
        std.debug.print("\n", .{});
        return;
    }
    for (cols) |c| if (c.lib != null and c.shown()) {
        std.debug.print(" rstd-{s} |", .{c.short});
        n += 1;
    };
    std.debug.print(" rstd err | rstd bits |\n| " ++ rule ++ " |", .{});
    for (0..n + 1) |_| std.debug.print(" -: |", .{});
    std.debug.print(" --- |\n", .{});
}

/// The column labels a merge (`merge.zig`) reads as rates, and each one's
/// short name in a ratio's header (`x/ob`).
pub const merge_labels = blk: {
    var l: []const merge.Label = &.{};
    for (lib_cols ++ ext_cols[4..] ++ lib3_cols[4..]) |c| l = l ++ [_]merge.Label{.{ .label = c.label, .short = c.short }};
    l = l ++ [_]merge.Label{.{ .label = "ob alone", .short = "ob alone" }};
    // the `l1`, `sparse`, `inline` and `compensated` tables' (#567)
    for ([_][2][]const u8{
        .{ "asm GB/s", "asm" },    .{ "core", "core" },     .{ "portable", "portable" },
        .{ "vec", "vec" },         .{ "vec x4", "vec x4" }, .{ "gemv N", "gemv N" },
        .{ "CSC A x", "CSC A x" }, .{ "gemv T", "gemv T" },
        .{ "CSC Aᵀ y", "CSC Aᵀ y" },
        .{ "4 sums", "4 sums" },   .{ "gather", "gather" }, .{ "inline", "inline" },
        .{ "library", "library" }, .{ "plain", "plain" },   .{ "comp.", "comp." },
        .{ "fast", "fast" },
    }) |c| l = l ++ [_]merge.Label{.{ .label = c[0], .short = c[1] }};
    break :blk l[0..l.len].*;
};

/// Marks the start of a section's output for `merge.zig`, which reruns
/// noisy sections by name.
pub fn sectionMark(comptime name: []const u8) void {
    section = name;
    std.debug.print("\n<!-- section " ++ name ++ " -->\n", .{});
}

/// A row after its lead: `printHeader`'s columns; `-` where one is missing,
/// `!` on a difference past its bound.
pub fn printRow(cols: []const Col, rates: []const ?f64, diffs: []const ?Diff, mine: Rstd) void {
    for (cols, rates) |c, r| {
        if (!c.shown()) continue;
        if (r) |v| std.debug.print("{d:.1} | ", .{v}) else std.debug.print("- | ", .{});
    }
    for (cols[1..], rates[1..]) |c, r| {
        if (!c.shown()) continue;
        if (r != null and rates[0] != null) std.debug.print("{d:.2} | ", .{r.? / rates[0].?}) else std.debug.print("- | ", .{});
    }
    if (fair()) {
        const o = obAlone();
        if (o) |v| std.debug.print("{d:.1} | ", .{v}) else std.debug.print("- | ", .{});
        if (o != null and rates[rstd_col] != null) std.debug.print("{d:.2} | ", .{rates[rstd_col].? / o.?}) else std.debug.print("- | ", .{});
    }
    if (bb.alone) return std.debug.print("\n", .{});
    for (cols, diffs) |c, d| {
        if (c.lib == null or !c.shown()) continue;
        if (d) |v| std.debug.print("{d:.1}{s} | ", .{ v.eps, if (v.eps > v.bound) " !" else "" }) else std.debug.print("- | ", .{});
    }
    if (mine.err) |e| std.debug.print("{d:.1} | ", .{e}) else std.debug.print("- | ", .{});
    if (mine.bits) |h| std.debug.print("{x:0>8} |\n", .{@as(u32, @truncate(h))}) else std.debug.print("- |\n", .{});
}

/// rstd's error from the f64 reference, in eps of the type relative to the
/// size of the terms, and a hash of its output.
pub const Rstd = struct { err: ?f64 = null, bits: ?u64 = null };

/// A library table row's contestants, `Ctx`s entered in checking order
/// (`rstdFirst`) for `n` columns, and what `printRow` shows beside their rates.
pub fn Entrants(comptime Ctx: type, comptime n: usize) type {
    return struct {
        ctxs: [n]Ctx = undefined,
        /// Each column's place among `ctxs`; null for one not entered.
        at: [n]?usize = @splat(null),
        diffs: [n]?Diff = @splat(null),
        mine: Rstd = .{},
        count: usize = 0,

        /// Enters column `col`'s contestant.
        pub fn add(self: *@This(), col: usize, ctx: Ctx) *Ctx {
            self.ctxs[self.count] = ctx;
            self.at[col] = self.count;
            self.count += 1;
            return &self.ctxs[self.count - 1];
        }

        /// Races them; each column's rate, `work` (bytes or flops) per call
        /// in billions per second, null where it was not entered.
        pub fn race(self: *@This(), io: Io, timing: bench.Timing, work: f64) ![n]?f64 {
            var calls: [n]bench.Call = undefined;
            for (calls[0..self.count], self.ctxs[0..self.count]) |*c, *x| c.* = .of(x);
            var secs: [n]f64 = undefined;
            try bench.race(io, timing, calls[0..self.count], &secs);
            if (lb.lib_err) |err| return err;
            var rates: [n]?f64 = undefined;
            for (self.at, &rates) |j, *rt| rt.* = if (j) |k| work / secs[k] / 1e9 else null;
            return rates;
        }

        /// The row after its lead.
        pub fn print(self: *const @This(), cols: []const Col, rates: []const ?f64) void {
            printRow(cols, rates, &self.diffs, self.mine);
        }
    };
}

/// What rstd's own table holds in slot `name`: converted `asm`, or
/// `portable`.
pub fn rstdRuns(comptime name: []const u8) []const u8 {
    if (bb.alone) return "-";
    return if (@field(blas.vtable().*, name).? == @field(bb.tables.portable, name).?) "portable" else "asm";
}

pub const diff_note = "`rstd-<lib>` is rstd's largest difference from that library's result, in eps of the type relative to the size of the terms; `!` marks one past what accumulation order explains: 0 for `copy`, `swap` and `i?amax`, 8 for the other element-wise routines, 2 (n + 8) for a sum of n terms. `rstd err` is rstd's own error from an f64 reference in the same units (gemm: at 256 entries), which tells whose result a difference comes from. `rstd bits` hashes rstd's output, the same on machines that run the same kernels";

/// The section being printed (`sectionMark`), the lead of its table's header
/// (`printHeader`) and the lead of its row (`printLead`): a row's key in the alone
/// run's rates.
var section: []const u8 = "";
var header_lead: []const u8 = "";
var row_lead_buf: [256]u8 = undefined;
var row_lead: []const u8 = "";

/// Prints a row's lead, the cells before `printRow`'s.
pub fn printLead(comptime fmt: []const u8, args: anytype) void {
    row_lead = std.mem.print(&row_lead_buf, fmt, args) catch "";
    std.debug.print("{s}", .{row_lead});
}

/// The race shows rstd against the alone run's OpenBLAS (issue #690).
fn fair() bool {
    return !bb.alone and alone.ob_rates != null;
}

/// This row's OpenBLAS rate in the alone run.
fn obAlone() ?f64 {
    var hc: [8][]const u8 = undefined;
    var rc: [8][]const u8 = undefined;
    const h = splitCells(&hc, header_lead);
    const r = splitCells(&rc, row_lead);
    var kb: [512]u8 = undefined;
    const k = alone.key(&kb, section, h, r) catch return null;
    return alone.ob_rates.?.get(k);
}

/// `line`'s cells between `|`s, trimmed, empty ones out.
fn splitCells(out: [][]const u8, line: []const u8) []const []const u8 {
    var n: usize = 0;
    var it = std.mem.splitScalar(u8, line, '|');
    while (it.next()) |c| {
        const t = std.mem.trim(u8, c, " ");
        if (t.len == 0 or n == out.len) continue;
        out[n] = t;
        n += 1;
    }
    return out[0..n];
}
