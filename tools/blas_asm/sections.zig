//! Cuts a long generated text into sections of at most `cap` lines, for the
//! generators' one-file-per-section output (`gen.zig`'s `text<nn>.zig`).
//!
//! An assembly text is cut only at a seam: a label or a `.macro` outside any
//! macro, with the alignment directives before it, so a section is a run of
//! whole labeled blocks or macros. The cuts are the fewest that keep every
//! section under the cap, then the most even; a run with no seam stays whole
//! and over the cap.

const std = @import("std");

/// The longest file the generators write in one piece (the layout rule).
pub const max_lines = 500;

/// One line of an assembly text (not `gen.zig`'s `__bk` wrapper); `top` when no
/// generated `.if \bpN` holds it.
pub const Line = struct { text: []const u8, top: bool = true };

fn startsWithAny(t: []const u8, list: []const []const u8) bool {
    for (list) |p| if (std.mem.startsWith(u8, t, p)) return true;
    return false;
}

/// Whether the line starts by defining a numeric label, as every label is
/// once `gen.zig` has renumbered them.
fn isLabel(t: []const u8) bool {
    var i: usize = 0;
    while (i < t.len and std.ascii.isDigit(t[i])) i += 1;
    return i > 0 and i < t.len and t[i] == ':';
}

/// The line indexes a section may start at, in order.
pub fn seams(a: std.mem.Allocator, lines: []const Line) ![]const usize {
    var out: std.ArrayList(usize) = .empty;
    var in_macro = false;
    for (lines, 0..) |l, i| {
        const t = std.mem.trim(u8, l.text, " ");
        const macro = std.mem.startsWith(u8, t, ".macro ");
        if (l.top and !in_macro and (macro or isLabel(t))) {
            // the block's alignment goes with it
            var at = i;
            while (at > 0 and lines[at - 1].top and startsWithAny(std.mem.trim(u8, lines[at - 1].text, " "), &.{ ".align ", ".balign ", ".p2align " })) at -= 1;
            if (at > 0 and (out.items.len == 0 or out.items[out.items.len - 1] < at)) try out.append(a, at);
        }
        if (macro) in_macro = true;
        if (std.mem.eql(u8, t, ".endm")) in_macro = false;
    }
    return out.items;
}

/// Where each section of `len` lines starts (the first at 0), cutting only at
/// `at` (ascending, within the text).
pub fn cut(a: std.mem.Allocator, len: usize, at: []const usize, cap: usize) ![]const usize {
    // the candidates' ends: every seam, then the text's end
    const ends = try a.alloc(usize, at.len + 1);
    @memcpy(ends[0..at.len], at);
    ends[at.len] = len;
    // `cost[j]`: the best cost of sections ending at `ends[j]`, from `from[j]`
    // (an index into `ends`, or none for the text's start)
    const cost = try a.alloc(u64, ends.len);
    const from = try a.alloc(?usize, ends.len);
    for (ends, 0..) |end, j| {
        cost[j] = sectionCost(end, cap);
        from[j] = null;
        for (ends[0..j], 0..) |start, i| {
            const c = cost[i] + sectionCost(end - start, cap);
            if (c < cost[j]) {
                cost[j] = c;
                from[j] = i;
            }
        }
    }
    var rev: std.ArrayList(usize) = .empty;
    var j = from[ends.len - 1];
    while (j) |i| : (j = from[i]) try rev.append(a, ends[i]);
    try rev.append(a, 0);
    std.mem.reverse(usize, rev.items);
    return rev.items;
}

/// A section over the cap costs more than any number of sections, a section
/// more than any evening out of their lengths.
fn sectionCost(n: usize, cap: usize) u64 {
    const over: u64 = if (n > cap) 1 << 48 else 0;
    return over + @as(u64, cap) * cap + @as(u64, n) * n;
}

test cut {
    var arena: std.heap.ArenaAllocator = .init(std.testing.allocator);
    defer arena.deinit();
    const a = arena.allocator();
    try std.testing.expectEqualSlices(usize, &.{0}, try cut(a, 90, &.{ 30, 60 }, 100));
    try std.testing.expectEqualSlices(usize, &.{ 0, 60 }, try cut(a, 120, &.{ 30, 60, 110 }, 100));
    try std.testing.expectEqualSlices(usize, &.{ 0, 90, 150 }, try cut(a, 240, &.{ 10, 90, 150 }, 100));
    // no seam in reach: the long run stays whole
    try std.testing.expectEqualSlices(usize, &.{ 0, 150 }, try cut(a, 200, &.{150}, 100));
}

test seams {
    var arena: std.heap.ArenaAllocator = .init(std.testing.allocator);
    defer arena.deinit();
    const lines = [_]Line{
        .{ .text = " mov 1, 2" },
        .{ .text = " .balign 16" },
        .{ .text = "20001:" },
        .{ .text = " .if \\bp0", .top = false },
        .{ .text = "20002:", .top = false },
        .{ .text = " .endif", .top = false },
        .{ .text = " .macro LOAD x" },
        .{ .text = "1: ld \\x" },
        .{ .text = " .endm" },
        .{ .text = "3: b 1b" },
    };
    try std.testing.expectEqualSlices(usize, &.{ 1, 6, 9 }, try seams(arena.allocator(), &lines));
}
