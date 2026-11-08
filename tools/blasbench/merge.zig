//! `blas_bench --merge <out.md> <log>...`: merges runs of one configuration by the timing rule (`readme.md`): a row agreeing within 3%
//! shows the pair's better rates, else is marked noisy; deterministic cells must match. Prints `NOISY` and `MISMATCH` lines for `flyoff.sh`.

const std = @import("std");
const Io = std.Io;
const Allocator = std.mem.Allocator;

/// A contestant's column header and its name in a ratio's (`x/ob`).
pub const Label = struct { label: []const u8, short: []const u8 };

/// Two runs agree when every rate is within this of the other.
const tolerance = 0.03;

const Table = struct {
    header: []const []const u8,
    sep: []const u8,
    /// the last `<!-- section <name> -->` before it
    section: []const u8,
    keys: std.ArrayList([]const u8) = .empty,
    rows: std.StringHashMapUnmanaged([]const []const u8) = .empty,
};

const Item = union(enum) { text: []const u8, table: *Table };

pub fn run(io: Io, gpa: Allocator, args: []const [:0]const u8, labels: []const Label) !void {
    if (args.len < 2) {
        std.debug.print("usage: blas_bench --merge <out.md> <log>...\n", .{});
        return error.Usage;
    }
    var arena_state: std.heap.ArenaAllocator = .init(gpa);
    defer arena_state.deinit();
    const arena = arena_state.allocator();
    var runs: std.ArrayList([]const Item) = .empty;
    for (args[1..]) |path| {
        const text = try Io.Dir.cwd().readFileAlloc(io, path, arena, .limited(64 << 20));
        try runs.append(arena, try parse(arena, text, labels));
    }

    var md: std.ArrayList(u8) = .empty;
    var report: std.ArrayList(u8) = .empty;
    for (runs.items[0]) |item| switch (item) {
        .text => |t| try md.print(arena, "{s}\n", .{t}),
        .table => |t| try mergeTable(arena, t, runs.items, labels, &md, &report),
    };
    try Io.Dir.cwd().writeFile(io, .{ .sub_path = args[0], .data = md.items });
    const out = if (report.items.len == 0) "all rows agree\n" else report.items;
    try Io.File.stdout().writeStreamingAll(io, out);
}

fn isRate(h: []const u8, labels: []const Label) bool {
    for (labels) |l| if (std.mem.eql(u8, h, l.label)) return true;
    return false;
}

fn isDeterministic(h: []const u8) bool {
    return std.mem.startsWith(u8, h, "rstd-") or std.mem.eql(u8, h, "rstd err") or std.mem.eql(u8, h, "rstd bits");
}

fn cells(arena: Allocator, line: []const u8) ![]const []const u8 {
    const inner = std.mem.trim(u8, std.mem.trim(u8, line, " \t\r"), "|");
    var list: std.ArrayList([]const u8) = .empty;
    var it = std.mem.splitScalar(u8, inner, '|');
    while (it.next()) |c| try list.append(arena, std.mem.trim(u8, c, " \t"));
    return list.items;
}

fn parse(arena: Allocator, text: []const u8, labels: []const Label) ![]const Item {
    var lines: std.ArrayList([]const u8) = .empty;
    var it = std.mem.splitScalar(u8, text, '\n');
    while (it.next()) |l| try lines.append(arena, std.mem.trimEnd(u8, l, "\r"));
    const ls = lines.items;
    var items: std.ArrayList(Item) = .empty;
    var section: []const u8 = "?";
    const mark = "<!-- section ";
    var i: usize = 0;
    while (i < ls.len) {
        const ln = ls[i];
        if (std.mem.startsWith(u8, ln, mark)) section = std.mem.trimEnd(u8, ln[mark.len..], " ->");
        if (std.mem.startsWith(u8, ln, "| ") and i + 1 < ls.len and std.mem.startsWith(u8, ls[i + 1], "| -")) {
            const t = try arena.create(Table);
            t.* = .{ .header = try cells(arena, ln), .sep = ls[i + 1], .section = section };
            const first = firstRate(t.header, labels) orelse t.header.len;
            i += 2;
            while (i < ls.len and std.mem.startsWith(u8, ls[i], "| ")) : (i += 1) {
                const c = try cells(arena, ls[i]);
                const key = try std.mem.join(arena, "|", c[0..@min(first, c.len)]);
                if (!t.rows.contains(key)) try t.keys.append(arena, key);
                try t.rows.put(arena, key, c);
            }
            try items.append(arena, .{ .table = t });
            continue;
        }
        try items.append(arena, .{ .text = ln });
        i += 1;
    }
    return items.items;
}

fn firstRate(header: []const []const u8, labels: []const Label) ?usize {
    for (header, 0..) |h, j| if (isRate(h, labels)) return j;
    return null;
}

fn sameHeader(a: []const []const u8, b: []const []const u8) bool {
    if (a.len != b.len) return false;
    for (a, b) |x, y| if (!std.mem.eql(u8, x, y)) return false;
    return true;
}

/// A rate cell's value; null for `-` and the like.
fn num(s: []const u8) ?f64 {
    const v = std.mem.trim(u8, s, " !");
    return std.fmt.parseFloat(f64, v) catch null;
}

fn mergeTable(arena: Allocator, base: *const Table, runs: []const []const Item, labels: []const Label, md: *std.ArrayList(u8), report: *std.ArrayList(u8)) !void {
    const h = base.header;
    try md.print(arena, "| {s} |\n{s}\n", .{ try std.mem.join(arena, " | ", h), base.sep });
    for (base.keys.items) |key| {
        // each run's row of this key in a table of this header
        var cand: std.ArrayList([]const []const u8) = .empty;
        for (runs) |r| for (r) |item| switch (item) {
            .table => |t| if (sameHeader(t.header, h)) if (t.rows.get(key)) |row| if (row.len == h.len) try cand.append(arena, row),
            .text => {},
        };
        for (h, 0..) |name, j| {
            if (!isDeterministic(name)) continue;
            for (cand.items[1..]) |c| if (!std.mem.eql(u8, c[j], cand.items[0][j])) {
                try report.print(arena, "MISMATCH {s} {s}: {s} vs {s}\n", .{ name, key, cand.items[0][j], c[j] });
                break;
            };
        }
        // the closest pair
        var best: ?[2]usize = null;
        var best_dev: f64 = 0;
        for (0..cand.items.len) |a| for (a + 1..cand.items.len) |b| {
            var dev: f64 = 0;
            for (h, 0..) |name, j| {
                if (!isRate(name, labels)) continue;
                const x = num(cand.items[a][j]) orelse continue;
                const y = num(cand.items[b][j]) orelse continue;
                if (x == 0 or y == 0) continue;
                dev = @max(dev, @abs(x - y) / @max(x, y));
            }
            if (best == null or dev < best_dev) {
                best = .{ a, b };
                best_dev = dev;
            }
        };
        const noisy = best == null or best_dev > tolerance;
        const row = try arena.dupe([]const u8, cand.items[0]);
        for (h, 0..) |name, j| {
            if (!isRate(name, labels)) continue;
            var top: ?f64 = null;
            for (cand.items, 0..) |c, k| {
                if (!noisy and k != best.?[0] and k != best.?[1]) continue;
                const v = num(c[j]) orelse continue;
                top = if (top) |t| @max(t, v) else v;
            }
            if (top) |t| row[j] = try std.fmt.allocPrint(arena, "{d:.1}", .{t});
        }
        // ratios `x/ob` over the merged rates
        for (h, 0..) |name, j| {
            const slash = std.mem.findScalar(u8, name, '/') orelse continue;
            const nj = column(h, labels, name[0..slash]) orelse continue;
            const dj = column(h, labels, name[slash + 1 ..]) orelse continue;
            const x = num(row[nj]) orelse continue;
            const y = num(row[dj]) orelse continue;
            if (y != 0) row[j] = try std.fmt.allocPrint(arena, "{d:.2}", .{x / y});
        }
        if (noisy) {
            row[0] = try std.fmt.allocPrint(arena, "{s} (noisy)", .{row[0]});
            if (best == null)
                try report.print(arena, "NOISY {s} {s}: dev n/a over {d} runs\n", .{ base.section, key, cand.items.len })
            else
                try report.print(arena, "NOISY {s} {s}: dev {d:.4} over {d} runs\n", .{ base.section, key, best_dev, cand.items.len });
        }
        try md.print(arena, "| {s} |\n", .{try std.mem.join(arena, " | ", row)});
    }
}

/// The column of the contestant whose short name is `short`.
fn column(header: []const []const u8, labels: []const Label, short: []const u8) ?usize {
    for (labels) |l| if (std.mem.eql(u8, l.short, short)) {
        for (header, 0..) |name, j| if (std.mem.eql(u8, name, l.label)) return j;
    };
    return null;
}

test "a pair within 3% merges to its better rates" {
    const labels = [_]Label{ .{ .label = "openblas", .short = "ob" }, .{ .label = "rstd", .short = "rstd" } };
    const a = "| p | n | openblas | rstd | rstd/ob | rstd err |\n| - | -: | -: | -: | -: | -: |\n| d | 64 | 10.0 | 20.0 | 2.00 | 1.0 |\n";
    const b = "| p | n | openblas | rstd | rstd/ob | rstd err |\n| - | -: | -: | -: | -: | -: |\n| d | 64 | 10.2 | 19.8 | 1.94 | 1.0 |\n";
    var arena_state: std.heap.ArenaAllocator = .init(std.testing.allocator);
    defer arena_state.deinit();
    const arena = arena_state.allocator();
    const runs = [_][]const Item{ try parse(arena, a, &labels), try parse(arena, b, &labels) };
    var md: std.ArrayList(u8) = .empty;
    var report: std.ArrayList(u8) = .empty;
    try mergeTable(arena, runs[0][0].table, &runs, &labels, &md, &report);
    try std.testing.expectEqualStrings("", report.items);
    try std.testing.expect(std.mem.find(u8, md.items, "| d | 64 | 10.2 | 20.0 | 1.96 | 1.0 |") != null);
}
