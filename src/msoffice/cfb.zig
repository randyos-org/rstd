//! Minimal [MS-CFB] compound file writer (version 3, 512-byte sectors).
//! Streams smaller than the 4096-byte cutoff go in the mini stream, as the
//! format requires. Timestamps and CLSIDs are zero so output is reproducible.
//! Entry names must be ASCII (at most 31 characters).
const std = @import("std");
const zigstd = @import("../zigstd.zig");
const bytes = @import("../bytes.zig");

const Allocator = zigstd.Allocator;
const ByteArrayList = bytes.ByteArrayList;

/// A storage (directory) or stream (file) in the compound file.
pub const Node = struct {
    name: []const u8,
    kind: union(enum) {
        stream: []const u8,
        storage: []const Node,
    },
};

const signature = [8]u8{ 0xD0, 0xCF, 0x11, 0xE0, 0xA1, 0xB1, 0x1A, 0xE1 };
const sector_size = 512;
const mini_sector_size = 64;
const mini_stream_cutoff = 4096;
const dir_entry_size = 128;
const ids_per_sector = sector_size / 4;
const header_difat_len = 109;

const max_reg_sect: u32 = 0xFFFFFFFA;
const fat_sect: u32 = 0xFFFFFFFD;
const end_of_chain: u32 = 0xFFFFFFFE;
const free_sect: u32 = 0xFFFFFFFF;
const no_stream: u32 = 0xFFFFFFFF;

const ObjectType = enum(u8) { unallocated = 0, storage = 1, stream = 2, root = 5 };
const black: u8 = 1;

const DirEntry = struct {
    name: []const u8,
    type: ObjectType,
    left: u32 = no_stream,
    right: u32 = no_stream,
    child: u32 = no_stream,
    data: []const u8 = &.{},
    start: u32 = end_of_chain,
    size: u32 = 0,
};

/// [MS-CFB] 2.6.4 sibling order: shorter names first, then by uppercased
/// code unit.
fn nameLessThan(_: void, a: DirEntry, b: DirEntry) bool {
    if (a.name.len != b.name.len) return a.name.len < b.name.len;
    for (a.name, b.name) |ca, cb| {
        const ua = std.ascii.toUpper(ca);
        const ub = std.ascii.toUpper(cb);
        if (ua != ub) return ua < ub;
    }
    return false;
}

const Builder = struct {
    gpa: Allocator,
    entries: std.ArrayList(DirEntry) = .empty,

    /// Appends `nodes` as siblings, recursing into storages, and returns the
    /// index of the root of their (balanced) sibling tree.
    fn addSiblings(self: *Builder, nodes: []const Node) !u32 {
        if (nodes.len == 0) return no_stream;
        const first: u32 = @intCast(self.entries.items.len);
        for (nodes) |n| {
            if (n.name.len == 0 or n.name.len > 31) return error.CfbNameLength;
            for (n.name) |c| if (c >= 0x80 or c == '/' or c == '\\' or c == ':' or c == '!') return error.CfbNameInvalid;
            try self.entries.append(self.gpa, switch (n.kind) {
                .stream => |data| .{ .name = n.name, .type = .stream, .data = data },
                .storage => .{ .name = n.name, .type = .storage },
            });
        }
        const siblings = self.entries.items[first..];
        // Sorting reorders the entries, so pair each with its node by name.
        std.mem.sort(DirEntry, siblings, {}, nameLessThan);
        for (siblings, 1..) |e, i| {
            if (i < siblings.len and !nameLessThan({}, e, siblings[i])) return error.CfbDuplicateName;
        }
        for (first..first + nodes.len) |idx| {
            if (self.entries.items[idx].type != .storage) continue;
            const name = self.entries.items[idx].name;
            const children = for (nodes) |n| {
                if (zigstd.equal(u8, n.name, name)) break n.kind.storage;
            } else unreachable;
            const child = try self.addSiblings(children);
            self.entries.items[idx].child = child;
        }
        return self.link(first, first + @as(u32, @intCast(nodes.len)));
    }

    /// Links sorted entries [lo, hi) into a balanced binary tree; returns its
    /// root. All nodes are black, which readers accept.
    fn link(self: *Builder, lo: u32, hi: u32) u32 {
        if (lo >= hi) return no_stream;
        const mid = lo + (hi - lo) / 2;
        self.entries.items[mid].left = self.link(lo, mid);
        self.entries.items[mid].right = self.link(mid + 1, hi);
        return mid;
    }
};

fn sectorsFor(len: usize, unit: usize) u32 {
    return @intCast((len + unit - 1) / unit);
}

/// Builds a complete compound file whose root storage holds `nodes`.
/// The result is owned by the caller.
pub fn write(gpa: Allocator, nodes: []const Node) ![]u8 {
    var arena_state: std.heap.ArenaAllocator = .init(gpa);
    defer arena_state.deinit();
    const arena = arena_state.allocator();

    var b: Builder = .{ .gpa = arena };
    try b.entries.append(arena, .{ .name = "Root Entry", .type = .root });
    // Not `entries.items[0].child = try ...`: adding siblings can reallocate.
    const root_child = try b.addSiblings(nodes);
    b.entries.items[0].child = root_child;
    const entries = b.entries.items;

    // Mini stream: every non-empty stream below the cutoff, 64-byte aligned.
    var mini_stream: ByteArrayList = .empty;
    var mini_fat: std.ArrayList(u32) = .empty;
    for (entries) |*e| {
        if (e.type != .stream) continue;
        e.size = std.math.cast(u32, e.data.len) orelse return error.CfbStreamTooLarge;
        if (e.data.len == 0 or e.data.len >= mini_stream_cutoff) continue;
        const first: u32 = @intCast(mini_fat.items.len);
        const n = sectorsFor(e.data.len, mini_sector_size);
        for (0..n) |i| try mini_fat.append(arena, if (i + 1 == n) end_of_chain else first + @as(u32, @intCast(i)) + 1);
        e.start = first;
        try mini_stream.appendSlice(arena, e.data);
        try mini_stream.appendNTimes(arena, 0, n * mini_sector_size - e.data.len);
    }

    // Regular sector layout: large streams, mini stream, mini FAT,
    // directory, then the FAT itself.
    var fat: std.ArrayList(u32) = .empty;
    const Chain = struct {
        fn add(list: *std.ArrayList(u32), a: Allocator, n: u32) !u32 {
            if (n == 0) return end_of_chain;
            const first: u32 = @intCast(list.items.len);
            for (0..n) |i| try list.append(a, if (i + 1 == n) end_of_chain else first + @as(u32, @intCast(i)) + 1);
            return first;
        }
    };
    for (entries) |*e| {
        if (e.type == .stream and e.data.len >= mini_stream_cutoff) {
            e.start = try Chain.add(&fat, arena, sectorsFor(e.data.len, sector_size));
        }
    }
    entries[0].start = try Chain.add(&fat, arena, sectorsFor(mini_stream.items.len, sector_size));
    entries[0].size = @intCast(mini_stream.items.len);
    const mini_fat_sectors = sectorsFor(mini_fat.items.len * 4, sector_size);
    const mini_fat_start = try Chain.add(&fat, arena, mini_fat_sectors);
    const dir_sectors = sectorsFor(entries.len * dir_entry_size, sector_size);
    const dir_start = try Chain.add(&fat, arena, dir_sectors);

    var fat_sectors: u32 = 0;
    while (true) {
        const needed = sectorsFor(fat.items.len + fat_sectors, ids_per_sector);
        if (needed == fat_sectors) break;
        fat_sectors = needed;
    }
    if (fat_sectors > header_difat_len) return error.CfbTooLarge;
    const fat_start: u32 = @intCast(fat.items.len);
    try fat.appendNTimes(arena, fat_sect, fat_sectors);
    if (fat.items.len > max_reg_sect) return error.CfbTooLarge;
    try fat.appendNTimes(arena, free_sect, fat_sectors * ids_per_sector - fat.items.len);

    var out: ByteArrayList = .empty;
    defer out.deinit(gpa);
    const total_sectors = fat_start + fat_sectors;
    try out.ensureTotalCapacity(gpa, sector_size * (1 + total_sectors));

    // Header.
    try out.appendSlice(gpa, &signature);
    try out.appendNTimes(gpa, 0, 16); // CLSID
    try bytes.appendIntLE(gpa, &out, u16, 0x003E); // minor version
    try bytes.appendIntLE(gpa, &out, u16, 0x0003); // major version 3
    try bytes.appendIntLE(gpa, &out, u16, 0xFFFE); // byte order
    try bytes.appendIntLE(gpa, &out, u16, 9); // sector shift
    try bytes.appendIntLE(gpa, &out, u16, 6); // mini sector shift
    try out.appendNTimes(gpa, 0, 6);
    try bytes.appendIntLE(gpa, &out, u32, 0); // directory sectors (0 for v3)
    try bytes.appendIntLE(gpa, &out, u32, fat_sectors);
    try bytes.appendIntLE(gpa, &out, u32, dir_start);
    try bytes.appendIntLE(gpa, &out, u32, 0); // transaction signature
    try bytes.appendIntLE(gpa, &out, u32, mini_stream_cutoff);
    try bytes.appendIntLE(gpa, &out, u32, mini_fat_start);
    try bytes.appendIntLE(gpa, &out, u32, mini_fat_sectors);
    try bytes.appendIntLE(gpa, &out, u32, end_of_chain); // first DIFAT sector
    try bytes.appendIntLE(gpa, &out, u32, 0); // DIFAT sectors
    for (0..header_difat_len) |i| {
        try bytes.appendIntLE(gpa, &out, u32, if (i < fat_sectors) fat_start + @as(u32, @intCast(i)) else free_sect);
    }
    zigstd.assert(out.items.len == sector_size);

    // Sectors, in the order they were allocated above.
    for (entries) |e| {
        if (e.type == .stream and e.data.len >= mini_stream_cutoff) try appendPadded(gpa, &out, e.data);
    }
    try appendPadded(gpa, &out, mini_stream.items);
    {
        const start = out.items.len;
        for (mini_fat.items) |id| try bytes.appendIntLE(gpa, &out, u32, id);
        try padIds(gpa, &out, start);
    }
    {
        const start = out.items.len;
        for (entries) |e| try appendDirEntry(gpa, &out, e);
        const unused = dir_sectors * (sector_size / dir_entry_size) - entries.len;
        for (0..unused) |_| try appendDirEntry(gpa, &out, .{ .name = "", .type = .unallocated });
        zigstd.assert(out.items.len - start == dir_sectors * sector_size);
    }
    for (fat.items) |id| try bytes.appendIntLE(gpa, &out, u32, id);
    zigstd.assert(out.items.len == sector_size * (1 + total_sectors));

    return out.toOwnedSlice(gpa);
}

fn appendPadded(gpa: Allocator, out: *ByteArrayList, data: []const u8) !void {
    try out.appendSlice(gpa, data);
    try out.appendNTimes(gpa, 0, sectorsFor(data.len, sector_size) * sector_size - data.len);
}

/// Pads a run of sector IDs that began at `start` with FREESECT to a whole
/// sector.
fn padIds(gpa: Allocator, out: *ByteArrayList, start: usize) !void {
    while ((out.items.len - start) % sector_size != 0) try bytes.appendIntLE(gpa, out, u32, free_sect);
}

fn appendDirEntry(gpa: Allocator, out: *ByteArrayList, e: DirEntry) !void {
    const start = out.items.len;
    for (e.name) |c| try bytes.appendIntLE(gpa, out, u16, c);
    try out.appendNTimes(gpa, 0, 64 - e.name.len * 2);
    const name_len: u16 = if (e.name.len == 0) 0 else @intCast((e.name.len + 1) * 2);
    try bytes.appendIntLE(gpa, out, u16, name_len);
    try out.append(gpa, @backingInt(e.type));
    try out.append(gpa, if (e.type == .unallocated) 0 else black);
    try bytes.appendIntLE(gpa, out, u32, e.left);
    try bytes.appendIntLE(gpa, out, u32, e.right);
    try bytes.appendIntLE(gpa, out, u32, e.child);
    try out.appendNTimes(gpa, 0, 16); // CLSID
    try bytes.appendIntLE(gpa, out, u32, 0); // state bits
    try bytes.appendIntLE(gpa, out, u64, 0); // creation time
    try bytes.appendIntLE(gpa, out, u64, 0); // modified time
    const is_stream = e.type == .stream or e.type == .root;
    try bytes.appendIntLE(gpa, out, u32, if (is_stream) e.start else 0);
    try bytes.appendIntLE(gpa, out, u64, if (is_stream) e.size else 0);
    zigstd.assert(out.items.len - start == dir_entry_size);
}

// ---------------------------------------------------------------------------
// Tests: a small reader that follows the same rules, for round-tripping.
// ---------------------------------------------------------------------------

const TestReader = struct {
    file: []const u8,
    fat: []const u32,
    mini_fat: []const u32,
    mini_stream: []const u8,
    dir: []const u8,
    arena: Allocator,

    fn sector(self: TestReader, id: u32) []const u8 {
        const off = (@as(usize, id) + 1) * sector_size;
        return self.file[off..][0..sector_size];
    }

    fn readChain(self: TestReader, start: u32) ![]u8 {
        var buf: ByteArrayList = .empty;
        var s = start;
        var guard: usize = 0;
        while (s != end_of_chain) : (guard += 1) {
            if (guard > self.fat.len) return error.ChainLoop;
            try buf.appendSlice(self.arena, self.sector(s));
            s = self.fat[s];
        }
        return buf.items;
    }

    fn init(arena: Allocator, file: []const u8) !TestReader {
        if (!zigstd.equal(u8, file[0..8], &signature)) return error.BadSignature;
        var r: TestReader = .{ .file = file, .fat = &.{}, .mini_fat = &.{}, .mini_stream = &.{}, .dir = &.{}, .arena = arena };
        const n_fat = bytes.readIntLE(u32, file[44..48]);
        var fat: std.ArrayList(u32) = .empty;
        for (0..n_fat) |i| {
            const sid = bytes.readIntLE(u32, file[76 + 4 * i ..][0..4]);
            const s = r.sector(sid);
            for (0..ids_per_sector) |j| try fat.append(arena, bytes.readIntLE(u32, s[4 * j ..][0..4]));
        }
        r.fat = fat.items;
        r.dir = try r.readChain(bytes.readIntLE(u32, file[48..52]));
        const mf = bytes.readIntLE(u32, file[60..64]);
        if (mf != end_of_chain) {
            const raw = try r.readChain(mf);
            var ids: std.ArrayList(u32) = .empty;
            for (0..raw.len / 4) |j| try ids.append(arena, bytes.readIntLE(u32, raw[4 * j ..][0..4]));
            r.mini_fat = ids.items;
        }
        const root_start = bytes.readIntLE(u32, r.dir[116..120]);
        if (root_start != end_of_chain) r.mini_stream = try r.readChain(root_start);
        return r;
    }

    fn entry(self: TestReader, idx: u32) []const u8 {
        return self.dir[@as(usize, idx) * dir_entry_size ..][0..dir_entry_size];
    }

    fn entryName(self: TestReader, buf: []u8, idx: u32) []const u8 {
        const e = self.entry(idx);
        const len = bytes.readIntLE(u16, e[64..66]) / 2 - 1;
        for (0..len) |i| buf[i] = @intCast(bytes.readIntLE(u16, e[2 * i ..][0..2]));
        return buf[0..len];
    }

    /// Finds `name` among the siblings rooted at `idx`.
    fn find(self: TestReader, idx: u32, name: []const u8) ?u32 {
        if (idx == no_stream) return null;
        var buf: [31]u8 = undefined;
        const here = self.entryName(&buf, idx);
        if (zigstd.equal(u8, here, name)) return idx;
        const e = self.entry(idx);
        const probe: DirEntry = .{ .name = name, .type = .stream };
        const cur: DirEntry = .{ .name = here, .type = .stream };
        const next = if (nameLessThan({}, probe, cur)) bytes.readIntLE(u32, e[68..72]) else bytes.readIntLE(u32, e[72..76]);
        return self.find(next, name);
    }

    fn stream(self: TestReader, path: []const []const u8) ![]const u8 {
        var idx: u32 = 0;
        for (path) |name| {
            const child = bytes.readIntLE(u32, self.entry(idx)[76..80]);
            idx = self.find(child, name) orelse return error.NotFound;
        }
        const e = self.entry(idx);
        const start = bytes.readIntLE(u32, e[116..120]);
        const size = bytes.readIntLE(u32, e[120..124]);
        if (size >= mini_stream_cutoff) return (try self.readChain(start))[0..size];
        var buf: ByteArrayList = .empty;
        var s = start;
        while (s != end_of_chain) : (s = self.mini_fat[s]) {
            try buf.appendSlice(self.arena, self.mini_stream[@as(usize, s) * mini_sector_size ..][0..mini_sector_size]);
        }
        return buf.items[0..size];
    }
};

test "round trip small, large, and nested streams" {
    var arena_state: std.heap.ArenaAllocator = .init(zigstd.testAllocator);
    defer arena_state.deinit();
    const arena = arena_state.allocator();

    const xs: [100]u8 = @splat('x');
    const big = try arena.alloc(u8, 70_000);
    for (big, 0..) |*c, i| c.* = @truncate(i *% 31);
    const nodes = [_]Node{
        .{ .name = "PROJECT", .kind = .{ .stream = "small stream" } },
        .{ .name = "Big", .kind = .{ .stream = big } },
        .{ .name = "Empty", .kind = .{ .stream = "" } },
        .{ .name = "VBA", .kind = .{ .storage = &.{
            .{ .name = "dir", .kind = .{ .stream = &xs } },
            .{ .name = "Module1", .kind = .{ .stream = big[0..5000] } },
            .{ .name = "_VBA_PROJECT", .kind = .{ .stream = "abc" } },
        } } },
    };
    const file = try write(arena, &nodes);
    try zigstd.testExpectEqual(@as(usize, 0), file.len % sector_size);

    const r = try TestReader.init(arena, file);
    try zigstd.testExpectEqualStrings("small stream", try r.stream(&.{"PROJECT"}));
    try zigstd.testExpectEqualSlices(u8, big, try r.stream(&.{"Big"}));
    try zigstd.testExpectEqualStrings("", try r.stream(&.{"Empty"}));
    try zigstd.testExpectEqualStrings(&xs, try r.stream(&.{ "VBA", "dir" }));
    try zigstd.testExpectEqualSlices(u8, big[0..5000], try r.stream(&.{ "VBA", "Module1" }));
    try zigstd.testExpectEqualStrings("abc", try r.stream(&.{ "VBA", "_VBA_PROJECT" }));
}

test "duplicate names are rejected" {
    const nodes = [_]Node{
        .{ .name = "dir", .kind = .{ .stream = "a" } },
        .{ .name = "DIR", .kind = .{ .stream = "b" } },
    };
    try zigstd.testExpectError(error.CfbDuplicateName, write(zigstd.testAllocator, &nodes));
}
