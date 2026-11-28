//! Tests: CFF2 blends (`VarStore.zig`, `dict.zig`, `charstring.zig`) on a
//! hand-made variation store of one axis and two regions.

const std = @import("std");
const VarStore = @import("../VarStore.zig");
const dict = @import("../dict.zig");
const private = @import("../private.zig");
const charstring = @import("../charstring.zig");
const Index = @import("../../../sfnt/tables/cff/Index.zig");

const testing = std.testing;

/// A 16-bit length, then a store: regions 0..1 peaking at 1 and -1..0
/// peaking at -1; data 0 uses both, data 1 only the first.
const store_bytes = [_]u8{
    0, 0, // length, unread
    0, 1, 0, 0, 0, 16, 0, 2, // format 1, region list at 16, two data
    0, 0, 0, 32, // data 0 at 32
    0, 0, 0, 42, // data 1 at 42
    0, 1, 0, 2, // one axis, two regions
    0x00, 0x00, 0x40, 0x00, 0x40, 0x00, // 0, 1, 1
    0xC0, 0x00, 0xC0, 0x00, 0x00, 0x00, // -1, -1, 0
    0, 0, 0, 0, 0, 2, 0, 0, 0, 1, // data 0: no items, regions 0 and 1
    0, 0, 0, 0, 0, 1, 0, 0, // data 1: region 0
};

test "cff2 blend vectors: one per region, the default instance's 1, 0, ..." {
    const al = testing.allocator;
    var vs = try VarStore.init(al, &store_bytes, 0);
    defer vs.deinit(al);
    // halfway to the maximum: the first region at half, the second off
    const half = try vs.blendVector(al, 0, &.{0x8000});
    defer al.free(half);
    try testing.expectEqualSlices(i32, &.{ 0x10000, 0x8000, 0 }, half);
    const def = try vs.blendVector(al, 1, &.{});
    defer al.free(def);
    try testing.expectEqualSlices(i32, &.{ 0x10000, 0 }, def);
    try testing.expectError(error.InvalidVariationStore, vs.blendVector(al, 2, &.{0}));
}

test "cff2 Private DICT: blended blue zones, rounded as cff_parse_num rounds" {
    const al = testing.allocator;
    var vs = try VarStore.init(al, &store_bytes, 0);
    defer vs.deinit(al);
    var vectors: VarStore.Vectors = .{ .allocator = al, .store = &vs, .coords = &.{-0x4000} };
    defer vectors.deinit();
    // BlueValues -10 blend(+0 +8) 20 blend(+0 -6), with the -0.25 vector 1, 0, 0.25:
    // -10 + 8 / 4 = -8 and 20 - 6 / 4 = 18.5, which rounds to 19
    const data = [_]u8{ 139 - 10, 139 + 20, 139, 139 + 8, 139, 139 - 6, 139 + 2, 23, 6 };
    const p = try private.parse(&data, &vectors);
    try testing.expectEqual(@as(usize, 1), p.blues.len);
    try testing.expectEqual([2]i32{ -8 << 16, 11 << 16 }, p.blues.values[0]);
}

/// A sink that keeps the points it is given.
const Points = struct {
    pts: [8][2]i32 = undefined,
    n: usize = 0,

    fn add(s: *Points, x: i32, y: i32) void {
        s.pts[s.n] = .{ x, y };
        s.n += 1;
    }
    pub fn moveTo(s: *Points, x: i32, y: i32) void {
        s.add(x, y);
    }
    pub fn lineTo(s: *Points, x: i32, y: i32) void {
        s.add(x, y);
    }
    pub fn curveTo(s: *Points, x1: i32, y1: i32, x2: i32, y2: i32, x3: i32, y3: i32) void {
        _ = .{ x1, y1, x2, y2 };
        s.add(x3, y3);
    }
    pub fn close(_: *Points) void {}
    pub fn hstem(_: *Points, _: i32, _: i32) void {}
    pub fn vstem(_: *Points, _: i32, _: i32) void {}
    pub fn hintMask(_: *Points, _: []const u8) void {}
    pub fn counterMask(_: *Points, _: []const u8) void {}
    pub fn closeOpenPath(_: *Points) void {}
};

test "cff2 charstrings: blend without a width, vsindex, and no endchar" {
    const al = testing.allocator;
    var vs = try VarStore.init(al, &store_bytes, 0);
    defer vs.deinit(al);
    var vectors: VarStore.Vectors = .{ .allocator = al, .store = &vs, .coords = &.{0x8000} };
    defer vectors.deinit();
    const empty: Index = try .init(&.{ 0, 0 }, true);
    const ctx: charstring.Context = .{ .global = empty, .local = null, .cff2 = .{ .vectors = &vectors, .vsindex = 0, .maxstack = 513 } };
    // vsindex 1 (one region), 100 blend(+40) 0 rmoveto: at half, 120; then
    // 10 blend(+20) 0 rlineto: 20 more; an endchar only clears the stack
    const code = [_]u8{ 139 + 1, 15, 139 + 100, 139 + 40, 139 + 1, 16, 139, 21, 139 + 10, 139 + 20, 139 + 1, 16, 139, 5, 14 };
    var sink: Points = .{};
    const r = try charstring.evaluate(ctx, &code, &sink, .{ 0, 0 });
    try testing.expectEqual(@as(?i32, null), r.width);
    try testing.expectEqualSlices([2]i32, &.{ .{ 120 << 16, 0 }, .{ 140 << 16, 0 } }, sink.pts[0..sink.n]);
    // vsindex after a blend fails
    var again: Points = .{};
    try testing.expectError(error.InvalidCharstring, charstring.evaluate(ctx, &[_]u8{ 139 + 10, 139 + 20, 139, 139 + 1, 16, 139 + 1, 15 }, &again, .{ 0, 0 }));
}
