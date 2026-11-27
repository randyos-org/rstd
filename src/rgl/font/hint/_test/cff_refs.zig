//! Tests: CFF glyphs loaded and hinted as FreeType 2.14.3 loads them,
//! against its own outlines (`refs.zig`), under `-Dtest-fonts`.

const std = @import("std");
const refs = @import("refs.zig");
const Sfnt = @import("../../sfnt/Sfnt.zig");
const cff = @import("../cff/__root__.zig");
const Glyph = @import("../Glyph.zig");
const test_fonts = @import("rstd_test_fonts");

const testing = std.testing;

const fonts = .{
    .{ "commit_mono", @embedFile("refs/commit_mono.bin.gz") },
    .{ "atkynson_mono", @embedFile("refs/atkynson_mono.bin.gz") },
};

/// FreeType's CFF engine ignores the light target.
const setups = [_]struct { mode: refs.Mode, hinted: bool }{
    .{ .mode = .unhinted, .hinted = false },
    .{ .mode = .native, .hinted = true },
    .{ .mode = .native_light, .hinted = true },
};

fn same(g: Glyph, o: refs.Outline) bool {
    if (o.err != 0 or g.advance != o.advance) return false;
    if (!std.mem.eql(u16, g.contours, o.ends) or g.points.len != o.points.len) return false;
    for (g.points, o.points, g.tags, o.tags) |p, q, t, u| {
        if (p[0] != q[0] or p[1] != q[1] or t != @backingInt(u)) return false;
    }
    return true;
}

fn expectMode(comptime name: []const u8, sfnt: *const Sfnt, r: refs.Refs, mode: refs.Mode, hinted: bool) !void {
    const al = testing.allocator;
    var bad: usize = 0;
    var compared: usize = 0;
    var first: ?struct { cp: u21, size: u16 } = null;
    for (r.sizes, 0..) |ppem, s| {
        const inst: cff.Instance = try .init(sfnt, ppem);
        for (r.glyphs) |g| {
            const o = r.outline(g, s, mode) orelse continue;
            var got = try cff.load(al, &inst, @fromBackingInt(g.gid), hinted);
            defer got.deinit(al);
            compared += 1;
            if (!same(got, o)) {
                bad += 1;
                if (first == null) first = .{ .cp = g.codepoint, .size = ppem };
            }
        }
    }
    try testing.expect(compared > 0);
    if (bad != 0) {
        std.debug.print("{s} {t}: {d} of {d} glyphs differ from FreeType, first U+{X:0>4} at {d} ppem\n", .{ name, mode, bad, compared, first.?.cp, first.?.size });
        return error.TestExpectedEqual;
    }
}

test "cff glyphs: FreeType's outlines and advances, unhinted and hinted (-Dtest-fonts)" {
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    var any = false;
    var failed = false;
    inline for (fonts) |f| if (@field(test_fonts, f[0])) |bytes| {
        any = true;
        const r = try refs.parse(arena.allocator(), f[1]);
        const sfnt = try Sfnt.init(arena.allocator(), bytes, .{});
        for (setups) |s| expectMode(f[0], &sfnt, r, s.mode, s.hinted) catch |e| switch (e) {
            error.TestExpectedEqual => failed = true,
            else => return e,
        };
    };
    if (!any) return error.SkipZigTest;
    if (failed) return error.TestExpectedEqual;
}
