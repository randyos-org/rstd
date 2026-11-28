//! Tests: CFF glyphs loaded and hinted as FreeType 2.14.3 loads them,
//! against its own outlines (`refs.zig`), under `-Dtest-fonts`.

const std = @import("std");
const refs = @import("refs.zig");
const Sfnt = @import("../../sfnt/Sfnt.zig");
const cff = @import("../cff/__root__.zig");
const vary = @import("../vary/__root__.zig");
const Glyph = @import("../Glyph.zig");
const test_fonts = @import("rstd_test_fonts");

const testing = std.testing;

/// FreeType's CFF engine ignores the light target.
const setups = [_]struct { mode: refs.Mode, hinted: bool }{
    .{ .mode = .unhinted, .hinted = false },
    .{ .mode = .native, .hinted = true },
    .{ .mode = .native_light, .hinted = true },
    .{ .mode = .native_mono, .hinted = true },
    .{ .mode = .native_lcd, .hinted = true },
    .{ .mode = .native_lcd_v, .hinted = true },
};

fn same(g: Glyph, o: refs.Outline) bool {
    if (o.err != 0 or g.advance != o.advance) return false;
    if (!std.mem.eql(u16, g.contours, o.ends) or g.points.len != o.points.len) return false;
    for (g.points, o.points, g.tags, o.tags) |p, q, t, u| {
        if (p[0] != q[0] or p[1] != q[1] or t != @backingInt(u)) return false;
    }
    return true;
}

fn expectMode(comptime name: []const u8, sfnt: *const Sfnt, r: refs.Refs, blend: ?*const vary.Blend, mode: refs.Mode, hinted: bool) !void {
    const al = testing.allocator;
    var bad: usize = 0;
    var compared: usize = 0;
    var first: ?struct { cp: u21, size: u16 } = null;
    for (r.sizes, 0..) |ppem, s| {
        var inst: cff.Instance = try .init(sfnt, ppem);
        inst.blend = blend;
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
    inline for (refs.files) |f| if (@field(test_fonts, f.font)) |bytes| {
        any = true;
        const r = try refs.parse(arena.allocator(), f.gz);
        const sfnt = try Sfnt.init(arena.allocator(), bytes, .{});
        var b = try refs.blend(testing.allocator, &sfnt, r);
        defer if (b) |*x| x.deinit();
        for (setups) |s| if (r.cff) expectMode(f.name, &sfnt, r, if (b) |*x| x else null, s.mode, s.hinted) catch |e| switch (e) {
            error.TestExpectedEqual => failed = true,
            else => return e,
        };
    };
    if (!any) return error.SkipZigTest;
    if (failed) return error.TestExpectedEqual;
}
