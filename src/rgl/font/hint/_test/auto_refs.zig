//! Tests: glyphs autohinted as FreeType 2.14.3 autohints them
//! (`FT_LOAD_FORCE_AUTOHINT`, normal and light), against its own outlines
//! (`refs.zig`), under `-Dtest-fonts`.

const std = @import("std");
const refs = @import("refs.zig");
const Sfnt = @import("../../sfnt/Sfnt.zig");
const Hinter = @import("../auto/Hinter.zig");
const Glyph = @import("../Glyph.zig");
const test_fonts = @import("rstd_test_fonts");

const testing = std.testing;

const setups = [_]struct { mode: refs.Mode, render: Hinter.RenderMode }{
    .{ .mode = .autohint, .render = .normal },
    .{ .mode = .autohint_light, .render = .light },
};

fn same(g: Glyph, o: refs.Outline) bool {
    if (o.err != 0 or g.advance != o.advance or g.lsb_delta != o.lsb_delta or g.rsb_delta != o.rsb_delta) return false;
    if (!std.mem.eql(u16, g.contours, o.ends) or g.points.len != o.points.len) return false;
    for (g.points, o.points, g.tags, o.tags) |p, q, t, u| {
        if (p[0] != q[0] or p[1] != q[1] or t != @backingInt(u)) return false;
    }
    return true;
}

fn expectMode(comptime name: []const u8, hinter: *Hinter, r: refs.Refs, mode: refs.Mode, render: Hinter.RenderMode) !void {
    const al = testing.allocator;
    var bad: usize = 0;
    var compared: usize = 0;
    var first: ?struct { cp: u21, size: u16 } = null;
    for (r.sizes, 0..) |ppem, s| {
        for (r.glyphs) |g| {
            const o = r.outline(g, s, mode) orelse continue;
            var got = try hinter.load(al, g.gid, ppem, render);
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

test "autohinted glyphs: FreeType's outlines, advances and side-bearing deltas (-Dtest-fonts)" {
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
        var hinter = try Hinter.init(testing.allocator, &sfnt, if (b) |*x| x else null);
        defer hinter.deinit();
        for (setups) |s| expectMode(f.name, &hinter, r, s.mode, s.render) catch |e| switch (e) {
            error.TestExpectedEqual => failed = true,
            else => return e,
        };
    };
    if (!any) return error.SkipZigTest;
    if (failed) return error.TestExpectedEqual;
}
