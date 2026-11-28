//! Tests: TrueType glyphs loaded and hinted as FreeType 2.14.3 loads them,
//! against its own outlines (`refs.zig`), under `-Dtest-fonts`.

const std = @import("std");
const refs = @import("refs.zig");
const Sfnt = @import("../../sfnt/Sfnt.zig");
const truetype = @import("../truetype/__root__.zig");
const vary = @import("../vary/__root__.zig");
const test_fonts = @import("rstd_test_fonts");

const testing = std.testing;
const Instance = truetype.Instance;
const Loader = truetype.Loader;

/// How `refs.zig`'s mode loads: hinted or not, and the instance's options.
const Setup = struct { mode: refs.Mode, hinted: bool, version: truetype.Engine.Version = .v40, target: truetype.Engine.Target = .normal };

const setups = [_]Setup{
    .{ .mode = .unhinted, .hinted = false },
    .{ .mode = .native, .hinted = true },
    .{ .mode = .native_v35, .hinted = true, .version = .v35 },
    .{ .mode = .native_light, .hinted = true, .target = .light },
    .{ .mode = .native_mono, .hinted = true, .target = .mono },
    .{ .mode = .native_lcd, .hinted = true, .target = .lcd },
    .{ .mode = .native_lcd_v, .hinted = true, .target = .lcd_v },
};

fn same(g: Loader.Glyph, o: refs.Outline) bool {
    if (o.err != 0 or g.advance != o.advance) return false;
    if (!std.mem.eql(u16, g.contours, o.ends) or g.points.len != o.points.len) return false;
    for (g.points, o.points, g.tags, o.tags) |p, q, t, u| {
        if (p[0] != q[0] or p[1] != q[1] or t != @backingInt(u)) return false;
    }
    return true;
}

fn expectMode(comptime name: []const u8, sfnt: *const Sfnt, r: refs.Refs, blend: ?*const vary.Blend, setup: Setup) !void {
    const al = testing.allocator;
    var bad: usize = 0;
    var compared: usize = 0;
    var first: ?struct { cp: u21, size: u16 } = null;
    for (r.sizes, 0..) |ppem, s| {
        var inst = try Instance.init(al, sfnt, .{ .ppem = ppem, .version = setup.version, .target = setup.target, .blend = blend });
        defer inst.deinit();
        for (r.glyphs) |g| {
            const o = r.outline(g, s, setup.mode) orelse continue;
            var got = try Loader.load(al, &inst, g.gid, setup.hinted);
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
        std.debug.print("{s} {t}: {d} glyphs differ from FreeType, first U+{X:0>4} at {d} ppem\n", .{ name, setup.mode, bad, first.?.cp, first.?.size });
        return error.TestExpectedEqual;
    }
}

test "truetype glyphs: FreeType's outlines and advances, unhinted and hinted (-Dtest-fonts)" {
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
        for (setups) |setup| if (!r.cff) expectMode(f.name, &sfnt, r, if (b) |*x| x else null, setup) catch |e| switch (e) {
            error.TestExpectedEqual => failed = true,
            else => return e,
        };
    };
    if (!any) return error.SkipZigTest;
    if (failed) return error.TestExpectedEqual;
}
