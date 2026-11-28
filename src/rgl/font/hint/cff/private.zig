//! The Private DICT values the CFF hinter and loader use, parsed as FreeType
//! 2.14.3 parses them (`src/cff/cffparse.c`, the `cfftoken.h` table,
//! `BlueScale` scaled by 1000) through `dict.zig`, after fontations
//! `read-fonts/src/ps/hinting.rs` (commit 4f115482, MIT). Values are 16.16.

const std = @import("std");
const dict = @import("dict.zig");
const VarStore = @import("VarStore.zig");

pub const Error = dict.Error;

/// Up to seven zones, each a bottom and a top edge (`Blues`).
pub const Blues = struct {
    values: [7][2]i32 = undefined,
    len: usize = 0,

    pub fn slice(b: *const Blues) []const [2]i32 {
        return b.values[0..b.len];
    }

    /// Pairs from `v`, which is already summed from deltas.
    fn init(v: []const i32) Blues {
        var b: Blues = .{};
        var i: usize = 0;
        while (i + 1 < v.len and b.len < 7) : (i += 2) {
            b.values[b.len] = .{ v[i], v[i + 1] };
            b.len += 1;
        }
        return b;
    }
};

pub const Private = struct {
    blues: Blues = .{},
    family_blues: Blues = .{},
    other_blues: Blues = .{},
    family_other_blues: Blues = .{},
    /// 0.039625, the spec's default.
    blue_scale: i32 = fromF64(0.039625),
    blue_shift: i32 = 7 << 16,
    blue_fuzz: i32 = 1 << 16,
    language_group: i32 = 0,
    /// The local Subrs' offset from the Private DICT's start.
    subrs: ?usize = null,
    default_width: ?i32 = null,
    nominal_width: i32 = 0,
    /// CFF2: the item variation data its blends, and its charstrings', start from.
    vsindex: usize = 0,
};

fn fromF64(comptime v: f64) i32 {
    return @intFromFloat(@round(v * 65536.0));
}

/// The DICT `data`; a CFF2 one blends with `vectors`.
pub fn parse(data: []const u8, vectors: ?*VarStore.Vectors) Error!Private {
    var p: Private = .{};
    var parser: dict.Parser = .{ .data = data, .vectors = vectors };
    while (try parser.next()) |item| {
        const args = item.args;
        const op = item.op;
        switch (op) {
            6, 7, 8, 9 => {
                var vals: [14]i32 = undefined;
                var sum: i32 = 0;
                // FreeType keeps 14 numbers of the blues, 10 of the others
                const k = @min(args.len, @as(usize, if (op == 7 or op == 9) 10 else 14));
                for (args[0..k], 0..) |a, i| {
                    sum +%= a.fixed();
                    vals[i] = sum;
                }
                const blues = Blues.init(vals[0..k]);
                switch (op) {
                    6 => p.blues = blues,
                    7 => p.other_blues = blues,
                    8 => p.family_blues = blues,
                    else => p.family_other_blues = blues,
                }
            },
            0x0C09 => if (args.len > 0) {
                const a = args[args.len - 1];
                p.blue_scale = if (a.bcd) |b| b.value(true) else a.fixed();
            },
            0x0C0A => if (args.len > 0) {
                p.blue_shift = args[args.len - 1].fixed() & ~@as(i32, 0xFFFF);
            },
            0x0C0B => if (args.len > 0) {
                p.blue_fuzz = args[args.len - 1].fixed() & ~@as(i32, 0xFFFF);
            },
            0x0C11 => if (args.len > 0) {
                p.language_group = args[args.len - 1].fixed() >> 16;
            },
            19 => if (args.len > 0) {
                const v = args[args.len - 1].fixed() >> 16;
                if (v >= 0) p.subrs = @intCast(v);
            },
            20 => if (args.len > 0) {
                p.default_width = args[args.len - 1].fixed() & ~@as(i32, 0xFFFF);
            },
            21 => if (args.len > 0) {
                p.nominal_width = args[args.len - 1].fixed() & ~@as(i32, 0xFFFF);
            },
            else => {},
        }
    }
    p.vsindex = parser.vsindex;
    return p;
}

test "cff private: a BlueScale real, deltas into zones, widths floored" {
    // BlueValues -20 20 480 20 (deltas), BlueScale 0.039625 as BCD, nominalWidthX 500.5
    const data = [_]u8{ 139 - 20, 139 + 20, 248, 116, 139 + 20, 6, 30, 0x0A, 0x03, 0x96, 0x25, 0xFF, 12, 9, 30, 0x50, 0x0A, 0x5F, 21 };
    const p = try parse(&data, null);
    try std.testing.expectEqual(@as(usize, 2), p.blues.len);
    try std.testing.expectEqual([2]i32{ -20 << 16, 0 }, p.blues.values[0]);
    try std.testing.expectEqual([2]i32{ 480 << 16, 500 << 16 }, p.blues.values[1]);
    try std.testing.expect(@abs(p.blue_scale - fromF64(0.039625)) <= 1);
    try std.testing.expectEqual(@as(i32, 500 << 16), p.nominal_width);
}
