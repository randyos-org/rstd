const std = @import("std");

const EdgeSegment = @import("EdgeSegment.zig");
const f64i = @import("math.zig").f64i;
const math = @import("math.zig");
const Shape = @import("Shape.zig");

const aux_types = @import("../sfnt/aux_types.zig");
const zigstd = @import("../../../zigstd.zig");
const Allocator = zigstd.Allocator;
const V2f64 = aux_types.V2f64;
const coloring_distance = @import("coloring_distance.zig");

pub const colorDistance = coloring_distance.colorDistance;

const edge_length_precision = 4;

pub const EdgeColor = enum(u8) {
    none = 0,
    red = 1,
    green = 2,
    blue = 3,
    red_green = 4,
    red_blue = 5,
    green_blue = 6,
    all = 7,

    pub fn init(rng: *std.Random.DefaultPrng) EdgeColor {
        const two_ch: [3]EdgeColor = .{ .red_green, .red_blue, .green_blue };
        return two_ch[rng.random().uintLessThan(usize, two_ch.len)];
    }

    fn cmyWithExcls(excl_1: EdgeColor, excl_2: EdgeColor) EdgeColor {
        const T = EdgeColor;
        const cmy_tuple = .{ T.red_green, T.red_blue, T.green_blue };

        var validation: [2]bool = @splat(false);
        for (cmy_tuple) |c| {
            if (c == excl_1) validation[0] = true;
            if (c == excl_2) validation[1] = true;
        }
        if (!std.meta.eql(validation, @splat(true)))
            @compileError("You must provide CMY exclusion values");

        for (cmy_tuple) |c| if (c != excl_1 and c != excl_2)
            return c;
    }

    fn colorSwitchTargets(self: EdgeColor) [2]EdgeColor {
        return switch (self) {
            .red_green => .{ .green_blue, .red_blue },
            .red_blue => .{ .red_green, .green_blue },
            .green_blue => .{ .red_blue, .red_green },
            else => @compileError("Switching is only allowed on two-channel edge colors"),
        };
    }

    pub fn change(self: *EdgeColor, banned: EdgeColor) void {
        switch (self.*) {
            inline .red_green, .red_blue, .green_blue => |c| switch (banned) {
                inline .red_green, .red_blue, .green_blue => |bc| if (comptime c != bc) {
                    self.* = comptime cmyWithExcls(c, bc);
                    return;
                },
                else => unreachable,
            },
            else => unreachable,
        }
    }

    pub fn random(self: *EdgeColor, rng: *std.Random.DefaultPrng) void {
        switch (self.*) {
            inline .red_green, .red_blue, .green_blue => |c| {
                const switch_targets = comptime colorSwitchTargets(c);
                self.* = switch_targets[rng.random().uintLessThan(usize, switch_targets.len)];
            },
            else => unreachable,
        }
    }

    pub fn hasChannel(self: EdgeColor, target_ch: EdgeColor) bool {
        return switch (target_ch) {
            .red => self == .all or self == .red_green or self == .red_blue or self == .red,
            .green => self == .all or self == .red_green or self == .green_blue or self == .green,
            .blue => self == .all or self == .red_blue or self == .green_blue or self == .blue,
            else => @panic("This function requires a single channel target"),
        };
    }
};

fn isCorner(a: V2f64, b: V2f64, cross_threshold: f64) bool {
    return math.dot(a, b) <= 0 or @abs(math.cross(a, b)) > cross_threshold;
}

pub fn findCorners(allocator: Allocator, contour: *const Shape.Contour, corners: *std.ArrayList(u32), cross_threshold: f64) !void {
    corners.clearRetainingCapacity();
    var prev_dir = math.normal(contour.edges.last().?.direction(1), true);
    for (contour.edges.items, 0..) |edge, i| {
        const edge_dir = math.normal(edge.direction(0), true);
        if (isCorner(prev_dir, edge_dir, cross_threshold))
            try corners.append(allocator, @intCast(i));

        prev_dir = math.normal(edge.direction(1), true);
    }
}

/// -1, 0 or 1 as `pos` falls in the first, middle or last third of `n`,
/// balanced so the results sum to zero (msdfgen 1c106ed).
pub fn symmetricalTrichotomy(pos: usize, n: usize) i32 {
    return @as(i32, @intCast((3 * pos + 1) / n)) - 1;
}

fn estimateEdgeLength(edge: *const EdgeSegment) f64 {
    var len: f64 = 0;
    var previous = edge.point(0);
    for (1..edge_length_precision) |i| {
        const current = edge.point(f64i(i) / @as(comptime_float, edge_length_precision));
        len += math.length(current - previous);
        previous = current;
    }
    return len;
}

/// colorSimple's and colorInkTrap's contour without corners: one random color.
fn colorNoCorners(contour: *Shape.Contour, color: *EdgeColor, rng: *std.Random.DefaultPrng) void {
    color.random(rng);
    for (contour.edges.items) |*edge| edge.color = color.*;
}

/// colorSimple's and colorInkTrap's contour with one corner at edge `corner`.
fn colorOneCorner(allocator: Allocator, contour: *Shape.Contour, corner: u32, rng: *std.Random.DefaultPrng) !void {
    const colors: [3]EdgeColor = .{ .init(rng), .all, .init(rng) };
    const corner_idx = 3 * corner;
    const edges_len = contour.edges.items.len;
    if (edges_len >= 3) {
        for (contour.edges.items, 0..) |*edge, i|
            edge.color = colors[@intCast(1 + symmetricalTrichotomy(i, edges_len))];
    } else if (edges_len >= 1) {
        var parts: [7]EdgeSegment = @splat(.{ .color = .none, .segment = undefined });

        @memcpy(parts[corner_idx..][0..3], &contour.edges.items[0].splitInThirds());
        if (edges_len >= 2) {
            @memcpy(parts[3 - corner_idx ..][0..3], &contour.edges.items[1].splitInThirds());
            for (0..6) |i| parts[i].color = colors[i / 2];
        } else for (0..3) |i| parts[i].color = colors[i];

        contour.edges.clearRetainingCapacity();
        var i: usize = 0;
        while (parts[i].color != .none) : (i += 1)
            try contour.edges.append(allocator, parts[i]);
    }
}

pub fn colorSimple(
    allocator: Allocator,
    rng_seed: u64,
    shape: *Shape,
    angle_threshold: f64,
) !void {
    const cross_threshold = @sin(angle_threshold);
    var rng: std.Random.DefaultPrng = .init(rng_seed);
    var color: EdgeColor = .init(&rng);

    var corners: std.ArrayList(u32) = .empty;
    defer corners.deinit(allocator);

    for (shape.contours.items) |*contour| {
        if (contour.edges.items.len == 0)
            continue;

        try findCorners(allocator, contour, &corners, cross_threshold);

        const corners_len = corners.items.len;
        switch (corners_len) {
            0 => colorNoCorners(contour, &color, &rng),
            1 => try colorOneCorner(allocator, contour, corners.items[0], &rng),
            else => {
                const start = corners.items[0];
                const edges_len = contour.edges.items.len;

                color.random(&rng);
                const initial_color = color;

                var spline: u32 = 0;
                for (0..edges_len) |i| {
                    const idx = (start + i) % edges_len;
                    if (spline + 1 < corners_len and corners.items[spline + 1] == idx) {
                        spline += 1;
                        if (spline == corners_len - 1)
                            color.change(initial_color)
                        else
                            color.random(&rng);
                    }
                    contour.edges.items[idx].color = color;
                }
            },
        }
    }
}

pub fn colorInkTrap(
    allocator: Allocator,
    rng_seed: u64,
    shape: *Shape,
    angle_threshold: f64,
) !void {
    const cross_threshold = @sin(angle_threshold);
    var rng: std.Random.DefaultPrng = .init(rng_seed);
    var color: EdgeColor = .init(&rng);

    var corners: std.ArrayList(struct {
        idx: u32,
        est_prev_edge_len: f64,
        minor: bool = false,
        color: EdgeColor = .none,
    }) = .empty;
    defer corners.deinit(allocator);

    for (shape.contours.items) |*contour| {
        if (contour.edges.items.len == 0) continue;

        var spline_length: f64 = 0.0;

        corners.clearRetainingCapacity();
        var prev_dir = math.normal(contour.edges.last().?.direction(1), true);
        for (contour.edges.items, 0..) |edge, i| {
            const edge_dir = math.normal(edge.direction(0), true);
            if (isCorner(prev_dir, edge_dir, cross_threshold)) {
                try corners.append(allocator, .{
                    .idx = @intCast(i),
                    .est_prev_edge_len = spline_length,
                });
                spline_length = 0;
            }

            spline_length += estimateEdgeLength(&edge);
            prev_dir = math.normal(edge.direction(1), true);
        }

        const corners_len = corners.items.len;
        switch (corners_len) {
            0 => colorNoCorners(contour, &color, &rng),
            1 => try colorOneCorner(allocator, contour, corners.items[0].idx, &rng),
            else => {
                var significant_corners = corners_len;
                if (corners_len > 3) {
                    corners.items[0].est_prev_edge_len += spline_length;
                    for (0..corners_len) |i| {
                        const current = corners.items[i].est_prev_edge_len;
                        const next = corners.items[(i + 1) % corners_len].est_prev_edge_len;
                        if (current > next and
                            next < corners.items[(i + 2) % corners_len].est_prev_edge_len)
                        {
                            corners.items[i].minor = true;
                            significant_corners -= 1;
                        }
                    }
                }

                var initial_color: EdgeColor = .none;
                for (0..corners_len) |i|
                    if (!corners.items[i].minor) {
                        significant_corners -= 1;
                        if (significant_corners == 0 or initial_color == .none)
                            color.random(&rng)
                        else
                            color.change(initial_color);
                        corners.items[i].color = color;
                        if (initial_color == .none)
                            initial_color = color;
                    };

                for (0..corners_len) |i| {
                    if (!corners.items[i].minor) {
                        color = corners.items[i].color;
                        continue;
                    }

                    color.change(corners.items[(i + 1) % corners_len].color);
                    corners.items[i].color = color;
                }

                const start = corners.items[0].idx;
                const edges_len = contour.edges.items.len;

                color = corners.items[0].color;

                var spline: u32 = 0;
                for (0..edges_len) |i| {
                    const idx = (start + i) % edges_len;
                    if (spline + 1 < corners_len and corners.items[spline + 1].idx == idx) {
                        spline += 1;
                        color = corners.items[spline].color;
                    }
                    contour.edges.items[idx].color = color;
                }
            },
        }
    }
}

test "symmetricalTrichotomy stays in -1..1, never decreases and balances" {
    for (1..33) |n| {
        var sum: i32 = 0;
        var prev: i32 = -1;
        for (0..n) |pos| {
            const t = symmetricalTrichotomy(pos, n);
            try zigstd.testExpectTrue(t >= prev and t <= 1);
            prev = t;
            sum += t;
        }
        try zigstd.testExpectEqual(@as(i32, 0), sum);
    }
}
