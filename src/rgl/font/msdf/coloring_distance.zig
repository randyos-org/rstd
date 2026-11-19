//! Edge coloring by the distances between splines, as a graph coloring.

const std = @import("std");
const coloring = @import("coloring.zig");
const EdgeSegment = @import("EdgeSegment.zig");
const math = @import("math.zig");
const Shape = @import("Shape.zig");
const zigstd = @import("../../../zigstd.zig");
const f64i = math.f64i;
const Allocator = zigstd.Allocator;
const EdgeColor = coloring.EdgeColor;
const findCorners = coloring.findCorners;
const symmetricalTrichotomy = coloring.symmetricalTrichotomy;

const max_distance_recolors = 32;
const edge_distance_precision = 32;

const ColorSelection = packed struct(u3) {
    red_green: bool = false,
    red_blue: bool = false,
    green_blue: bool = false,

    pub const all: ColorSelection = .{
        .red_green = true,
        .red_blue = true,
        .green_blue = true,
    };

    pub fn set(self: *ColorSelection, color: EdgeColor) void {
        switch (color) {
            .red_green => self.red_green = true,
            .red_blue => self.red_blue = true,
            .green_blue => self.green_blue = true,
            else => {},
        }
    }
};

fn splineToSplineDistance(
    edge_segments: []const *EdgeSegment,
    a_from: u32,
    a_to: u32,
    b_from: u32,
    b_to: u32,
) f64 {
    var min_dist = std.math.floatMax(f64);
    for (a_from..a_to) |a_idx| for (b_from..b_to) |b_idx| {
        const a = edge_segments[a_idx];
        const b = edge_segments[b_idx];

        const eql = std.meta.eql;
        if (eql(a.point(0), b.point(0)) or
            eql(a.point(0), b.point(1)) or
            eql(a.point(1), b.point(0)) or
            eql(a.point(1), b.point(1)))
            return 0.0;

        min_dist = @min(min_dist, math.length(b.point(0) - a.point(0)));
        for (0..edge_distance_precision) |i| {
            const t = f64i(i) / @as(comptime_float, edge_distance_precision);
            const ab_dist, _ = a.signedDistance(b.point(t));
            const ba_dist, _ = b.signedDistance(a.point(t));
            min_dist = @min(min_dist, @abs(ab_dist.distance));
            min_dist = @min(min_dist, @abs(ba_dist.distance));
        }
    };

    return min_dist;
}

fn colorGraph(colors: []EdgeColor, edge_matrix: []const []const bool, rng: *std.Random.DefaultPrng) void {
    for (colors, 0..) |*color, i| {
        var banned_colors: ColorSelection = .{};

        for (0..i) |j| if (edge_matrix[i][j])
            banned_colors.set(colors[j]);

        switch (banned_colors) {
            .{} => {
                const rem_colors: [3]EdgeColor = .{ .red_green, .red_blue, .green_blue };
                color.* = rem_colors[rng.random().uintLessThan(usize, rem_colors.len)];
            },
            .{ .red_green = true } => {
                const rem_colors: [2]EdgeColor = .{ .red_blue, .green_blue };
                color.* = rem_colors[rng.random().uintLessThan(usize, rem_colors.len)];
            },
            .{ .red_blue = true } => {
                const rem_colors: [2]EdgeColor = .{ .red_green, .green_blue };
                color.* = rem_colors[rng.random().uintLessThan(usize, rem_colors.len)];
            },
            .{ .green_blue = true } => {
                const rem_colors: [2]EdgeColor = .{ .red_green, .red_blue };
                color.* = rem_colors[rng.random().uintLessThan(usize, rem_colors.len)];
            },
            .{ .red_green = true, .red_blue = true } => color.* = .green_blue,
            .{ .red_green = true, .green_blue = true } => color.* = .red_blue,
            .{ .red_blue = true, .green_blue = true } => color.* = .red_green,
            .all => color.* = .red_green,
        }
    }
}

fn vertexPossibleColors(colors: []EdgeColor, edge_vec: []const bool) ?EdgeColor {
    var banned_colors: ColorSelection = .{};
    for (colors, edge_vec) |color, edge|
        if (edge) banned_colors.set(color);

    return switch (banned_colors) {
        .{} => .red_green,
        .{ .red_green = true } => .red_blue,
        .{ .red_blue = true } => .green_blue,
        .{ .green_blue = true } => .red_green,
        .{ .red_green = true, .red_blue = true } => .green_blue,
        .{ .red_green = true, .green_blue = true } => .red_blue,
        .{ .red_blue = true, .green_blue = true } => .red_green,
        .all => null,
    };
}

fn uncolorSameNeighbors(
    allocator: Allocator,
    uncolored: *std.Deque(u32),
    colors: []EdgeColor,
    edge_matrix: []const []const bool,
    idx: usize,
) !void {
    const size = colors.len;
    for (idx + 1..size) |i|
        if (edge_matrix[idx][i] and colors[i] == colors[idx]) {
            colors[i] = .none;
            try uncolored.pushFront(allocator, @intCast(i));
        };

    for (0..idx) |i|
        if (edge_matrix[idx][i] and colors[i] == colors[idx]) {
            colors[i] = .none;
            try uncolored.pushFront(allocator, @intCast(i));
        };
}

fn addEdge(allocator: Allocator, front_buf: []EdgeColor, back_buf: []EdgeColor, edge_matrix: [][]bool, a: usize, b: usize) !bool {
    edge_matrix[a][b] = true;
    edge_matrix[b][a] = true;

    if (front_buf[a] != front_buf[b])
        return true;

    if (vertexPossibleColors(front_buf, edge_matrix[b])) |color| {
        front_buf[b] = color;
        return true;
    }

    @memcpy(back_buf, front_buf);

    back_buf[b] = switch (back_buf[a]) {
        .red_green => .red_blue,
        .red_blue => .green_blue,
        .green_blue => .red_green,
        else => unreachable,
    };

    var uncolored: std.Deque(u32) = .empty;
    defer uncolored.deinit(allocator);

    try uncolorSameNeighbors(allocator, &uncolored, back_buf, edge_matrix, b);

    var step: u16 = 0;
    while (step < max_distance_recolors) {
        const i = uncolored.popFront() orelse break;
        if (vertexPossibleColors(back_buf, edge_matrix[i])) |color| {
            back_buf[i] = color;
            continue;
        }

        step += 1;
        back_buf[i] = switch (step % 3) {
            0 => .red_green,
            1 => .red_blue,
            2 => .green_blue,
            else => unreachable,
        };

        while (edge_matrix[i][a] and
            back_buf[i] == back_buf[a]) : (step += 1)
            back_buf[i] = switch (step % 3) {
                0 => .red_green,
                1 => .red_blue,
                2 => .green_blue,
                else => unreachable,
            };

        try uncolorSameNeighbors(allocator, &uncolored, back_buf, edge_matrix, i);
    }

    if (uncolored.len > 0) {
        edge_matrix[a][b] = false;
        edge_matrix[b][a] = false;
        return false;
    }

    @memcpy(front_buf, back_buf);
    return true;
}

pub fn colorDistance(
    allocator: Allocator,
    rng_seed: u64,
    shape: *Shape,
    angle_threshold: f64,
) !void {
    const cross_threshold = @sin(angle_threshold);
    var rng: std.Random.DefaultPrng = .init(rng_seed);

    var edge_segments: std.ArrayList(*EdgeSegment) = .empty;
    defer edge_segments.deinit(allocator);

    var spline_starts: std.ArrayList(u32) = .empty;
    defer spline_starts.deinit(allocator);

    var corners: std.ArrayList(u32) = .empty;
    defer corners.deinit(allocator);

    for (shape.contours.items) |*contour| {
        if (contour.edges.items.len == 0)
            continue;

        try findCorners(allocator, contour, &corners, cross_threshold);
        try spline_starts.append(allocator, @intCast(edge_segments.items.len));

        const corners_len = corners.items.len;
        switch (corners_len) {
            0 => {
                for (contour.edges.items) |*edge|
                    try edge_segments.append(allocator, edge);
            },
            1 => {
                const corner = corners.items[0];
                const edges_len = contour.edges.items.len;
                if (edges_len >= 3) {
                    for (0..edges_len) |i| {
                        if (i == edges_len / 2)
                            try spline_starts.append(allocator, @intCast(edge_segments.items.len));

                        const edge = &contour.edges.items[(corner + i) % edges_len];
                        if (symmetricalTrichotomy(i, edges_len) != 0)
                            try edge_segments.append(allocator, edge)
                        else
                            edge.color = .all;
                    }
                } else if (edges_len >= 1) {
                    var parts: [7]EdgeSegment = @splat(.{ .color = .none, .segment = undefined });
                    var append: [7]bool = @splat(false);

                    const corner_idx = 3 * corner;
                    @memcpy(parts[corner_idx..][0..3], &contour.edges.items[0].splitInThirds());
                    if (edges_len >= 2) {
                        @memcpy(parts[3 - corner_idx ..][0..3], &contour.edges.items[1].splitInThirds());
                        append[0] = true;
                        append[1] = true;
                        parts[2].color = .all;
                        parts[3].color = .all;
                        try spline_starts.append(allocator, @intCast(edge_segments.items.len + 2));
                        append[4] = true;
                        append[5] = true;
                    } else {
                        append[0] = true;
                        parts[1].color = .all;
                        try spline_starts.append(allocator, @intCast(edge_segments.items.len + 1));
                        append[2] = true;
                    }

                    contour.edges.clearRetainingCapacity();
                    var i: usize = 0;
                    while (parts[i].color != .none) : (i += 1) {
                        try contour.edges.append(allocator, parts[i]);
                        if (append[i])
                            try edge_segments.append(allocator, &contour.edges.items[contour.edges.items.len - 1]);
                    }
                }
            },
            else => {
                const start = corners.items[0];
                const edges_len = contour.edges.items.len;

                var spline: u32 = 0;
                for (0..edges_len) |i| {
                    const idx = (start + i) % edges_len;
                    if (spline + 1 < corners_len and corners.items[spline + 1] == idx) {
                        try spline_starts.append(allocator, @intCast(edge_segments.items.len));
                        spline += 1;
                    }
                    try edge_segments.append(allocator, &contour.edges.items[idx]);
                }
            },
        }
    }
    try spline_starts.append(allocator, @intCast(edge_segments.items.len));

    const spline_count = spline_starts.items.len - 1;
    if (spline_count == 0) return;

    const distance_matrix = try allocator.alloc(f64, spline_count * spline_count);
    defer allocator.free(distance_matrix);
    @memset(distance_matrix, 0.0);

    for (0..spline_count) |i| {
        distance_matrix[i * spline_count + i] = -1.0;
        for (i + 1..spline_count) |j| {
            const dist = splineToSplineDistance(
                edge_segments.items,
                spline_starts.items[i],
                spline_starts.items[i + 1],
                spline_starts.items[j],
                spline_starts.items[j + 1],
            );
            distance_matrix[i * spline_count + j] = dist;
            distance_matrix[j * spline_count + i] = dist;
        }
    }

    const graph_edge_distances = try allocator.alloc(u32, spline_count * (spline_count - 1) / 2);
    defer allocator.free(graph_edge_distances);
    @memset(graph_edge_distances, 0);

    var graph_idx: usize = 0;
    for (0..spline_count) |i| for (i + 1..spline_count) |j| {
        graph_edge_distances[graph_idx] = @intCast(i * spline_count + j);
        graph_idx += 1;
    };

    const SortCtx = struct {
        distance_matrix: []const f64,

        pub fn lessThan(self: @This(), a: u32, b: u32) bool {
            return self.distance_matrix[a] < self.distance_matrix[b];
        }
    };
    std.sort.pdq(
        u32,
        graph_edge_distances[0..graph_idx],
        SortCtx{ .distance_matrix = distance_matrix },
        SortCtx.lessThan,
    );

    const edge_matrix_buf = try allocator.alloc(bool, spline_count * spline_count);
    defer allocator.free(edge_matrix_buf);
    @memset(edge_matrix_buf, false);

    const edge_matrix = try allocator.alloc([]bool, spline_count);
    defer allocator.free(edge_matrix);
    for (edge_matrix, 0..) |*vec, i|
        vec.* = edge_matrix_buf[i * spline_count ..][0..spline_count];

    var next_edge: usize = 0;
    while (next_edge < graph_idx) : (next_edge += 1) {
        const dist_idx = graph_edge_distances[next_edge];
        if (distance_matrix[dist_idx] != 0.0) break;

        const row = dist_idx / spline_count;
        const col = dist_idx % spline_count;
        edge_matrix[row][col] = true;
        edge_matrix[col][row] = true;
    }

    const colors = try allocator.alloc(EdgeColor, spline_count * 2);
    defer allocator.free(colors);
    @memset(colors, .none);

    const front_buf = colors[0..spline_count];
    const back_buf = colors[spline_count..][0..spline_count];

    colorGraph(front_buf, edge_matrix, &rng);
    while (next_edge < graph_idx) : (next_edge += 1) {
        const dist_idx = graph_edge_distances[next_edge];
        const row = dist_idx / spline_count;
        const col = dist_idx % spline_count;
        _ = try addEdge(
            allocator,
            front_buf,
            back_buf,
            edge_matrix,
            row,
            col,
        );
    }

    var spline: usize = 0;
    for (spline_starts.items[0]..edge_segments.items.len) |i| {
        if (spline_starts.items[spline + 1] == i)
            spline += 1;

        const new_color = front_buf[spline];
        if (new_color != .none) edge_segments.items[i].color = new_color;
    }
}
