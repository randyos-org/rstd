//! Distance field pixels of a shape: the distances, their sign correction and the output bytes.

const std = @import("std");
const zigstd = @import("../../../zigstd.zig");
const aux_types = @import("../sfnt/aux_types.zig");
const coloring = @import("coloring.zig");
const EdgeSegment = @import("EdgeSegment.zig");
const error_correction = @import("error_correction.zig");
const math = @import("math.zig");
const Scanline = @import("Scanline.zig");
const Shape = @import("Shape.zig");
const types = @import("types.zig");
const Image = @import("../../image.zig").Image;
const Allocator = zigstd.Allocator;
const EdgeColor = coloring.EdgeColor;
const V2f64 = aux_types.V2f64;
const f64i = math.f64i;
const SdfType = types.SdfType;
const Options = types.Options;
const Msdf10Pixel = types.Msdf10Pixel;
const DistPixel = types.DistPixel;
const OutPixel = types.OutPixel;

pub fn getSdfPixels(
    allocator: Allocator,
    opts: *const Options,
    w: u16,
    h: u16,
    shape: *Shape,
    tfm: V2f64,
) ![]const u8 {
    const px_size = f64i(opts.px_size);
    const px_range = f64i(opts.px_range) / px_size;

    switch (opts.sdf_type) {
        inline else => |ty| {
            const n = comptime ty.numChannels();
            var dist: Image(DistPixel(ty)) = try .init(allocator, w, h, null);
            defer dist.deinit(allocator);

            if (ty.requiresColoring()) switch (opts.coloring_method) {
                .simple => try coloring.colorSimple(allocator, opts.coloring_rng_seed, shape, opts.corner_angle_threshold),
                .ink_trap => try coloring.colorInkTrap(allocator, opts.coloring_rng_seed, shape, opts.corner_angle_threshold),
                .distance => try coloring.colorDistance(allocator, opts.coloring_rng_seed, shape, opts.corner_angle_threshold),
            };
            generate(ty, dist, px_size, shape.*, px_range, tfm);

            if (!opts.orient_contours)
                if (opts.scanline_fill_rule) |fill_rule| {
                    if (comptime n == 1)
                        try sdfSignCorrection(allocator, dist, px_size, shape.*, tfm, fill_rule)
                    else
                        try msdfSignCorrection(n, allocator, dist, px_size, shape.*, tfm, fill_rule);
                };

            if (comptime ty.requiresColoring())
                if (opts.error_correction_opts) |ec_opts|
                    try error_correction.correct(n, allocator, shape, ec_opts, opts.scanline_fill_rule != null, px_size, px_range, tfm, dist);

            const out: Image(OutPixel(ty)) = try .init(allocator, w, h, null);
            for (out.data, dist.data) |*o, d| {
                if (comptime ty == .msdf10)
                    o.* = std.mem.toBytes(Msdf10Pixel{
                        .r = @trunc(std.math.maxInt(u10) * std.math.clamp(d[0], 0.0, 1.0)),
                        .g = @trunc(std.math.maxInt(u10) * std.math.clamp(d[1], 0.0, 1.0)),
                        .b = @trunc(std.math.maxInt(u10) * std.math.clamp(d[2], 0.0, 1.0)),
                        .a = std.math.maxInt(u2),
                    })
                else for (o, d) |*v, c|
                    v.* = @trunc(std.math.maxInt(u8) * std.math.clamp(c, 0.0, 1.0));
            }
            return std.mem.sliceAsBytes(out.data);
        },
    }
}

fn sdfSignCorrection(
    allocator: Allocator,
    dist: Image([1]f64),
    scale: f64,
    shape: Shape,
    tfm: V2f64,
    fill_rule: Scanline.FillRule,
) !void {
    var scanline: Scanline = .{};
    defer scanline.intersections.deinit(allocator);
    for (0..dist.height) |y| {
        const row: u32 = @intCast(dist.height - y - 1);
        try shape.scanline(&scanline, (f64i(y) + 0.5) / scale + tfm[1], allocator);
        for (dist.row(row), 0..) |*d, x| {
            const distance = d[0];
            if ((distance > 0.5) != scanline.filled((f64i(x) + 0.5) / scale + tfm[0], fill_rule))
                d[0] = 1.0 - distance;
        }
    }
}

fn msdfSignCorrection(
    comptime n: usize,
    allocator: Allocator,
    dist: Image([n]f64),
    scale: f64,
    shape: Shape,
    tfm: V2f64,
    fill_rule: Scanline.FillRule,
) !void {
    var scanline: Scanline = .{};
    defer scanline.intersections.deinit(allocator);

    const w = dist.width;
    const h = dist.height;
    const match_map = try allocator.alloc(i32, @as(usize, w) * h);
    defer allocator.free(match_map);
    @memset(match_map, 0);

    var ambiguous = false;
    var match_idx: usize = 0;
    for (0..h) |y| {
        const row: u32 = @intCast(h - y - 1);
        try shape.scanline(&scanline, (f64i(y) + 0.5) / scale + tfm[1], allocator);
        for (dist.row(row), 0..) |*px, x| {
            const filled = scanline.filled((f64i(x) + 0.5) / scale + tfm[0], fill_rule);
            const distance = math.median(px[0], px[1], px[2]);
            if (distance == 0.5) {
                ambiguous = true;
            } else if ((distance > 0.5) != filled) {
                for (px[0..3]) |*v| v.* = 1.0 - v.*;
                match_map[match_idx] = -1;
            } else match_map[match_idx] = 1;
            if (n >= 4 and (px[3] > 0.5) != filled)
                px[3] = 1.0 - px[3];
            match_idx += 1;
        }
    }

    if (!ambiguous) return;

    match_idx = 0;
    for (0..h) |y| {
        const row: u32 = @intCast(h - y - 1);
        for (dist.row(row), 0..) |*px, x| {
            const match = match_map[match_idx];
            if (match == 0) {
                var neighbor_match: i32 = 0;
                if (x > 0) neighbor_match += match_map[match_idx - 1];
                if (x < w - 1) neighbor_match += match_map[match_idx + 1];
                if (y > 0) neighbor_match += match_map[match_idx - w];
                if (y < h - 1) neighbor_match += match_map[match_idx + w];
                if (neighbor_match < 0) {
                    for (px[0..3]) |*v| v.* = 1.0 - v.*;
                }
            }
            match_idx += 1;
        }
    }
}

fn pxRangeNorm(dist: f64, px_range: f64) f64 {
    return (dist + px_range / 2.0) / px_range;
}

pub fn findDistanceAt(
    comptime sdf_type: SdfType,
    shape: Shape,
    p: V2f64,
    px_range: f64,
) switch (sdf_type) {
    .sdf, .psdf => f64,
    inline .msdf, .msdf10, .mtsdf => |ty| [ty.numChannels()]f64,
} {
    const PsdfData = struct {
        dist: EdgeSegment.SignedDist = .init,
        edge: ?*const EdgeSegment = null,
        point_pos: EdgeSegment.PointPosition = .within_segment,
    };

    var true_ch: EdgeSegment.SignedDist = .init;
    var perp_ch: [if (sdf_type == .psdf) 1 else 3]PsdfData = @splat(.{});
    for (shape.contours.items) |contour| for (contour.edges.items) |*edge| {
        const dist, const point_pos = edge.signedDistance(p);

        switch (sdf_type) {
            .sdf, .mtsdf => {
                if (dist.lessThan(true_ch))
                    true_ch = dist;
            },
            .psdf => {
                if (dist.lessThan(perp_ch[0].dist)) perp_ch[0] = .{
                    .dist = dist,
                    .edge = edge,
                    .point_pos = point_pos,
                };
            },
            else => {},
        }

        if (sdf_type != .sdf and sdf_type != .psdf)
            for ([_]struct { channel: EdgeColor, target: *PsdfData }{
                .{ .channel = .red, .target = &perp_ch[0] },
                .{ .channel = .green, .target = &perp_ch[1] },
                .{ .channel = .blue, .target = &perp_ch[2] },
            }) |params|
                if (edge.color.hasChannel(params.channel) and dist.lessThan(params.target.dist)) {
                    params.target.* = .{
                        .dist = dist,
                        .edge = edge,
                        .point_pos = point_pos,
                    };
                };
    };

    if (sdf_type != .sdf) for (&perp_ch) |*psdf| {
        if (psdf.edge) |edge|
            edge.perpDistConvert(&psdf.dist, p, psdf.point_pos);
    };

    return switch (sdf_type) {
        .sdf => pxRangeNorm(true_ch.distance, px_range),
        .psdf => pxRangeNorm(perp_ch[0].dist.distance, px_range),
        .msdf, .msdf10 => .{
            pxRangeNorm(perp_ch[0].dist.distance, px_range),
            pxRangeNorm(perp_ch[1].dist.distance, px_range),
            pxRangeNorm(perp_ch[2].dist.distance, px_range),
        },
        .mtsdf => .{
            pxRangeNorm(perp_ch[0].dist.distance, px_range),
            pxRangeNorm(perp_ch[1].dist.distance, px_range),
            pxRangeNorm(perp_ch[2].dist.distance, px_range),
            pxRangeNorm(true_ch.distance, px_range),
        },
    };
}

fn generate(
    comptime sdf_type: SdfType,
    out: Image(DistPixel(sdf_type)),
    scale: f64,
    shape: Shape,
    px_range: f64,
    tfm: V2f64,
) void {
    for (0..out.height) |y| {
        const row: u32 = @intCast(out.height - y - 1);
        for (out.row(row), 0..) |*px, x| {
            const p = V2f64{
                (f64i(x) + 0.5),
                (f64i(y) + 0.5),
            } / math.splat(scale) + tfm;

            const dist = findDistanceAt(sdf_type, shape, p, px_range);
            px.* = if (comptime sdf_type.numChannels() == 1) .{dist} else dist;
        }
    }
}
