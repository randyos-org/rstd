const std = @import("std");

const f64i = @import("__root__.zig").f64i;
const findDistanceAt = @import("__root__.zig").findDistanceAt;
const math = @import("math.zig");
const Shape = @import("Shape.zig");
const aux_types = @import("../aux_types.zig");
const Image = @import("../../../image.zig").Image;

const V2f64 = aux_types.V2f64;

const f64_nan = std.math.nan(f64);

const artifact_t_epsilon = 0.01;
const protection_radius_tolerance = 1.001;

pub const Mode = enum { indiscriminate, edge_priority, edge_only };
pub const Options = struct {
    mode: Mode = .edge_priority,
    /// Will be forcefully turned off if the scanline pass is enabled
    check_distance: bool = true,
    min_deviation_ratio: f64 = 10.0 / 9.0,
    min_improve_ratio: f64 = 10.0 / 9.0,
};

const StencilFlags = packed struct {
    err: bool = false,
    protected: bool = false,
};

const ClassifierFlags = packed struct {
    candidate: bool = false,
    artifact: bool = false,

    pub fn merge(self: ClassifierFlags, other: ClassifierFlags) ClassifierFlags {
        return @bitCast(@as(u2, @bitCast(self)) | @as(u2, @bitCast(other)));
    }
};

const ColorFlags = packed struct {
    red: bool = false,
    green: bool = false,
    blue: bool = false,
};

const DistanceEvaluationInfo = struct {
    shape: *const Shape,
    px_range: f64 = f64_nan,
    dir: V2f64 = @splat(f64_nan),
    sdf_coord: V2f64 = @splat(f64_nan),
    shape_coord: V2f64 = @splat(f64_nan),
    texel_size: V2f64 = @splat(f64_nan),
    sdf: ?*const [3]f64 = null,
};

/// Flattens MSDF/MTSDF texels whose channel median would draw an artifact.
/// Only the first three channels of each texel are read or corrected.
/// `scanline_pass` disables the distance check, which it would invalidate.
pub fn correct(
    comptime N: usize,
    allocator: std.mem.Allocator,
    shape: *const Shape,
    options: Options,
    scanline_pass: bool,
    scale: f64,
    px_range: f64,
    tfm: V2f64,
    sdf: Image([N]f64),
) !void {
    var opts = options;
    if (scanline_pass and opts.check_distance) {
        std.log.warn("Attempted to use distance-based error correction with a non-null scanline pass, disabling", .{});
        opts.check_distance = false;
    }

    var stencil: Image(StencilFlags) = try .init(allocator, sdf.width, sdf.height, .{});
    defer stencil.destroy(allocator);

    var self: Corrector(N) = .{
        .options = opts,
        .dist_eval = .{ .shape = shape },
        .stencil = stencil,
        .sdf = sdf,
    };

    switch (opts.mode) {
        .edge_priority => {
            self.protectCorners(shape, scale, tfm);
            self.protectEdges(px_range / scale * protection_radius_tolerance);
        },
        .edge_only => {
            for (stencil.data) |*mask| mask.protected = true;
        },
        .indiscriminate => {},
    }

    self.findErrors(scale, px_range, tfm);

    for (stencil.data, sdf.data) |mask, *texel| if (mask.err) {
        const msdf = texel[0..3];
        @memset(msdf, median(msdf));
    };
}

fn Corrector(comptime N: usize) type {
    return struct {
        const Self = @This();

        options: Options,
        dist_eval: DistanceEvaluationInfo,
        stencil: Image(StencilFlags),
        sdf: Image([N]f64),

        fn texel(self: *const Self, x: anytype, y: anytype) *const [3]f64 {
            return self.sdf.at(@intCast(x), @intCast(y))[0..3];
        }

        fn mask(self: *const Self, x: anytype, y: anytype) *StencilFlags {
            return self.stencil.at(@intCast(x), @intCast(y));
        }

        fn protectCorners(self: *Self, shape: *const Shape, scale: f64, tfm: V2f64) void {
            const w = self.stencil.width;
            const h = self.stencil.height;
            for (shape.contours.items) |contour| {
                if (contour.edges.items.len == 0) continue;

                var prev_edge = contour.edges.last().?;
                for (contour.edges.items) |edge| {
                    const common_color = @backingInt(prev_edge.color) & @backingInt(edge.color);
                    prev_edge = edge;
                    if ((common_color & (common_color - 1)) == 0) continue;
                    const base_point = edge.point(0) * math.splatV2f64(scale) + tfm;
                    const left: i32 = @trunc(base_point[0] - 0.5);
                    const bottom: i32 = @trunc(f64i(h) - @floor(base_point[1] - 0.5) - 2.0);
                    const right = left + 1;
                    const top = bottom + 1;
                    if (left < w and bottom < h and right >= 0 and top >= 0) {
                        if (left >= 0 and bottom >= 0) self.mask(left, bottom).protected = true;
                        if (right < w and bottom >= 0) self.mask(right, bottom).protected = true;
                        if (left >= 0 and top < h) self.mask(left, top).protected = true;
                        if (right < w and top < h) self.mask(right, top).protected = true;
                    }
                }
            }
        }

        fn protectEdges(self: *Self, radius: f64) void {
            const sdf_w = self.sdf.width;
            const sdf_h = self.sdf.height;
            for (0..sdf_h) |y| {
                const left = self.texel(0, y);
                const right = self.texel(1, y);
                const median_left = median(left);
                const median_right = median(right);
                for (0..sdf_w - 1) |x| if (@abs(median_left - 0.5) + @abs(median_right - 0.5) < radius) {
                    const m = edgeBetweenTexels(left, right);
                    protectExtremeChannels(self.mask(x, y), left, median_left, m);
                    protectExtremeChannels(self.mask(x + 1, y), right, median_right, m);
                };
            }

            for (0..sdf_h - 1) |y| {
                const bottom = self.texel(0, y);
                const top = self.texel(0, y + 1);
                const median_bottom = median(bottom);
                const median_top = median(top);
                for (0..sdf_w) |x| if (@abs(median_bottom - 0.5) + @abs(median_top - 0.5) < radius) {
                    const m = edgeBetweenTexels(bottom, top);
                    protectExtremeChannels(self.mask(x, y), bottom, median_bottom, m);
                    protectExtremeChannels(self.mask(x, y + 1), top, median_top, m);
                };
            }

            const diag_radius = radius * @sqrt(2.0);
            for (0..sdf_h - 1) |y| {
                const bottom_left = self.texel(0, y);
                const bottom_right = self.texel(1, y);
                const top_left = self.texel(0, y + 1);
                const top_right = self.texel(1, y + 1);
                const median_bottom_left = median(bottom_left);
                const median_bottom_right = median(bottom_right);
                const median_top_left = median(top_left);
                const median_top_right = median(top_right);
                for (0..sdf_w - 1) |x| {
                    if (@abs(median_bottom_left - 0.5) + @abs(median_top_right - 0.5) < diag_radius) {
                        const m = edgeBetweenTexels(bottom_left, top_right);
                        protectExtremeChannels(self.mask(x, y), bottom_left, median_bottom_left, m);
                        protectExtremeChannels(self.mask(x + 1, y + 1), top_right, median_top_right, m);
                    }
                    if (@abs(median_bottom_right - 0.5) + @abs(median_top_left - 0.5) < diag_radius) {
                        const m = edgeBetweenTexels(bottom_right, top_left);
                        protectExtremeChannels(self.mask(x + 1, y), bottom_right, median_bottom_right, m);
                        protectExtremeChannels(self.mask(x, y + 1), top_left, median_top_left, m);
                    }
                }
            }
        }

        /// mist's clamps here swap w and h, and `lr_px` reads (left, right):
        /// either can step past a row's end, so keep flat indexing until fixed.
        fn texelFlat(self: *const Self, x: u32, y: u32) *const [3]f64 {
            return self.sdf.data[@as(usize, y) * self.sdf.width + x][0..3];
        }

        fn evaluateArtifact(self: *const Self, flags: ClassifierFlags, t: f64) bool {
            if (flags.artifact) return true;
            if (!self.options.check_distance or !flags.candidate) return false;

            const t_vec = math.splatV2f64(t) * self.dist_eval.dir;

            const tx_sdf_coord = self.dist_eval.sdf_coord + t_vec;
            const tx_x = tx_sdf_coord[0] - 0.5;
            const tx_y = tx_sdf_coord[1] - 0.5;
            const floor_x = @floor(tx_x);
            const floor_y = @floor(tx_y);
            const lr = tx_x - floor_x;
            const bt = tx_y - floor_y;
            const sdf_w = self.sdf.width;
            const sdf_h = self.sdf.height;
            const left = @min(@as(u32, @trunc(floor_x)), sdf_w - 1);
            const bottom = @min(@as(u32, @trunc(floor_y)), sdf_w - 1);
            const right = @min(left + 1, sdf_h - 1);
            const top = @min(bottom + 1, sdf_h - 1);
            const lb_px = self.texelFlat(left, bottom);
            const rb_px = self.texelFlat(right, bottom);
            const lt_px = self.texelFlat(left, top);
            const lr_px = self.texelFlat(left, right);
            var old_sdf: [3]f64 = undefined;
            for (&old_sdf, 0..) |*c, i|
                c.* = math.mix(
                    math.mix(lb_px[i], rb_px[i], lr),
                    math.mix(lt_px[i], lr_px[i], lr),
                    bt,
                );

            const wt = (1 - @abs(t_vec[0])) * (1 - @abs(t_vec[1]));
            const sdf = self.dist_eval.sdf.?;
            const m = median(sdf);
            const new_sdf: [3]f64 = .{
                old_sdf[0] + wt * (m - sdf[0]),
                old_sdf[1] + wt * (m - sdf[1]),
                old_sdf[2] + wt * (m - sdf[2]),
            };
            const om = median(&old_sdf);
            const nm = median(&new_sdf);
            const dist = findDistanceAt(
                .psdf,
                self.dist_eval.shape.*,
                self.dist_eval.shape_coord + t_vec * self.dist_eval.texel_size,
                self.dist_eval.px_range,
            );
            return self.options.min_improve_ratio * @abs(nm - dist) < @abs(om - dist);
        }

        fn hasLinearArtifact(
            self: *const Self,
            span: f64,
            protected: bool,
            am: f64,
            a: *const [3]f64,
            b: *const [3]f64,
        ) bool {
            const bm = median(b);
            if (@abs(am - 0.5) < @abs(bm - 0.5)) return false;
            for (0..3) |idx| {
                const next_idx = (idx + 1) % 3;
                const da = a[next_idx] - a[idx];
                const db = b[next_idx] - b[idx];
                const delta = da - db;
                if (delta == 0) continue;
                const t = da / (da - db);
                if (t > artifact_t_epsilon and t < 1 - artifact_t_epsilon) {
                    const xm = math.median(
                        math.mix(a[0], b[0], t),
                        math.mix(a[1], b[1], t),
                        math.mix(a[2], b[2], t),
                    );
                    if (self.evaluateArtifact(rangeTest(span, protected, 0, 1, t, am, bm, xm), t))
                        return true;
                }
            }
            return false;
        }

        fn hasDiagonalArtifact(
            self: *const Self,
            span: f64,
            protected: bool,
            am: f64,
            a: *const [3]f64,
            b: *const [3]f64,
            c: *const [3]f64,
            d: *const [3]f64,
        ) bool {
            const dm = median(d);
            if (@abs(am - 0.5) < @abs(dm - 0.5)) return false;

            const abc: [3]f64 = .{
                a[0] - b[0] - c[0],
                a[1] - b[1] - c[1],
                a[2] - b[2] - c[2],
            };
            const q: [3]f64 = .{
                d[0] + abc[0],
                d[1] + abc[1],
                d[2] + abc[2],
            };
            const l: [3]f64 = .{
                -a[0] - abc[0],
                -a[1] - abc[1],
                -a[2] - abc[2],
            };
            const t_ex: [3]f64 = .{
                -0.5 * l[0] / q[0],
                -0.5 * l[1] / q[1],
                -0.5 * l[2] / q[2],
            };

            for (0..3) |idx| {
                const next_idx = (idx + 1) % 3;
                const d_a = a[next_idx] - a[idx];
                const d_bc = b[next_idx] - b[idx] + (c[next_idx] - c[idx]);
                const d_d = d[next_idx] - d[idx];
                const t_ex_0 = t_ex[idx];
                const t_ex_1 = t_ex[next_idx];

                var buf: [2]f64 = undefined;
                for (math.solveEquation(&buf, &.{
                    .a = d_d - d_bc + d_a,
                    .b = d_bc - d_a - d_a,
                    .c = d_a,
                })) |root|
                    if (root > artifact_t_epsilon and root < 1 - artifact_t_epsilon) {
                        const xm = interpolatedMedianBilinear(a, &l, &q, root);
                        var flags = rangeTest(span, protected, 0, 1, root, am, dm, xm);
                        var t_end: [2]f64 = undefined;
                        var em: [2]f64 = undefined;

                        if (t_ex_0 > 0 and t_ex_0 < 1) {
                            t_end[0] = 0;
                            t_end[1] = 1;
                            em[0] = am;
                            em[1] = dm;
                            t_end[@intFromBool(t_ex_0 > root)] = t_ex_0;
                            em[@intFromBool(t_ex_0 > root)] = interpolatedMedianBilinear(a, &l, &q, t_ex_0);
                            flags = flags.merge(rangeTest(span, protected, t_end[0], t_end[1], root, em[0], em[1], xm));
                        }

                        if (t_ex_1 > 0 and t_ex_1 < 1) {
                            t_end[0] = 0;
                            t_end[1] = 1;
                            em[0] = am;
                            em[1] = dm;
                            t_end[@intFromBool(t_ex_1 > root)] = t_ex_1;
                            em[@intFromBool(t_ex_1 > root)] = interpolatedMedianBilinear(a, &l, &q, t_ex_1);
                            flags = flags.merge(rangeTest(span, protected, t_end[0], t_end[1], root, em[0], em[1], xm));
                        }

                        if (self.evaluateArtifact(flags, root))
                            return true;
                    };
            }

            return false;
        }

        fn findErrors(self: *Self, scale: f64, px_range: f64, tfm: V2f64) void {
            const min_deviation_ratio = self.options.min_deviation_ratio;
            const vscale = math.splatV2f64(scale);
            const span = px_range / scale * min_deviation_ratio;
            const diag_span = span * @sqrt(2.0);
            const sdf_w = self.sdf.width;
            const sdf_h = self.sdf.height;

            self.dist_eval.px_range = px_range;
            self.dist_eval.texel_size = math.splatV2f64(1.0) / vscale;

            for (0..sdf_h) |y| for (0..sdf_w) |x| {
                const current = self.texel(x, y);
                const median_current = median(current);
                const current_stencil = self.mask(x, y);
                const is_protected = current_stencil.protected;

                const fx: f64 = @floatFromInt(x);
                const fy: f64 = @floatFromInt(y);
                self.dist_eval.sdf = current;
                self.dist_eval.sdf_coord = .{ fx + 0.5, fy + 0.5 };
                self.dist_eval.shape_coord = self.dist_eval.sdf_coord / vscale + tfm;

                if (x > 0) {
                    const left = self.texel(x - 1, y);
                    self.dist_eval.dir = .{ -1, 0 };
                    if (self.hasLinearArtifact(span, is_protected, median_current, current, left)) {
                        current_stencil.err = true;
                        continue;
                    }

                    if (y > 0) {
                        const top = self.texel(x, y - 1);
                        const top_left = self.texel(x - 1, y - 1);
                        self.dist_eval.dir = .{ -1, -1 };
                        if (self.hasDiagonalArtifact(diag_span, is_protected, median_current, current, left, top, top_left)) {
                            current_stencil.err = true;
                            continue;
                        }
                    }

                    if (y < sdf_h - 1) {
                        const bottom = self.texel(x, y + 1);
                        const bottom_left = self.texel(x - 1, y + 1);
                        self.dist_eval.dir = .{ -1, 1 };
                        if (self.hasDiagonalArtifact(diag_span, is_protected, median_current, current, left, bottom, bottom_left)) {
                            current_stencil.err = true;
                            continue;
                        }
                    }
                }

                if (y > 0) {
                    const top = self.texel(x, y - 1);
                    self.dist_eval.dir = .{ 0, -1 };
                    if (self.hasLinearArtifact(span, is_protected, median_current, current, top)) {
                        current_stencil.err = true;
                        continue;
                    }
                }

                if (x < sdf_w - 1) {
                    const right = self.texel(x + 1, y);
                    self.dist_eval.dir = .{ 1, 0 };
                    if (self.hasLinearArtifact(span, is_protected, median_current, current, right)) {
                        current_stencil.err = true;
                        continue;
                    }

                    if (y > 0) {
                        const top = self.texel(x, y - 1);
                        const top_right = self.texel(x + 1, y - 1);
                        self.dist_eval.dir = .{ 1, -1 };
                        if (self.hasDiagonalArtifact(diag_span, is_protected, median_current, current, right, top, top_right)) {
                            current_stencil.err = true;
                            continue;
                        }
                    }

                    if (y < sdf_h - 1) {
                        const bottom = self.texel(x, y + 1);
                        const bottom_right = self.texel(x + 1, y + 1);
                        self.dist_eval.dir = .{ 1, 1 };
                        if (self.hasDiagonalArtifact(diag_span, is_protected, median_current, current, right, bottom, bottom_right)) {
                            current_stencil.err = true;
                            continue;
                        }
                    }
                }

                if (y < sdf_h - 1) {
                    const bottom = self.texel(x, y + 1);
                    self.dist_eval.dir = .{ 0, 1 };
                    if (self.hasLinearArtifact(span, is_protected, median_current, current, bottom)) {
                        current_stencil.err = true;
                        continue;
                    }
                }
            };
        }
    };
}

fn edgeBetweenTexels(a: *const [3]f64, b: *const [3]f64) ColorFlags {
    var mask: ColorFlags = .{};
    for (0..3) |channel| {
        const delta = a[channel] - b[channel];
        if (delta == 0.0)
            continue;

        const t = (a[channel] - 0.5) / delta;
        if (t < 0.0 or t > 1.0)
            continue;

        const c: [3]f64 = .{
            math.mix(a[0], b[0], t),
            math.mix(a[1], b[1], t),
            math.mix(a[2], b[2], t),
        };
        if (median(&c) == c[channel])
            switch (channel) {
                0 => mask.red = true,
                1 => mask.green = true,
                2 => mask.blue = true,
                else => unreachable,
            };
    }

    return mask;
}

fn protectExtremeChannels(stencil_point: *StencilFlags, msd: *const [3]f64, m: f64, mask: ColorFlags) void {
    if ((mask.red and msd[0] != m) or
        (mask.green and msd[1] != m) or
        (mask.blue and msd[2] != m))
        stencil_point.protected = true;
}

fn median(arr: *const [3]f64) f64 {
    return math.median(arr[0], arr[1], arr[2]);
}

fn interpolatedMedianBilinear(a: *const [3]f64, l: *const [3]f64, q: *const [3]f64, t: f64) f64 {
    return math.median(
        t * (t * q[0] + l[0]) + a[0],
        t * (t * q[1] + l[1]) + a[1],
        t * (t * q[2] + l[2]) + a[2],
    );
}

fn rangeTest(span: f64, protected: bool, at: f64, bt: f64, xt: f64, am: f64, bm: f64, xm: f64) ClassifierFlags {
    if (!((am > 0.5 and bm > 0.5 and xm <= 0.5) or
        (am < 0.5 and bm < 0.5 and xm >= 0.5) or
        (!protected and math.median(am, bm, xm) != xm)))
        return .{};

    const ax_span = (xt - at) * span;
    const bx_span = (bt - xt) * span;
    return .{
        .candidate = true,
        .artifact = !(xm >= am - ax_span and xm <= am + ax_span and xm >= bm - bx_span and xm <= bm + bx_span),
    };
}
