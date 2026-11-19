//! The signed distance from a point to an edge segment: linear, quadratic and cubic.

const EdgeSegment = @import("EdgeSegment.zig");
const math = @import("math.zig");
const aux_types = @import("../sfnt/aux_types.zig");
const dot = math.dot;
const cross = math.cross;
const normal = math.normal;
const splat = math.splat;
const V2f64 = aux_types.V2f64;
const SignedDist = EdgeSegment.SignedDist;
const PointPosition = EdgeSegment.PointPosition;

const cubic_starts = 4;
const cubic_steps = 4;

pub fn signedDistance(self: EdgeSegment, origin: V2f64) struct { SignedDist, PointPosition } {
    switch (self.segment) {
        .linear => |p| {
            const aq = origin - p[0];
            const ab = p[1] - p[0];
            const t = dot(aq, ab) / dot(ab, ab);
            const eq = p[@intFromBool(t > 0.5)] - origin;
            const point_pos: PointPosition = .fromNormalized(t);
            const endpoint_dist = math.length(eq);
            if (point_pos == .within_segment) {
                const ortho_dist = dot(-math.orthonormal(ab, true), aq);
                if (@abs(ortho_dist) < endpoint_dist)
                    return .{ .{ .distance = ortho_dist }, point_pos };
            }

            return .{
                .{
                    .distance = math.nonZeroSign(cross(aq, ab)) * endpoint_dist,
                    .dot = @abs(dot(normal(ab, true), normal(eq, true))),
                },
                point_pos,
            };
        },
        .quadratic_bezier => |p| {
            const qa = p[0] - origin;
            const qc = p[2] - origin;
            const ab = p[1] - p[0];
            const br = p[2] - p[1] - ab;
            const d0 = self.direction(0);
            const d1 = self.direction(1);

            var min_dist = math.nonZeroSign(cross(d0, qa)) * math.length(qa);
            var point_pos: PointPosition = .fromNormalized(-dot(qa, d0) / dot(d0, d0));

            var dist = math.length(qc);
            if (dist < @abs(min_dist)) {
                min_dist = math.nonZeroSign(cross(d1, qc)) * dist;
                point_pos = .fromNormalized(dot(origin - p[1], d1) / dot(d1, d1));
            }

            var buf: [3]f64 = undefined;
            for (math.solveEquation(&buf, &.{
                .a = dot(br, br),
                .b = 3.0 * dot(ab, br),
                .c = 2.0 * dot(ab, ab) + dot(qa, br),
                .d = dot(qa, ab),
            })) |root|
                if (root > 0.0 and root < 1.0) {
                    const qe = qa + ab * splat(root * 2.0) + br * splat(root * root);
                    dist = math.length(qe);
                    if (dist < @abs(min_dist)) {
                        min_dist = math.nonZeroSign(cross(ab + br * splat(root), qe)) * dist;
                        point_pos = .fromNormalized(root);
                    }
                };

            return .{
                .{
                    .distance = min_dist,
                    .dot = switch (point_pos) {
                        .within_segment => 0.0,
                        .outside_start => @abs(dot(
                            normal(d0, true),
                            normal(qa, true),
                        )),
                        .outside_end => @abs(dot(
                            normal(d1, true),
                            normal(qc, true),
                        )),
                    },
                },
                point_pos,
            };
        },
        .cubic_bezier => |p| {
            const qa = p[0] - origin;
            const qd = p[3] - origin;
            const ab = p[1] - p[0];
            const br = p[2] - p[1] - ab;
            const as = (p[3] - p[2]) - (p[2] - p[1]) - br;
            const d0 = self.direction(0);
            const d1 = self.direction(1);

            var min_dist = math.nonZeroSign(cross(d0, qa)) * math.length(qa);
            var point_pos: PointPosition = .fromNormalized(-dot(qa, d0) / dot(d0, d0));

            const qd_dist = math.length(qd);
            if (qd_dist < @abs(min_dist)) {
                min_dist = math.nonZeroSign(cross(d1, qd)) * qd_dist;
                point_pos = .fromNormalized(dot(d1 - qd, d1) / dot(d1, d1));
            }

            // cubic_starts + 1 seeds, t = 0..1 inclusive, as msdfgen
            for (0..cubic_starts + 1) |i| {
                const fi: f64 = @floatFromInt(i);
                var t = fi / cubic_starts;

                var qe = qa +
                    ab * splat(t * 3.0) +
                    br * splat(t * t * 3.0) +
                    as * splat(t * t * t);
                var cd1 = ab * splat(3.0) +
                    br * splat(t * 6.0) +
                    as * splat(t * t * 3.0);
                var cd2 = br * splat(6.0) +
                    as * splat(t * 6.0);

                var next_t = t - dot(qe, cd1) / (dot(cd1, cd1) + dot(qe, cd2));
                if (next_t > 0.0 and next_t < 1.0) {
                    for (0..cubic_steps) |_| {
                        t = next_t;

                        qe = qa +
                            ab * splat(t * 3.0) +
                            br * splat(t * t * 3.0) +
                            as * splat(t * t * t);
                        cd1 = ab * splat(3.0) +
                            br * splat(t * 6.0) +
                            as * splat(t * t * 3.0);
                        cd2 = br * splat(6.0) +
                            as * splat(t * 6.0);

                        next_t = t - dot(qe, cd1) / (dot(cd1, cd1) + dot(qe, cd2));
                        if (next_t <= 0.0 or next_t >= 1.0) break;
                    }

                    const qe_dist = math.length(qe);
                    if (qe_dist < @abs(min_dist)) {
                        min_dist = math.nonZeroSign(cross(cd1, qe)) * qe_dist;
                        point_pos = .fromNormalized(t);
                    }
                }
            }

            return .{
                .{
                    .distance = min_dist,
                    .dot = switch (point_pos) {
                        .within_segment => 0.0,
                        .outside_start => @abs(dot(
                            normal(d0, true),
                            normal(qa, true),
                        )),
                        .outside_end => @abs(dot(
                            normal(d1, true),
                            normal(qd, true),
                        )),
                    },
                },
                point_pos,
            };
        },
    }
}
