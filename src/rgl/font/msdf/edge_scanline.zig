//! Where a horizontal scanline crosses an edge segment.

const std = @import("std");
const EdgeSegment = @import("EdgeSegment.zig");
const math = @import("math.zig");
const mix = math.mix;

pub fn scanlineIntersections(self: EdgeSegment, x: *[3]f64, dy: *[3]i32, y: f64) u32 {
    switch (self.segment) {
        .linear => |p| {
            if ((y >= p[0][1] and y < p[1][1]) or (y >= p[1][1] and y < p[0][1])) {
                const param = (y - p[0][1]) / (p[1][1] - p[0][1]);
                x[0] = mix(p[0][0], p[1][0], param);
                dy[0] = std.math.sign(p[1][1] - p[0][1]);
                return 1;
            }
            return 0;
        },
        .quadratic_bezier => |p| {
            var total: u32 = 0;
            var next_dy: i32 = math.boolSign(y > p[0][1]);
            x[total] = p[0][0];
            if (p[0][1] == y) {
                if (p[0][1] < p[1][1] or (p[0][1] == p[1][1] and p[0][1] < p[2][1])) {
                    dy[total] = 1;
                    total += 1;
                } else next_dy = 1;
            }

            const ab = p[1] - p[0];
            const br = p[2] - p[1] - ab;

            var buf: [2]f64 = undefined;
            const roots = math.solveEquation(&buf, &.{
                .a = br[1],
                .b = 2 * ab[1],
                .c = p[0][1] - y,
            });
            std.sort.pdq(f64, roots, {}, struct {
                pub fn lessThan(_: void, a: f64, b: f64) bool {
                    return a < b;
                }
            }.lessThan);

            for (roots) |root|
                if (total < 2 and root >= 0 and root <= 1) {
                    x[total] = p[0][0] + 2 * root * ab[0] + root * root * br[0];
                    if (@as(f64, @floatFromInt(next_dy)) * (ab[1] + root * br[1]) >= 0) {
                        dy[total] = next_dy;
                        total += 1;
                        next_dy = -next_dy;
                    }
                };

            if (p[2][1] == y) {
                if (next_dy > 0 and total > 0) {
                    total -= 1;
                    next_dy = -1;
                }
                if ((p[2][1] < p[1][1] or (p[2][1] == p[1][1] and p[2][1] < p[0][1])) and total < 2) {
                    x[total] = p[2][0];
                    if (next_dy < 0) {
                        dy[total] = -1;
                        total += 1;
                        next_dy = 1;
                    }
                }
            }

            if (next_dy != math.boolSign(y >= p[2][1])) {
                if (total > 0)
                    total -= 1
                else {
                    if (@abs(p[2][1] - y) < @abs(p[0][1] - y)) x[total] = p[2][0];
                    dy[total] = next_dy;
                    total += 1;
                }
            }

            return total;
        },
        .cubic_bezier => |p| {
            var total: u32 = 0;
            var next_dy: i32 = math.boolSign(y > p[0][1]);
            x[total] = p[0][0];
            if (p[0][1] == y) {
                if (p[0][1] < p[1][1] or (p[0][1] == p[1][1] and (p[0][1] < p[2][1] or (p[0][1] == p[2][1] and p[0][1] < p[3][1])))) {
                    dy[total] = 1;
                    total += 1;
                } else next_dy = 1;
            }

            const ab = p[1] - p[0];
            const br = p[2] - p[1] - ab;
            const as = (p[3] - p[2]) - (p[2] - p[1]) - br;

            var buf: [3]f64 = undefined;
            const roots = math.solveEquation(&buf, &.{
                .a = as[1],
                .b = 3 * br[1],
                .c = 3 * ab[1],
                .d = p[0][1] - y,
            });
            std.sort.pdq(f64, roots, {}, struct {
                pub fn lessThan(_: void, a: f64, b: f64) bool {
                    return a < b;
                }
            }.lessThan);

            for (roots) |root|
                if (total < 3 and root >= 0 and root <= 1) {
                    x[total] = p[0][0] +
                        3 * root * ab[0] +
                        3 * root * root * br[0] +
                        root * root * root * as[0];
                    if (@as(f64, @floatFromInt(next_dy)) * (ab[1] + 2 * root * br[1] + root * root * as[1]) >= 0) {
                        dy[total] = next_dy;
                        total += 1;
                        next_dy = -next_dy;
                    }
                };

            if (p[3][1] == y) {
                if (next_dy > 0 and total > 0) {
                    total -= 1;
                    next_dy = -1;
                }

                if ((p[3][1] < p[2][1] or (p[3][1] == p[2][1] and (p[3][1] < p[1][1] or (p[3][1] == p[1][1] and p[3][1] < p[0][1])))) and total < 3) {
                    x[total] = p[3][0];
                    if (next_dy < 0) {
                        dy[total] = -1;
                        total += 1;
                        next_dy = 1;
                    }
                }
            }

            if (next_dy != math.boolSign(y >= p[3][1])) {
                if (total > 0)
                    total -= 1
                else {
                    if (@abs(p[3][1] - y) < @abs(p[0][1] - y)) x[total] = p[3][0];
                    dy[total] = next_dy;
                    total += 1;
                }
            }

            return total;
        },
    }
}
