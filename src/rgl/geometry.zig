//! Shared, backend-agnostic geometry types; used by both the sw and vulkan
//! rgl backends, and by the public render.zig API.

const zigstd = @import("../zigstd.zig");

pub const Point = struct {
    const Self = @This();
    x: u32,
    y: u32,

    pub const zero: Self = .{ .x = 0, .y = 0 };

    pub fn toPointF32(self: Self) PointF32 {
        return .{ .x = @floatFromInt(self.x), .y = @floatFromInt(self.y) };
    }

    pub fn toPointI32(self: Self) PointI32 {
        return .{ .x = @intCast(self.x), .y = @intCast(self.y) };
    }
};

pub const PointI32 = struct {
    const Self = @This();
    x: i32,
    y: i32,

    pub const zero: Self = .{ .x = 0, .y = 0 };

    pub fn toPointF32(self: Self) PointF32 {
        return .{ .x = @floatFromInt(self.x), .y = @floatFromInt(self.y) };
    }

    pub fn toPoint(self: Self) Point {
        return .{ .x = @intCast(self.x), .y = @intCast(self.y) };
    }

    pub fn add(self: Self, other: Self) Self {
        return .{ .x = self.x + other.x, .y = self.y + other.y };
    }
};

pub const PointF32 = struct {
    const Self = @This();
    x: f32,
    y: f32,

    pub const zero: Self = .{ .x = 0, .y = 0 };

    pub fn toPoint(self: Self) Point {
        return .{ .x = @round(self.x), .y = @round(self.y) };
    }

    pub fn toPointI32(self: Self) PointI32 {
        return .{ .x = @round(self.x), .y = @round(self.y) };
    }

    pub fn scaled(self: Self, s: f32) Self {
        return .{ .x = self.x * s, .y = self.y * s };
    }

    pub fn add(self: Self, other: Self) Self {
        return .{ .x = self.x + other.x, .y = self.y + other.y };
    }

    pub fn sub(self: Self, other: Self) Self {
        return .{ .x = self.x - other.x, .y = self.y - other.y };
    }
};

/// Axis-aligned, origin top-left, integer pixel rectangle; the unit every
/// render.zig box primitive (fillRect/strokeRect/fillRoundRect/...) takes.
pub const Rect = struct {
    const Self = @This();
    x: i32,
    y: i32,
    w: u32,
    h: u32,

    pub fn topLeft(self: Self) PointI32 {
        return .{ .x = self.x, .y = self.y };
    }

    pub fn topLeftF32(self: Self) PointF32 {
        return self.topLeft().toPointF32();
    }

    pub fn bottomRight(self: Self) PointI32 {
        return .{ .x = self.x + @as(i32, @intCast(self.w)), .y = self.y + @as(i32, @intCast(self.h)) };
    }

    pub fn bottomRightF32(self: Self) PointF32 {
        return self.bottomRight().toPointF32();
    }

    pub fn center(self: Self) PointF32 {
        return .{
            .x = @as(f32, @floatFromInt(self.x)) + @as(f32, @floatFromInt(self.w)) / 2,
            .y = @as(f32, @floatFromInt(self.y)) + @as(f32, @floatFromInt(self.h)) / 2,
        };
    }
};

const testExpectEqual = zigstd.testExpectEqual;

test "Point/PointF32/PointI32 conversions and arithmetic" {
    const p: Point = .{ .x = 3, .y = 5 };
    try testExpectEqual(PointF32{ .x = 3, .y = 5 }, p.toPointF32());
    try testExpectEqual(PointI32{ .x = 3, .y = 5 }, p.toPointI32());

    const pf: PointF32 = .{ .x = 2.6, .y = 4.4 };
    try testExpectEqual(Point{ .x = 3, .y = 4 }, pf.toPoint());
    try testExpectEqual(PointI32{ .x = 3, .y = 4 }, pf.toPointI32());

    const pi: PointI32 = .{ .x = 2, .y = 7 };
    try testExpectEqual(PointF32{ .x = 2, .y = 7 }, pi.toPointF32());
    try testExpectEqual(Point{ .x = 2, .y = 7 }, pi.toPoint());
    try testExpectEqual(PointI32{ .x = 7, .y = 12 }, pi.add(.{ .x = 5, .y = 5 }));

    try testExpectEqual(PointF32{ .x = 5, .y = 10 }, (PointF32{ .x = 1, .y = 2 }).scaled(5));
    try testExpectEqual(PointF32{ .x = 4, .y = 6 }, (PointF32{ .x = 1, .y = 2 }).add(.{ .x = 3, .y = 4 }));
    try testExpectEqual(PointF32{ .x = -2, .y = -2 }, (PointF32{ .x = 1, .y = 2 }).sub(.{ .x = 3, .y = 4 }));
}

test "Rect corners and center" {
    const r: Rect = .{ .x = 10, .y = 20, .w = 30, .h = 40 };
    try testExpectEqual(PointI32{ .x = 10, .y = 20 }, r.topLeft());
    try testExpectEqual(PointI32{ .x = 40, .y = 60 }, r.bottomRight());
    try testExpectEqual(PointF32{ .x = 25, .y = 40 }, r.center());
}
