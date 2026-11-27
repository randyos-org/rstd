//! The outline the CFF engine draws into: FreeType 2.14.3
//! `src/psaux/psobjs.c` (`ps_builder_add_point`, `ps_builder_add_contour`,
//! `ps_builder_start_point`, `ps_builder_close_contour`) and the
//! callbacks of `psft.c` (`cf2_builder_moveTo`, `lineTo`, `cubeTo`).
//! Points arrive in 16.16 and are kept in 26.6.

const std = @import("std");

const Allocator = std.mem.Allocator;
const Builder = @This();

pub const tag_on: u8 = 1;
pub const tag_cubic: u8 = 2;

allocator: Allocator,
points: std.ArrayList([2]i32) = .empty,
tags: std.ArrayList(u8) = .empty,
contours: std.ArrayList(u16) = .empty,
path_begun: bool = false,

pub fn deinit(b: *Builder) void {
    b.points.deinit(b.allocator);
    b.tags.deinit(b.allocator);
    b.contours.deinit(b.allocator);
}

fn addPoint(b: *Builder, p: [2]i32, tag: u8) Allocator.Error!void {
    try b.points.append(b.allocator, .{ p[0] >> 10, p[1] >> 10 });
    try b.tags.append(b.allocator, tag);
}

fn addContour(b: *Builder) Allocator.Error!void {
    if (b.contours.items.len > 0) b.contours.items[b.contours.items.len - 1] = @intCast(b.points.items.len -% 1);
    try b.contours.append(b.allocator, 0);
}

/// Opens a contour at `p` unless one is open.
fn startPoint(b: *Builder, p: [2]i32) Allocator.Error!void {
    if (b.path_begun) return;
    b.path_begun = true;
    try b.addContour();
    try b.addPoint(p, tag_on);
}

/// Ends the open contour; a last point on the first one is dropped, and a
/// contour of one point or none with it.
pub fn closeContour(b: *Builder) void {
    const nc = b.contours.items.len;
    const n = b.points.items.len;
    const first: usize = if (nc <= 1) 0 else @as(usize, b.contours.items[nc - 2]) + 1;
    if (nc > 0 and first == n) {
        b.contours.items.len -= 1;
        return;
    }
    if (n > 1) {
        const p1 = b.points.items[first];
        const p2 = b.points.items[n - 1];
        if (p1[0] == p2[0] and p1[1] == p2[1] and b.tags.items[n - 1] == tag_on) b.dropPoint();
    }
    if (nc > 0) {
        if (first == b.points.items.len - 1) {
            b.contours.items.len -= 1;
            b.dropPoint();
        } else {
            b.contours.items[nc - 1] = @intCast(b.points.items.len - 1);
        }
    }
}

fn dropPoint(b: *Builder) void {
    b.points.items.len -= 1;
    b.tags.items.len -= 1;
}

pub fn moveTo(b: *Builder) void {
    b.closeContour();
    b.path_begun = false;
}

pub fn lineTo(b: *Builder, p0: [2]i32, p1: [2]i32) Allocator.Error!void {
    try b.startPoint(p0);
    try b.addPoint(p1, tag_on);
}

pub fn cubeTo(b: *Builder, p0: [2]i32, p1: [2]i32, p2: [2]i32, p3: [2]i32) Allocator.Error!void {
    try b.startPoint(p0);
    try b.addPoint(p1, tag_cubic);
    try b.addPoint(p2, tag_cubic);
    try b.addPoint(p3, tag_on);
}
