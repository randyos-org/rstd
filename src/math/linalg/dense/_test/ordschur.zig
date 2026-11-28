//! Tests: `ordschur`.

const std = @import("std");
const matrix = @import("../../matrix.zig");
const schur_form = @import("../schur_form.zig");
const ordschur_ = @import("../ordschur.zig");
const similarity = @import("similarity.zig");

const testing = std.testing;
const ConstView = matrix.ConstView;
const ordschur = ordschur_.ordschur;
const Region = ordschur_.Region;
const residual = similarity.residual;
const eye4 = [_]f64{ 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };

test "ordschur: the selected eigenvalues lead T, U T Uᴴ unchanged" {
    const al = testing.allocator;
    const io = testing.io;
    // triangular with U = I: 3 and 4 move to the front, in their order
    const a = [_]f64{ 1, 0, 0, 0, 0.5, 2, 0, 0, -1, 0.25, 3, 0, 2, 1, -0.5, 4 };
    const v: ConstView(f64) = try .init(&a, 4, 4, 4);
    const id: ConstView(f64) = try .init(&eye4, 4, 4, 4);
    const s = try ordschur(al, io, id, v, .{ .mask = &.{ false, false, true, true } });
    defer s.deinit(al);
    try testing.expect(residual(&a, s.u, s.t) < 1e-13 * 8);
    try testing.expectApproxEqAbs(3, s.t.at(0, 0), 1e-13);
    try testing.expectApproxEqAbs(4, s.t.at(1, 1), 1e-13);
    for (0..4) |j| for (j + 1..4) |i| try testing.expectEqual(0, s.t.at(i, j));

    // a real pair ahead of 2: selecting 2 moves it past the 2 x 2 block
    const r = [_]f64{ 0, -1, 0, 1, 0, 0, 0, 0, 2 };
    const rv = try ConstView(f64).init(&r, 3, 3, 3);
    const sr = try schur_form.schur(al, io, rv);
    defer sr.deinit(al);
    var at2: usize = 0;
    for (0..3) |i| {
        if (@abs(sr.t.at(i, i) - 2) < 1e-12) at2 = i;
    }
    var pick = [_]bool{ false, false, false };
    pick[at2] = true;
    const o = try ordschur(al, io, sr.u, sr.t, .{ .mask = &pick });
    defer o.deinit(al);
    try testing.expect(residual(&r, o.u, o.t) < 1e-13 * 8);
    try testing.expectApproxEqAbs(2, o.t.at(0, 0), 1e-13);

    // complex: the same with schurComplex's triangular form
    const sc = try schur_form.schurComplex(al, io, rv);
    defer sc.deinit(al);
    var cpick = [_]bool{ false, false, false };
    for (0..3) |i| cpick[i] = @abs(sc.t.at(i, i).re - 2) < 1e-12;
    const oc = try ordschur(al, io, sc.u, sc.t, .{ .mask = &cpick });
    defer oc.deinit(al);
    try testing.expectApproxEqAbs(2, oc.t.at(0, 0).re, 1e-13);
    try testing.expectApproxEqAbs(0, oc.t.at(0, 0).im, 1e-13);

    try testing.expectError(error.DimensionMismatch, ordschur(al, io, id, v, .{ .mask = &.{ true, false } }));
    try testing.expectError(error.DimensionMismatch, ordschur(al, io, id, v, .{ .clusters = &.{ 1, 2, 3 } }));
}

test "ordschur: region keywords, MATLAB's boundaries" {
    const al = testing.allocator;
    const io = testing.io;
    // eigenvalues 0.5, -2, 0 and 1 on the diagonal
    const a = [_]f64{ 0.5, 0, 0, 0, 1, -2, 0, 0, 0.5, 0.25, 0, 0, -1, 2, 0.75, 1 };
    const v: ConstView(f64) = try .init(&a, 4, 4, 4);
    const id: ConstView(f64) = try .init(&eye4, 4, 4, 4);
    const cases = [_]struct { r: Region, lead: []const f64 }{
        .{ .r = .lhp, .lead = &.{-2} },
        .{ .r = .rhp, .lead = &.{ 0.5, 0, 1 } },
        .{ .r = .udi, .lead = &.{ 0.5, 0 } },
        .{ .r = .udo, .lead = &.{ -2, 1 } },
    };
    for (cases) |c| {
        const s = try ordschur(al, io, id, v, .{ .region = c.r });
        defer s.deinit(al);
        try testing.expect(residual(&a, s.u, s.t) < 1e-13 * 8);
        for (c.lead, 0..) |e, i| try testing.expectApproxEqAbs(e, s.t.at(i, i), 1e-13);
    }

    // a real pair 0.5 ± 2i (modulus above 1) behind -0.25: udo picks the pair
    const p = [_]f64{ -0.25, 0, 0, 1, 0.5, -2, 0.5, 2, 0.5 };
    const pv = try ConstView(f64).init(&p, 3, 3, 3);
    const eye3 = [_]f64{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    const s = try ordschur(al, io, try ConstView(f64).init(&eye3, 3, 3, 3), pv, .{ .region = .udo });
    defer s.deinit(al);
    try testing.expect(residual(&p, s.u, s.t) < 1e-13 * 8);
    try testing.expect(s.t.at(1, 0) != 0);
    try testing.expectApproxEqAbs(-0.25, s.t.at(2, 2), 1e-13);
}

test "ordschur: clusters, the highest index first" {
    const al = testing.allocator;
    const io = testing.io;
    const a = [_]f64{ 1, 0, 0, 0, 0.5, 2, 0, 0, -1, 0.25, 3, 0, 2, 1, -0.5, 4 };
    const v: ConstView(f64) = try .init(&a, 4, 4, 4);
    const id: ConstView(f64) = try .init(&eye4, 4, 4, 4);
    // 1 and 3 in cluster 2, 4 in cluster 5, 2 in cluster 0: 4, then 1 and 3, then 2
    const s = try ordschur(al, io, id, v, .{ .clusters = &.{ 2, 0, 2, 5 } });
    defer s.deinit(al);
    try testing.expect(residual(&a, s.u, s.t) < 1e-13 * 8);
    for ([_]f64{ 4, 1, 3, 2 }, 0..) |e, i| try testing.expectApproxEqAbs(e, s.t.at(i, i), 1e-13);

    // one cluster: nothing moves
    const one = try ordschur(al, io, id, v, .{ .clusters = &.{ 7, 7, 7, 7 } });
    defer one.deinit(al);
    for (0..4) |i| try testing.expectEqual(a[i + i * 4], one.t.at(i, i));

    // a real pair takes its rows' higher index: it leads, then 2
    const p = [_]f64{ 2, 0, 0, 1, 0.5, -2, 0.5, 2, 0.5 };
    const pv = try ConstView(f64).init(&p, 3, 3, 3);
    const eye3 = [_]f64{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    const sp = try ordschur(al, io, try ConstView(f64).init(&eye3, 3, 3, 3), pv, .{ .clusters = &.{ 1, 0, 3 } });
    defer sp.deinit(al);
    try testing.expect(residual(&p, sp.u, sp.t) < 1e-13 * 8);
    try testing.expect(sp.t.at(1, 0) != 0);
    try testing.expectApproxEqAbs(2, sp.t.at(2, 2), 1e-13);
}
