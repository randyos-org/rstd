//! `TrigTransformPlan` and `ModifiedCosineTransform` against their
//! definitions (SciPy's `dct` and `dst` docs, unnormalized), orthogonal
//! under `.ortho`, inverted as documented, in batches.

const std = @import("std");
const fft = @import("../__root__.zig");

const TrigTransformPlan = fft.TrigTransformPlan;
const Kind = fft.TrigTransformKind;
const ModifiedCosineTransform = fft.ModifiedCosineTransform;

const testing = std.testing;
const pi = std.math.pi;

const kinds = std.enums.values(Kind);

/// A test signal with no symmetry to hide behind.
fn signal(comptime T: type, out: []T, seed: f64) void {
    for (out, 0..) |*x, i| {
        const t: f64 = @floatFromInt(i);
        x.* = @floatCast(@sin(seed + 0.7 * t) + 0.25 * @cos(0.3 * t * t));
    }
}

/// `kind`'s unnormalized transform of `x` by its definition, in f64.
fn direct(kind: Kind, out: []f64, x: []const f64) void {
    const n = x.len;
    const nf: f64 = @floatFromInt(n);
    for (out, 0..) |*o, k| {
        const kf: f64 = @floatFromInt(k);
        var s: f64 = 0;
        for (x, 0..) |v, j| {
            const jf: f64 = @floatFromInt(j);
            s += v * switch (kind) {
                .dct1 => if (j == 0) 1 else if (j == n - 1) (if (k % 2 == 0) @as(f64, 1) else -1) else 2 * @cos(pi * kf * jf / (nf - 1)),
                .dct2 => 2 * @cos(pi * kf * (2 * jf + 1) / (2 * nf)),
                .dct3 => if (j == 0) 1 else 2 * @cos(pi * (2 * kf + 1) * jf / (2 * nf)),
                .dct4 => 2 * @cos(pi * (2 * kf + 1) * (2 * jf + 1) / (4 * nf)),
                .dst1 => 2 * @sin(pi * (kf + 1) * (jf + 1) / (nf + 1)),
                .dst2 => 2 * @sin(pi * (kf + 1) * (2 * jf + 1) / (2 * nf)),
                .dst3 => if (j == n - 1) (if (k % 2 == 0) @as(f64, 1) else -1) else 2 * @sin(pi * (2 * kf + 1) * (jf + 1) / (2 * nf)),
                .dst4 => 2 * @sin(pi * (2 * kf + 1) * (2 * jf + 1) / (4 * nf)),
            };
        }
        o.* = s;
    }
}

fn norm(x: anytype) f64 {
    var s: f64 = 0;
    for (x) |v| s += @as(f64, v) * @as(f64, v);
    return @sqrt(s);
}

fn checkKind(comptime T: type, kind: Kind, n: usize) !void {
    const al = testing.allocator;
    if (kind == .dct1 and n < 2) {
        try testing.expectError(error.LengthTooShort, TrigTransformPlan(T).init(al, n, kind));
        return;
    }
    var plan = try TrigTransformPlan(T).init(al, n, kind);
    defer plan.deinit(al);
    const scratch = try al.alloc(T, plan.scratchLen());
    defer al.free(scratch);
    const x = try al.alloc(T, n);
    defer al.free(x);
    const y = try al.alloc(T, n);
    defer al.free(y);
    const x64 = try al.alloc(f64, n);
    defer al.free(x64);
    const want = try al.alloc(f64, n);
    defer al.free(want);
    signal(T, x, 1);
    for (x64, x) |*a, b| a.* = b;
    direct(kind, want, x64);
    @memcpy(y, x);
    try plan.transform(testing.io, y, scratch, .{}, .backward);
    const eps: f64 = std.math.floatEps(T);
    const nf: f64 = @floatFromInt(n);
    // the transform's rounding, and the direct sums' in f64
    var sum: f64 = 0;
    for (x64) |v| sum += @abs(v);
    const tol = 32 * eps * (@log2(nf) + 2) * norm(x64) + 4 * nf * std.math.floatEps(f64) * sum;
    for (want, y, 0..) |a, b, k| if (@abs(a - b) > tol) {
        std.debug.print("{s} {t} n={d} k={d}: {d} against {d}\n", .{ @typeName(T), kind, n, k, b, a });
        return error.TestUnexpectedResult;
    };
    // orthogonal under `.ortho`
    @memcpy(y, x);
    try plan.transform(testing.io, y, scratch, .{}, .ortho);
    try testing.expectApproxEqRel(norm(x), norm(y), 64 * eps);
    // the documented inverse under `.forward`
    const inv: Kind = switch (kind) {
        .dct2 => .dct3,
        .dct3 => .dct2,
        .dst2 => .dst3,
        .dst3 => .dst2,
        else => kind,
    };
    var back = try TrigTransformPlan(T).init(al, n, inv);
    defer back.deinit(al);
    const s2 = try al.alloc(T, back.scratchLen());
    defer al.free(s2);
    @memcpy(y, x);
    try plan.transform(testing.io, y, scratch, .{}, .backward);
    try back.transform(testing.io, y, s2, .{}, .forward);
    for (x, y) |a, b| try testing.expectApproxEqAbs(@as(f64, a), @as(f64, b), 64 * eps * (norm(x) + 1));
}

test "TrigTransformPlan matches each kind's definition, n = 1 to 40 and past, f32 and f64; ortho is orthogonal, the inverse as documented" {
    inline for (.{ f32, f64 }) |T| for (kinds) |kind| {
        for (1..41) |n| try checkKind(T, kind, n);
        // pffft's, Bluestein's and the generic pass's lengths
        for ([_]usize{ 64, 96, 97, 128, 210, 1009 }) |n| try checkKind(T, kind, n);
    };
}

test "TrigTransformPlan: a strided batch has each vector's bits alone, and short slices" {
    const al = testing.allocator;
    const n = 12;
    const count = 3;
    for (kinds) |kind| {
        var plan = try TrigTransformPlan(f64).init(al, n, kind);
        defer plan.deinit(al);
        const scratch = try al.alloc(f64, 2 * plan.scratchLen());
        defer al.free(scratch);
        // interleaved channels
        var data: [n * count]f64 = undefined;
        signal(f64, &data, 2);
        var want: [n * count]f64 = undefined;
        for (0..count) |c| {
            var v: [n]f64 = undefined;
            for (&v, 0..) |*e, i| e.* = data[i * count + c];
            try plan.transform(testing.io, &v, scratch, .{}, .ortho);
            for (v, 0..) |e, i| want[i * count + c] = e;
        }
        try plan.transform(testing.io, &data, scratch, .{ .count = count, .stride = count, .distance = 1 }, .ortho);
        try testing.expectEqualSlices(f64, &want, &data);
        try testing.expectError(error.SliceTooShort, plan.transform(testing.io, &data, scratch[0 .. plan.scratchLen() - 1], .{}, .ortho));
        try testing.expectError(error.SliceTooShort, plan.transform(testing.io, data[0 .. n - 1], scratch, .{}, .ortho));
    }
}

test "ModifiedCosineTransform matches its definition and its inverse's; odd lengths refused" {
    const al = testing.allocator;
    try testing.expectError(error.OddLength, ModifiedCosineTransform(f64).init(al, 7));
    inline for (.{ f32, f64 }) |T| for ([_]usize{ 2, 4, 6, 10, 16, 48, 100 }) |n| {
        var m = try ModifiedCosineTransform(T).init(al, n);
        defer m.deinit(al);
        const scratch = try al.alloc(T, m.scratchLen());
        defer al.free(scratch);
        const x = try al.alloc(T, 2 * n);
        defer al.free(x);
        const big = try al.alloc(T, 2 * n);
        defer al.free(big);
        const coef = try al.alloc(T, n);
        defer al.free(coef);
        signal(T, x, 3);
        try m.forward(coef, x, scratch);
        const nf: f64 = @floatFromInt(n);
        const eps: f64 = std.math.floatEps(T);
        const tol = 32 * eps * (@log2(nf) + 2) * norm(x);
        for (coef, 0..) |got, k| {
            var s: f64 = 0;
            for (x, 0..) |v, j| s += @as(f64, v) * @cos(pi / nf * (@as(f64, @floatFromInt(j)) + 0.5 + nf / 2) * (@as(f64, @floatFromInt(k)) + 0.5));
            try testing.expectApproxEqAbs(s, @as(f64, got), tol);
        }
        try m.backward(big, coef, scratch);
        for (big, 0..) |got, j| {
            var s: f64 = 0;
            for (coef, 0..) |v, k| s += @as(f64, v) * @cos(pi / nf * (@as(f64, @floatFromInt(j)) + 0.5 + nf / 2) * (@as(f64, @floatFromInt(k)) + 0.5));
            try testing.expectApproxEqAbs(2 * s / nf, @as(f64, got), tol);
        }
        try testing.expectError(error.SliceTooShort, m.forward(coef[0 .. n - 1], x, scratch));
        try testing.expectError(error.SliceTooShort, m.backward(big[0 .. 2 * n - 1], coef, scratch));
    };
}
