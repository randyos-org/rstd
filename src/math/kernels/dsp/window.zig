//! Window functions: Hann, Hamming and Blackman as MATLAB's `hann`,
//! `hamming` and `blackman`, Kaiser as its `kaiser`; the sine window and the
//! Kaiser-Bessel-derived one, `ModifiedCosineTransform`'s Princen-Bradley windows, as
//! SciPy's `windows.cosine` and `windows.kaiser_bessel_derived`. Computed in
//! f64 (wider types in their own), written to the caller's slice.

const std = @import("std");

/// MATLAB's `'symmetric'` (filter design) and `'periodic'` (spectral
/// analysis: the symmetric window one longer, its last point dropped).
pub const Symmetry = enum { symmetric, periodic };

pub const Window = union(enum) {
    hann,
    hamming,
    blackman,
    /// Shape `beta` (MATLAB's default 0.5).
    kaiser: f64,
    /// `sin(π (n + ½) / L)`.
    sine,
    /// Kaiser-Bessel-derived, shape `beta` (AAC's `π α`, α 4 or 6): even
    /// lengths, symmetric only.
    kbd: f64,
};

pub const Error = error{
    /// KBD of an odd length.
    OddLength,
    /// KBD asked for periodic.
    UnsupportedSymmetry,
};

/// `out[n]` for `n < out.len`; an empty `out` is nothing to do, one point
/// is 1 (MATLAB's).
pub fn window(comptime T: type, out: []T, kind: Window, symmetry: Symmetry) Error!void {
    comptime if (@typeInfo(T) != .float) @compileError("rstd.math.kernels.dsp: " ++ @typeName(T) ++ " is not a float type");
    const H = if (@bitSizeOf(T) > 64) T else f64;
    if (kind == .kbd) {
        if (symmetry == .periodic) return error.UnsupportedSymmetry;
        if (out.len % 2 != 0) return error.OddLength;
    }
    if (out.len == 0) return;
    if (out.len == 1) {
        out[0] = 1;
        return;
    }
    // the symmetric window of `len` points, the periodic one's first `out.len`
    const len = if (symmetry == .periodic) out.len + 1 else out.len;
    const m: H = @floatFromInt(len - 1);
    const pi = std.math.pi;
    switch (kind) {
        .hann => for (out, 0..) |*w, n| {
            w.* = @floatCast(0.5 - 0.5 * @cos(2 * pi * f(H, n) / m));
        },
        .hamming => for (out, 0..) |*w, n| {
            w.* = @floatCast(0.54 - 0.46 * @cos(2 * pi * f(H, n) / m));
        },
        .blackman => for (out, 0..) |*w, n| {
            const x = 2 * pi * f(H, n) / m;
            w.* = @floatCast(0.42 - 0.5 * @cos(x) + 0.08 * @cos(2 * x));
        },
        .kaiser => |beta| {
            const b: H = @floatCast(beta);
            const den = besselI0(H, b);
            for (out, 0..) |*w, n| {
                const r = 2 * f(H, n) / m - 1;
                w.* = @floatCast(besselI0(H, b * @sqrt(@max(0, 1 - r * r))) / den);
            }
        },
        .sine => for (out, 0..) |*w, n| {
            w.* = @floatCast(@sin(pi * (f(H, n) + 0.5) / f(H, len)));
        },
        .kbd => |beta| {
            // the cumulative Kaiser window of `len / 2 + 1` points, square-rooted
            // and mirrored
            const half = out.len / 2;
            const b: H = @floatCast(beta);
            const hm: H = @floatFromInt(half);
            var total: H = 0;
            for (0..half + 1) |n| total += kaiserPoint(H, b, n, hm);
            var csum: H = 0;
            for (0..half) |n| {
                csum += kaiserPoint(H, b, n, hm);
                const v: T = @floatCast(@sqrt(csum / total));
                out[n] = v;
                out[out.len - 1 - n] = v;
            }
        },
    }
}

/// Point `n` of the symmetric Kaiser window over `0..m`, unnormalized.
fn kaiserPoint(comptime H: type, beta: H, n: usize, m: H) H {
    const r = 2 * f(H, n) / m - 1;
    return besselI0(H, beta * @sqrt(@max(0, 1 - r * r)));
}

/// `I0(x) = Σ ((x / 2)^k / k!)²`, summed until a term no longer changes it.
fn besselI0(comptime H: type, x: H) H {
    const q = x * x / 4;
    var term: H = 1;
    var sum: H = 1;
    var k: H = 1;
    while (true) : (k += 1) {
        term *= q / (k * k);
        const next = sum + term;
        if (next == sum) return sum;
        sum = next;
    }
}

inline fn f(comptime H: type, n: usize) H {
    return @floatFromInt(n);
}

// ---- tests ----

const testing = std.testing;

fn expectWindow(want: []const f64, kind: Window, symmetry: Symmetry) !void {
    var got: [16]f64 = undefined;
    try window(f64, got[0..want.len], kind, symmetry);
    for (want, got[0..want.len]) |a, b| try testing.expectApproxEqAbs(a, b, 1e-4);
}

test "window: MATLAB's and SciPy's values" {
    try expectWindow(&.{ 0, 0.5, 1, 0.5, 0 }, .hann, .symmetric);
    try expectWindow(&.{ 0, 0.5, 1, 0.5 }, .hann, .periodic);
    try expectWindow(&.{ 0.08, 0.54, 1, 0.54, 0.08 }, .hamming, .symmetric);
    try expectWindow(&.{ 0, 0.34, 1, 0.34, 0 }, .blackman, .symmetric);
    try expectWindow(&.{ 1, 1, 1, 1, 1 }, .{ .kaiser = 0 }, .symmetric);
    try expectWindow(&.{ 0.0367, 0.5529, 1, 0.5529, 0.0367 }, .{ .kaiser = 5 }, .symmetric);
    try expectWindow(&.{ 0.3827, 0.9239, 0.9239, 0.3827 }, .sine, .symmetric);
    try expectWindow(&.{1}, .hann, .periodic);
    var none: [0]f64 = .{};
    try window(f64, &none, .hann, .symmetric);
}

test "window: the sine and KBD windows meet Princen-Bradley, w[n]² + w[n + N]² = 1" {
    inline for (.{ f16, f32, f64, f80, f128 }) |T| {
        // checked in f64
        const tol: f64 = 64 * @max(@as(f64, @floatCast(std.math.floatEps(T))), std.math.floatEps(f64));
        for ([_]Window{ .sine, .{ .kbd = 4 * std.math.pi }, .{ .kbd = 6 * std.math.pi } }) |kind| {
            var w: [32]T = undefined;
            try window(T, &w, kind, .symmetric);
            for (0..16) |n| {
                const a: f64 = @floatCast(w[n]);
                const b: f64 = @floatCast(w[n + 16]);
                try testing.expectApproxEqAbs(1, a * a + b * b, tol);
                try testing.expectApproxEqAbs(a, @as(f64, @floatCast(w[31 - n])), tol);
            }
        }
    }
}

test "window: KBD takes even lengths, symmetric only" {
    var w: [5]f32 = undefined;
    try testing.expectError(error.OddLength, window(f32, &w, .{ .kbd = 4 }, .symmetric));
    try testing.expectError(error.UnsupportedSymmetry, window(f32, w[0..4], .{ .kbd = 4 }, .periodic));
}
