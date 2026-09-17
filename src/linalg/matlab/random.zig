//! Uniform/normal random number generation, ported from `matlab88.f90`'s
//! `URAND` (a portable linear-congruential generator after D.E. Knuth,
//! *The Art of Computer Programming* vol. 2) and the `RAND` command body
//! in `MATFN6` (uniform fill, and the Marsaglia polar method for a normal
//! fill).
//!
//! `URAND` is a 1970s-vintage generator, carried here for REPL-command
//! parity (`RAND(N)` and friends) and reproducible sequences, not because
//! it's a good source of randomness by modern standards - reach for
//! `std.Random` elsewhere if statistical quality matters.
//!
//! The original self-detected the host machine's integer word length at
//! runtime (many machines of the era had different `INTEGER` widths); this
//! port fixes that at the canonical 32-bit word size the constants below
//! were derived for, computing the linear-congruential step in a wider
//! (`i64`) accumulator so the same overflow-correction logic the original
//! spelled out by hand (rather than relying on wraparound, which Fortran
//! did not guarantee portably) runs without tripping Zig's overflow checks.

const std = @import("std");
const zigstd = @import("../../zigstd.zig");

/// URAND's generator state (the original's `IY` parameter).
pub const Urand = struct {
    // A genuine 32-bit `INTEGER`, matching the original: this generator
    // relies on `IY*IA` and `IY+IC` silently wrapping around on overflow,
    // which is what every real machine the original ran on did for
    // `INTEGER` arithmetic (undefined by the Fortran standard, universal
    // in practice) - so the multiply/add below use Zig's explicit
    // wrapping operators (`*%`/`+%`) to reproduce that on purpose, with
    // the same manual corrections the original layered on top for
    // portability to hardware whose overflow behavior differed slightly.
    iy: i32,

    // The word-length self-detection loop ("double M until doubling stops
    // increasing it") finds the largest power of 2 a 32-bit *signed*
    // INTEGER can hold before the next doubling overflows the sign bit:
    // 2^30, not 2^31.
    const m2: i32 = 1 << 30;
    const half_m: f64 = @floatFromInt(m2);
    // atan(1.0) == pi/4, written directly since std.math.atan won't
    // const-evaluate here.
    const ia: i32 = 8 * @as(i32, @intFromFloat(@floor(half_m * (std.math.pi / 4.0) / 8.0))) + 5;
    const ic: i32 = 2 * @as(i32, @intFromFloat(@floor(half_m * (0.5 - @sqrt(3.0) / 6.0)))) + 1;
    const mic: i32 = (m2 -% ic) +% m2;
    const scale: f64 = 0.5 / half_m;

    /// Seed with an arbitrary nonzero integer (matching the original's
    /// requirement that `IY` start as "an arbitrary integer").
    pub fn init(seed: i32) Urand {
        return .{ .iy = seed };
    }

    /// URAND: the next uniform value in `(0, 1)`.
    pub fn next(self: *Urand) f64 {
        var iy: i32 = self.iy *% ia;
        if (iy > mic) iy = (iy -% m2) -% m2;
        iy +%= ic;
        if (@divTrunc(iy, 2) > m2) iy = (iy -% m2) -% m2;
        if (iy < 0) iy = (iy +% m2) +% m2;
        self.iy = iy;
        return @as(f64, @floatFromInt(iy)) * scale;
    }

    /// RAND's normal-fill formula: reject-sample a point uniform in the
    /// unit disk, then rescale it by the Marsaglia polar transform. Only
    /// one of the transform's two independent normal deviates is used per
    /// call, matching the original (which drew a fresh disk point per
    /// matrix element rather than caching the second deviate).
    pub fn nextNormal(self: *Urand) f64 {
        while (true) {
            const sr = 2.0 * self.next() - 1.0;
            const si = 2.0 * self.next() - 1.0;
            const t = sr * sr + si * si;
            if (t <= 1.0) return sr * @sqrt(-2.0 * @log(t) / t);
        }
    }

    /// RAND(M, N): fill an `m`x`n` real matrix (column-major, leading
    /// dimension `m`) with independent uniform values in `(0, 1)`.
    pub fn fillUniform(self: *Urand, a: []f64, m: usize, n: usize) void {
        for (a[0 .. m * n]) |*v| v.* = self.next();
    }

    /// RAND(M, N) in normal mode: fill an `m`x`n` real matrix (column-major,
    /// leading dimension `m`) with independent standard-normal values.
    pub fn fillNormal(self: *Urand, a: []f64, m: usize, n: usize) void {
        for (a[0 .. m * n]) |*v| v.* = self.nextNormal();
    }
};


test "Urand produces a deterministic, repeatable, bounded sequence" {
    var a = Urand.init(12345);
    var b = Urand.init(12345);
    for (0..20) |_| {
        const va = a.next();
        const vb = b.next();
        try zigstd.testExpectEqual(va, vb);
        try zigstd.testExpectTrue(va > 0.0 and va < 1.0);
    }
}

test "Urand.next has roughly the right mean over many draws" {
    var rng = Urand.init(1);
    var sum: f64 = 0.0;
    const n = 20000;
    for (0..n) |_| sum += rng.next();
    const mean = sum / @as(f64, n);
    try zigstd.testExpectApproxEqAbs(0.5, mean, 0.02);
}

test "Urand.nextNormal has roughly zero mean and unit variance" {
    var rng = Urand.init(7);
    var sum: f64 = 0.0;
    var sum2: f64 = 0.0;
    const n = 20000;
    for (0..n) |_| {
        const v = rng.nextNormal();
        sum += v;
        sum2 += v * v;
    }
    const mean = sum / @as(f64, n);
    const variance = sum2 / @as(f64, n) - mean * mean;
    try zigstd.testExpectApproxEqAbs(0.0, mean, 0.05);
    try zigstd.testExpectApproxEqAbs(1.0, variance, 0.1);
}
