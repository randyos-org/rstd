//! `fft.TrigTransformPlan`: the DCT and DST of types 1 to 4, PocketFFT's
//! `T_dct1`, `T_dst1`, `T_dcst23` and `T_dcst4` (`pocketfft_hdronly.h`, cpp
//! branch, commit 4d865cba, BSD-3, LICENSE) over `RealPlan`'s halfcomplex
//! transform, an even DCT-IV over `Plan`'s; PocketFFT's odd DCT-IV is FFTW3's
//! `apply_re11`, BSD-3 by its authors' leave. Scaled as SciPy's `dct` and
//! `dst` (`pypocketfft.cxx`'s `norm_fct`).

const std = @import("std");
const fft = @import("__root__.zig");
const cmplx = @import("cmplx.zig");
const complex_plan = @import("plan.zig");
const real_plan = @import("real_plan.zig");
const twiddle = @import("twiddle.zig");
const batch = @import("batch.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const C = cmplx.C;
const BatchLayout = fft.BatchLayout;
const Normalization = fft.Normalization;

/// DCT and DST types 1 to 4, as SciPy's `dct` and `dst`.
pub const Kind = enum { dct1, dct2, dct3, dct4, dst1, dst2, dst3, dst4 };

/// `fft.InitError`, and `LengthTooShort` for a DCT-I of one point, which
/// SciPy refuses too.
pub const InitError = fft.InitError || error{LengthTooShort};

const sqrt2 = 1.414213562373095048801688724209698;

pub fn TrigTransformPlan(comptime T: type) type {
    if (@typeInfo(T) != .float) @compileError("rstd.math.kernels.fft: " ++ @typeName(T) ++ " is not a float type");
    return struct {
        n: usize,
        kind: Kind,
        scratch_len: usize,
        inner: Inner,
        /// `T_dcst23`'s twiddles (`n`), else empty.
        twiddle: []T,
        /// `T_dcst4`'s `C2` for an even `n` (`n / 2`), else empty.
        c2: []C(T),

        const Self = @This();
        const Inner = union(enum) { real: real_plan.RealPlan(T), complex: complex_plan.Plan(T) };

        pub fn init(allocator: Allocator, n: usize, kind: Kind) InitError!Self {
            if (n == 0) return error.ZeroLength;
            var s: Self = .{ .n = n, .kind = kind, .scratch_len = 0, .inner = undefined, .twiddle = &.{}, .c2 = &.{} };
            switch (kind) {
                .dct1 => {
                    if (n < 2) return error.LengthTooShort;
                    s.inner = .{ .real = try .init(allocator, 2 * (n - 1)) };
                },
                .dst1 => s.inner = .{ .real = try .init(allocator, 2 * (n + 1)) },
                .dct2, .dct3, .dst2, .dst3 => {
                    s.inner = .{ .real = try .init(allocator, n) };
                    errdefer s.inner.real.deinit(allocator);
                    var tw = try twiddle.Twiddles(T).init(allocator, 4 * n);
                    defer tw.deinit(allocator);
                    s.twiddle = try allocator.alloc(T, n);
                    for (s.twiddle, 1..) |*t, i| t.* = tw.at(i).re;
                },
                .dct4, .dst4 => if (n % 2 == 1) {
                    s.inner = .{ .real = try .init(allocator, n) };
                } else {
                    s.inner = .{ .complex = try .init(allocator, n / 2) };
                    errdefer s.inner.complex.deinit(allocator);
                    var tw = try twiddle.Twiddles(T).init(allocator, 16 * n);
                    defer tw.deinit(allocator);
                    s.c2 = try allocator.alloc(C(T), n / 2);
                    for (s.c2, 0..) |*c, i| {
                        const w = tw.at(8 * i + 1);
                        c.* = .{ .re = w.re, .im = -w.im };
                    }
                },
            }
            // the gathered data, then the inner plan's buffer and space
            s.scratch_len = n + switch (s.inner) {
                .real => |*r| r.scratchLen(),
                .complex => |*c| 2 * c.scratchLen(),
            };
            return s;
        }

        pub fn deinit(self: *Self, allocator: Allocator) void {
            switch (self.inner) {
                inline else => |*p| p.deinit(allocator),
            }
            allocator.free(self.twiddle);
            allocator.free(self.c2);
            self.* = undefined;
        }

        /// Scratch per transform, as `Plan`'s.
        pub fn scratchLen(self: Self) usize {
            return self.scratch_len;
        }

        /// The kind's transform of each of `layout`'s vectors, in place.
        /// `.backward` scales none, `.forward` divides by `N` and `.ortho`
        /// by `√N`, weighting the ends so the matrix is orthogonal, `N` being
        /// `2 (n - 1)` for the DCT-I, `2 (n + 1)` for the DST-I, else `2 n`:
        /// a type 2's inverse is the type 3 under `.forward` (or both
        /// `.ortho`), and the reverse; a type 1's or 4's, itself.
        pub fn transform(self: *const Self, io: Io, data: []T, scratch: []T, layout: BatchLayout, normalization: Normalization) fft.TransformError!void {
            if (scratch.len < self.scratch_len) return error.SliceTooShort;
            try complex_plan.covers(data.len, self.n, layout);
            const Ctx = struct {
                plan: *const Self,
                data: []T,
                scratch: []T,
                stride: usize,
                dist: usize,
                f: T,
                ortho: bool,

                fn one(c: *const @This(), k: usize, slot: usize) void {
                    const p = c.plan;
                    const s = c.scratch[slot * p.scratch_len ..][0..p.scratch_len];
                    const v = s[0..p.n];
                    const d = c.data[k * c.dist ..];
                    for (v, 0..) |*x, i| x.* = d[i * c.stride];
                    p.exec(v, s[p.n..], c.f, c.ortho);
                    for (v, 0..) |x, i| d[i * c.stride] = x;
                }
            };
            const ctx: Ctx = .{
                .plan = self,
                .data = data,
                .scratch = scratch,
                .stride = layout.stride,
                .dist = layout.distance orelse (self.n * layout.stride),
                .f = self.factor(normalization),
                .ortho = normalization == .ortho,
            };
            try batch.run(io, layout.count, scratch.len / self.scratch_len, &ctx, Ctx.one);
        }

        /// `norm_fct`'s scale for `normalization`.
        fn factor(self: *const Self, normalization: Normalization) T {
            const big: T = @floatFromInt(switch (self.kind) {
                .dct1 => 2 * (self.n - 1),
                .dst1 => 2 * (self.n + 1),
                else => 2 * self.n,
            });
            return switch (normalization) {
                .backward => 1,
                .forward => 1 / big,
                .ortho => 1 / @sqrt(big),
            };
        }

        /// One transform of the contiguous `c` (`n`), scaled by `f`, `rest`
        /// (at least `scratchLen() - n`) its space: PocketFFT's `exec`.
        pub fn exec(self: *const Self, c: []T, rest: []T, f: T, ortho: bool) void {
            switch (self.kind) {
                .dct1 => self.dct1(c, rest, f, ortho),
                .dst1 => self.dst1(c, rest, f),
                .dct2 => self.dcst23(c, rest, f, ortho, 2, true),
                .dct3 => self.dcst23(c, rest, f, ortho, 3, true),
                .dst2 => self.dcst23(c, rest, f, ortho, 2, false),
                .dst3 => self.dcst23(c, rest, f, ortho, 3, false),
                .dct4 => self.dcst4(c, rest, f, true),
                .dst4 => self.dcst4(c, rest, f, false),
            }
        }

        /// The halfcomplex transform of `tmp` (the plan's length) by the
        /// real plan, `rest` after it.
        fn rfft(self: *const Self, tmp: []T, rest: []T, f: T) void {
            self.inner.real.run(tmp, rest, f, true);
        }

        /// `T_dct1::exec`.
        fn dct1(self: *const Self, c: []T, rest: []T, f: T, ortho: bool) void {
            const big = self.inner.real.n;
            const n = big / 2 + 1;
            const s2: T = sqrt2;
            if (ortho) {
                c[0] *= s2;
                c[n - 1] *= s2;
            }
            const tmp = rest[0..big];
            tmp[0] = c[0];
            for (1..n) |i| {
                tmp[i] = c[i];
                tmp[big - i] = c[i];
            }
            self.rfft(tmp, rest[big..], f);
            c[0] = tmp[0];
            for (1..n) |i| c[i] = tmp[2 * i - 1];
            if (ortho) {
                c[0] *= s2 * 0.5;
                c[n - 1] *= s2 * 0.5;
            }
        }

        /// `T_dst1::exec`.
        fn dst1(self: *const Self, c: []T, rest: []T, f: T) void {
            const big = self.inner.real.n;
            const n = big / 2 - 1;
            const tmp = rest[0..big];
            tmp[0] = 0;
            tmp[n + 1] = 0;
            for (0..n) |i| {
                tmp[i + 1] = c[i];
                tmp[big - 1 - i] = -c[i];
            }
            self.rfft(tmp, rest[big..], f);
            for (0..n) |i| c[i] = -tmp[2 * i + 2];
        }

        /// `T_dcst23::exec`: types 2 and 3, cosine or sine.
        fn dcst23(self: *const Self, c: []T, rest: []T, f: T, ortho: bool, comptime kind: u2, comptime cosine: bool) void {
            const big = self.n;
            const ns2 = (big + 1) / 2;
            const tw = self.twiddle;
            const s2: T = sqrt2;
            if (kind == 2) {
                if (!cosine) {
                    var k: usize = 1;
                    while (k < big) : (k += 2) c[k] = -c[k];
                }
                c[0] *= 2;
                if (big % 2 == 0) c[big - 1] *= 2;
                var k: usize = 1;
                while (k + 1 < big) : (k += 2) mpInPlace(&c[k + 1], &c[k]);
                self.inner.real.run(c, rest, f, false);
                k = 1;
                var kc = big - 1;
                while (k < ns2) : ({
                    k += 1;
                    kc -= 1;
                }) {
                    const t1 = tw[k - 1] * c[kc] + tw[kc - 1] * c[k];
                    const t2 = tw[k - 1] * c[k] - tw[kc - 1] * c[kc];
                    c[k] = 0.5 * (t1 + t2);
                    c[kc] = 0.5 * (t1 - t2);
                }
                if (big % 2 == 0) c[ns2] *= tw[ns2 - 1];
                if (!cosine) std.mem.reverse(T, c[0..big]);
                if (ortho) {
                    if (cosine) c[0] *= s2 * 0.5 else c[big - 1] *= s2 * 0.5;
                }
            } else {
                if (ortho) {
                    if (cosine) c[0] *= s2 else c[big - 1] *= s2;
                }
                // `k < NS2` swaps: a reversal, an odd middle swapped with itself
                if (!cosine) std.mem.reverse(T, c[0..big]);
                var k: usize = 1;
                var kc = big - 1;
                while (k < ns2) : ({
                    k += 1;
                    kc -= 1;
                }) {
                    const t1 = c[k] + c[kc];
                    const t2 = c[k] - c[kc];
                    c[k] = tw[k - 1] * t2 + tw[kc - 1] * t1;
                    c[kc] = tw[k - 1] * t1 - tw[kc - 1] * t2;
                }
                if (big % 2 == 0) c[ns2] *= 2 * tw[ns2 - 1];
                self.inner.real.run(c, rest, f, true);
                k = 1;
                while (k + 1 < big) : (k += 2) mpInPlace(&c[k], &c[k + 1]);
                if (!cosine) {
                    k = 1;
                    while (k < big) : (k += 2) c[k] = -c[k];
                }
            }
        }

        /// `T_dcst4::exec`: the type 4, cosine or sine.
        fn dcst4(self: *const Self, c: []T, rest: []T, f: T, cosine: bool) void {
            const big = self.n;
            // `k < n2` swaps: a reversal, an odd middle staying
            if (!cosine) std.mem.reverse(T, c[0..big]);
            if (big % 2 == 1) self.odd4(c, rest, f) else self.even4(c, rest, f);
            if (!cosine) {
                var k: usize = 1;
                while (k < big) : (k += 2) c[k] = -c[k];
            }
        }

        /// The odd DCT-IV: FFTW3's `apply_re11`, as PocketFFT has it.
        fn odd4(self: *const Self, c: []T, rest: []T, f: T) void {
            const big = self.n;
            const n2 = big / 2;
            const y = rest[0..big];
            {
                var i: usize = 0;
                var m = n2;
                while (m < big) : ({
                    i += 1;
                    m += 4;
                }) y[i] = c[m];
                while (m < 2 * big) : ({
                    i += 1;
                    m += 4;
                }) y[i] = -c[2 * big - m - 1];
                while (m < 3 * big) : ({
                    i += 1;
                    m += 4;
                }) y[i] = -c[m - 2 * big];
                while (m < 4 * big) : ({
                    i += 1;
                    m += 4;
                }) y[i] = c[4 * big - m - 1];
                while (i < big) : ({
                    i += 1;
                    m += 4;
                }) y[i] = c[m - 4 * big];
            }
            self.rfft(y, rest[big..], f);
            c[n2] = y[0] * sgn(T, n2 + 1);
            var i: usize = 0;
            var ia: usize = 1;
            var k: usize = 1;
            while (k < n2) : ({
                i += 1;
                ia += 1;
                k += 2;
            }) {
                c[i] = y[2 * k - 1] * sgn(T, ia) + y[2 * k] * sgn(T, i);
                c[big - ia] = y[2 * k - 1] * sgn(T, big - i) - y[2 * k] * sgn(T, big - ia);
                c[n2 - ia] = y[2 * k + 1] * sgn(T, n2 - i) - y[2 * k + 2] * sgn(T, n2 - ia);
                c[n2 + ia] = y[2 * k + 1] * sgn(T, n2 + i + 2) + y[2 * k + 2] * sgn(T, n2 + ia);
            }
            if (k == n2) {
                c[i] = y[2 * k - 1] * sgn(T, i + 1) + y[2 * k] * sgn(T, i);
                c[big - ia] = y[2 * k - 1] * sgn(T, i + 2) + y[2 * k] * sgn(T, ia);
            }
        }

        /// The even DCT-IV, by the complex FFT of half its length
        /// (<https://www.appletonaudio.com/blog/2013/derivation-of-fast-dct-4-algorithm-based-on-dft/>,
        /// as PocketFFT).
        fn even4(self: *const Self, c: []T, rest: []T, f: T) void {
            const big = self.n;
            const n2 = big / 2;
            const cs = asComplex(T, rest);
            const y = cs[0..n2];
            const c2 = self.c2;
            for (y, c2, 0..) |*v, w, i| v.* = mul(T, .{ .re = c[2 * i], .im = c[big - 1 - 2 * i] }, w);
            self.inner.complex.exec(y, cs[n2..], f, true);
            for (0..n2) |i| {
                const ic = n2 - 1 - i;
                c[2 * i] = 2 * (y[i].re * c2[i].re - y[i].im * c2[i].im);
                c[2 * i + 1] = -2 * (y[ic].im * c2[ic].re + y[ic].re * c2[ic].im);
            }
        }
    };
}

/// `MPINPLACE`: `a := a - b`, `b := a + b`.
inline fn mpInPlace(a: anytype, b: @TypeOf(a)) void {
    const t = a.*;
    a.* -= b.*;
    b.* = t + b.*;
}

/// `cmplx::operator*=`, in its order.
inline fn mul(comptime T: type, a: C(T), b: C(T)) C(T) {
    return .{ .re = a.re * b.re - a.im * b.im, .im = a.re * b.im + a.im * b.re };
}

/// `T_dcst4`'s `SGN`: `-√2` when bit 1 of `i` is set, else `√2`.
inline fn sgn(comptime T: type, i: usize) T {
    const s2: T = sqrt2;
    return if (i & 2 != 0) -s2 else s2;
}

/// Reals' storage as complex values, as `RealPlan` lends Bluestein's.
fn asComplex(comptime T: type, rest: []T) []C(T) {
    const bytes = std.mem.sliceAsBytes(rest);
    return @alignCast(std.mem.bytesAsSlice(C(T), bytes[0 .. bytes.len / @sizeOf(C(T)) * @sizeOf(C(T))]));
}
