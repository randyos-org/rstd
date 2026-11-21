//! pffft's f32 FFT (marton78/pffft efc1dd2d, `src/pffft_priv_impl.h` as
//! `src/pffft.c` builds it; UCAR, LICENSE): its setup, transform and
//! convolution, the kernels in `cfft.zig`, `rfftf.zig`, `rfftb.zig` and
//! `finalize.zig` over `simd.zig`'s macro layer. Sizes are pffft's: real
//! `n` a multiple of 32, complex a multiple of 16, `n / 32` or `n / 16`
//! of the form `2^a 3^b 5^c`, `n` at most `2^26`. Transforms are unscaled.
//! With AVX, sizes 8 lanes take (`validFor`) run 8 wide (`finalize8.zig`,
//! N35), so the unordered order depends on the size and the core.

const std = @import("std");
const sysinfo = @import("builtin");
const simd = @import("simd.zig");
const cfft = @import("cfft.zig");
const rfftf = @import("rfftf.zig");
const rfftb = @import("rfftb.zig");
const finalize = @import("finalize.zig");
const finalize8 = @import("finalize8.zig");

const Allocator = std.mem.Allocator;

/// Whether the 8-lane steps build: with AVX, and in tests on any target.
const eight = simd.wide or sysinfo.is_test;

pub const Kind = enum { real, complex };

/// `Setup` for f32, else a placeholder no plan builds (f64 SIMD is later).
pub fn Simd(comptime T: type) type {
    return if (T == f32) Setup else struct {
        pub fn deinit(_: *@This(), _: Allocator) void {}
    };
}

comptime {
    // `Plan(f32)` hands its values to pffft as interleaved floats
    const C = std.math.Complex(f32);
    std.debug.assert(@sizeOf(C) == 8 and @offsetOf(C, "re") == 0 and @offsetOf(C, "im") == 4);
}

/// `IFAC_MAX_SIZE`: the factors, plus `n` and their count.
const ifac_max = 25;

/// `pffft_min_fft_size`.
pub fn minSize(kind: Kind) usize {
    return minFor(kind, 4);
}

fn minFor(kind: Kind, lanes: usize) usize {
    return switch (kind) {
        .real => 2 * lanes * lanes,
        .complex => lanes * lanes,
    };
}

/// `pffft_is_valid_size`, and `pffft_new_setup`'s bound.
pub fn validSize(n: usize, kind: Kind) bool {
    return validFor(n, kind, 4);
}

/// `validSize` for `lanes`-wide steps.
fn validFor(n: usize, kind: Kind, lanes: usize) bool {
    if (n == 0 or n > 1 << 26) return false;
    const n_min = minFor(kind, lanes);
    var r = n;
    while (r >= 5 * n_min and r % 5 == 0) r /= 5;
    while (r >= 3 * n_min and r % 3 == 0) r /= 3;
    while (r >= 2 * n_min and r % 2 == 0) r /= 2;
    return r == n_min;
}

/// `PFFFT_Setup`: immutable once built, so one serves every thread.
pub const Setup = struct {
    n: usize,
    /// the vectors' width: 8 with AVX at the sizes it takes, else 4
    lanes: usize,
    /// complex vectors: `n / 2 lanes` real, `n / lanes` complex
    ncvec: usize,
    ifac: [ifac_max]usize,
    kind: Kind,
    /// `e` (`2 (lanes - 1) ncvec` floats), then the passes' twiddles (`2 ncvec`)
    data: []align(32) f32,

    /// `pffft_new_setup`; `n` must pass `validSize`.
    pub fn init(allocator: Allocator, n: usize, kind: Kind) Allocator.Error!Setup {
        std.debug.assert(validSize(n, kind));
        return initLanes(allocator, n, kind, if (simd.wide and validFor(n, kind, 8)) 8 else 4);
    }

    /// `init` at a given width, 4 or 8 (which every target runs, if slowly), for tests.
    pub fn initLanes(allocator: Allocator, n: usize, kind: Kind, lanes: usize) Allocator.Error!Setup {
        std.debug.assert(validFor(n, kind, lanes));
        const ncvec = (if (kind == .real) n / 2 else n) / lanes;
        const data = try allocator.alignedAlloc(f32, .@"32", 2 * lanes * ncvec);
        @memset(data, 0); // the real twiddles leave gaps
        var s: Setup = .{ .n = n, .lanes = lanes, .ncvec = ncvec, .ifac = @splat(0), .kind = kind, .data = data };
        const rows = lanes - 1;
        const ev = data[0 .. 2 * rows * ncvec];
        const pi = simd.lit(simd.m_pi);
        const nf: f32 = @floatFromInt(n);
        for (0..ncvec) |k| {
            const i = k / lanes;
            const j = k % lanes;
            const kf: f32 = @floatFromInt(k);
            for (0..rows) |m| {
                const mf: f32 = @floatFromInt(m + 1);
                const a = -2 * pi * mf * kf / nf;
                ev[(2 * (i * rows + m) + 0) * lanes + j] = @cos(a);
                ev[(2 * (i * rows + m) + 1) * lanes + j] = @sin(a);
            }
        }
        const wa = data[2 * rows * ncvec ..];
        switch (kind) {
            .real => rffti1(n / lanes, wa, &s.ifac),
            .complex => cffti1(n / lanes, wa, &s.ifac),
        }
        // `validSize` leaves only the factors 2, 3, 4, 5 and 8
        var m: usize = 1;
        for (s.ifac[2..][0..s.ifac[1]]) |f| m *= f;
        std.debug.assert(m == n / lanes);
        return s;
    }

    pub fn deinit(self: *Setup, allocator: Allocator) void {
        allocator.free(self.data);
        self.* = undefined;
    }

    fn e(self: *const Setup, comptime lanes: usize) simd.For(lanes).ConstPtr {
        return @ptrCast(self.data.ptr);
    }

    fn twiddle(self: *const Setup) [*]const f32 {
        return self.data[2 * (self.lanes - 1) * self.ncvec ..].ptr;
    }

    /// `pffft_transform_internal`: `n` reals (complex: `2 n` floats,
    /// interleaved) from `input` to `output`, which may be the same;
    /// `scratch` holds as many and overlaps neither. `ordered` is `pffft_transform_ordered`'s
    /// natural order, else pffft's own (`pffft_transform`).
    pub fn transform(self: *const Setup, finput: [*]const f32, foutput: [*]f32, scratch: [*]f32, forward: bool, ordered: bool) void {
        if (eight and self.lanes == 8)
            self.transformFor(8, finput, foutput, scratch, forward, ordered)
        else
            self.transformFor(4, finput, foutput, scratch, forward, ordered);
    }

    fn transformFor(self: *const Setup, comptime lanes: usize, finput: [*]const f32, foutput: [*]f32, scratch: [*]f32, forward: bool, ordered: bool) void {
        const S = simd.For(lanes);
        const C = cfft.For(lanes);
        const RF = rfftf.For(lanes);
        const RB = rfftb.For(lanes);
        const F = if (lanes == 8) finalize8 else finalize;
        const ConstPtr = S.ConstPtr;
        const Ptr = S.Ptr;
        const ncvec = self.ncvec;
        const nf_odd = self.ifac[1] & 1 == 1;
        var vinput: ConstPtr = @ptrCast(finput);
        const voutput: Ptr = @ptrCast(foutput);
        const buff = [2]Ptr{ voutput, @ptrCast(scratch) };
        var ib: u1 = @intFromBool(nf_odd != ordered);
        const e_ = self.e(lanes);
        const ifac = self.ifac[0..];

        if (forward) {
            ib = ~ib;
            if (self.kind == .real) {
                ib = @intFromBool(RF.rfftf1(ncvec * 2, vinput, buff[ib], buff[~ib], self.twiddle(), ifac) != buff[0]);
                F.realFinalize(ncvec, buff[ib], buff[~ib], e_);
            } else {
                const tmp = buff[ib];
                for (0..ncvec) |k| {
                    tmp[k * 2 ..][0..2].* = S.uninterleave2(vinput[k * 2], vinput[k * 2 + 1]);
                }
                ib = @intFromBool(C.cfftf1(ncvec, buff[ib], buff[~ib], buff[ib], self.twiddle(), ifac, -1) != buff[0]);
                F.cplxFinalize(ncvec, buff[ib], buff[~ib], e_);
            }
            if (ordered) {
                F.zreorder(self.n, ncvec, self.kind == .real, @ptrCast(buff[~ib]), @ptrCast(buff[ib]), true);
            } else ib = ~ib;
        } else {
            if (vinput == buff[ib]) ib = ~ib; // finput == foutput
            if (ordered) {
                F.zreorder(self.n, ncvec, self.kind == .real, @ptrCast(vinput), @ptrCast(buff[ib]), false);
                vinput = buff[ib];
                ib = ~ib;
            }
            if (self.kind == .real) {
                F.realPreprocess(ncvec, vinput, buff[ib], e_);
                ib = @intFromBool(RB.rfftb1(ncvec * 2, buff[ib], buff[0], buff[1], self.twiddle(), ifac) != buff[0]);
            } else {
                F.cplxPreprocess(ncvec, vinput, buff[ib], e_);
                ib = @intFromBool(C.cfftf1(ncvec, buff[ib], buff[0], buff[1], self.twiddle(), ifac, 1) != buff[0]);
                for (0..ncvec) |k| {
                    buff[ib][k * 2 ..][0..2].* = S.interleave2(buff[ib][k * 2], buff[ib][k * 2 + 1]);
                }
            }
        }

        if (buff[ib] != voutput) {
            // only when finput == foutput
            std.debug.assert(@intFromPtr(finput) == @intFromPtr(foutput));
            for (0..ncvec) |k| {
                const a = buff[ib][2 * k];
                const b = buff[ib][2 * k + 1];
                voutput[2 * k] = a;
                voutput[2 * k + 1] = b;
            }
            ib = ~ib;
        }
        std.debug.assert(buff[ib] == voutput);
    }

    /// `pffft_zconvolve_accumulate`: `ab += a b scaling` on spectra in
    /// pffft's order. The three may alias.
    pub fn convolveAccumulate(self: *const Setup, a: [*]const f32, b: [*]const f32, ab: [*]f32, scaling: f32) void {
        if (eight and self.lanes == 8)
            self.convolveFor(8, a, b, ab, scaling)
        else
            self.convolveFor(4, a, b, ab, scaling);
    }

    fn convolveFor(self: *const Setup, comptime lanes: usize, a: [*]const f32, b: [*]const f32, ab: [*]f32, scaling: f32) void {
        const S = simd.For(lanes);
        const ncvec = self.ncvec;
        const va: S.ConstPtr = @ptrCast(a);
        const vb: S.ConstPtr = @ptrCast(b);
        const vab: S.Ptr = @ptrCast(ab);
        const vscal = S.ld(scaling);

        // a real spectrum's `X[0]` and `X[n/2]`, both real, share lane 0
        const ar0 = a[0];
        const ai0 = a[lanes];
        const br0 = b[0];
        const bi0 = b[lanes];
        const abr0 = ab[0];
        const abi0 = ab[lanes];

        var i: usize = 0;
        while (i < ncvec) : (i += 2) {
            var ar = va[2 * i + 0];
            var ai = va[2 * i + 1];
            var br = vb[2 * i + 0];
            var bi = vb[2 * i + 1];
            S.cplxMul(&ar, &ai, br, bi);
            vab[2 * i + 0] = S.madd(ar, vscal, vab[2 * i + 0]);
            vab[2 * i + 1] = S.madd(ai, vscal, vab[2 * i + 1]);
            ar = va[2 * i + 2];
            ai = va[2 * i + 3];
            br = vb[2 * i + 2];
            bi = vb[2 * i + 3];
            S.cplxMul(&ar, &ai, br, bi);
            vab[2 * i + 2] = S.madd(ar, vscal, vab[2 * i + 2]);
            vab[2 * i + 3] = S.madd(ai, vscal, vab[2 * i + 3]);
        }
        if (self.kind == .real) {
            ab[0] = abr0 + ar0 * br0 * scaling;
            ab[lanes] = abi0 + ai0 * bi0 * scaling;
        }
    }
};

/// `decompose`: `n`'s factors from `ntryh` into `ifac[2..]`, a 2 first.
fn decompose(n: usize, ifac: *[ifac_max]usize, ntryh: []const usize) usize {
    var nl = n;
    var nf: usize = 0;
    for (ntryh) |ntry| {
        while (nl != 1) {
            const nq = nl / ntry;
            const nr = nl - ntry * nq;
            if (nr != 0) break;
            std.debug.assert(2 + nf < ifac_max);
            ifac[2 + nf] = ntry;
            nf += 1;
            nl = nq;
            if (ntry == 2 and nf != 1) {
                var i: usize = 2;
                while (i <= nf) : (i += 1) {
                    const ib = nf - i + 2;
                    ifac[ib + 1] = ifac[ib];
                }
                ifac[2] = 2;
            }
        }
    }
    ifac[0] = n;
    ifac[1] = nf;
    return nf;
}

/// `rffti1_ps`: the real passes' factors and twiddles for `n` vectors.
fn rffti1(n: usize, wa: []f32, ifac: *[ifac_max]usize) void {
    const nf = decompose(n, ifac, &.{ 4, 2, 3, 5 });
    const argh = (2 * simd.lit(simd.m_pi)) / @as(f32, @floatFromInt(n));
    var is: usize = 0;
    var l1: usize = 1;
    for (1..nf) |k1| {
        const ip = ifac[k1 + 1];
        var ld: usize = 0;
        const l2 = l1 * ip;
        const ido = n / l2;
        for (1..ip) |_| {
            var i = is;
            var fi: usize = 0;
            ld += l1;
            const argld = @as(f32, @floatFromInt(ld)) * argh;
            var ii: usize = 3;
            while (ii <= ido) : (ii += 2) {
                i += 2;
                fi += 1;
                const arg = @as(f32, @floatFromInt(fi)) * argld;
                wa[i - 2] = @cos(arg);
                wa[i - 1] = @sin(arg);
            }
            is += ido;
        }
        l1 = l2;
    }
}

/// `cffti1_ps`: the complex passes' factors and twiddles for `n` vectors;
/// radix 8 only for `n == 8` (pffft's measurements).
fn cffti1(n: usize, wa: []f32, ifac: *[ifac_max]usize) void {
    const ntryh: []const usize = if (n == 8) &.{ 5, 3, 8, 4, 2 } else &.{ 5, 3, 4, 2 };
    const nf = decompose(n, ifac, ntryh);
    const argh = (2 * simd.lit(simd.m_pi)) / @as(f32, @floatFromInt(n));
    var i: usize = 1;
    var l1: usize = 1;
    for (1..nf + 1) |k1| {
        const ip = ifac[k1 + 1];
        var ld: usize = 0;
        const l2 = l1 * ip;
        const ido = n / l2;
        const idot = ido + ido + 2;
        for (1..ip) |_| {
            const ia = i;
            var fi: usize = 0;
            wa[i - 1] = 1;
            wa[i] = 0;
            ld += l1;
            const argld = @as(f32, @floatFromInt(ld)) * argh;
            var ii: usize = 4;
            while (ii <= idot) : (ii += 2) {
                i += 2;
                fi += 1;
                const arg = @as(f32, @floatFromInt(fi)) * argld;
                wa[i - 1] = @cos(arg);
                wa[i] = @sin(arg);
            }
            if (ip > 16) {
                wa[ia - 1] = wa[i - 1];
                wa[ia] = wa[i];
            }
        }
        l1 = l2;
    }
}

test "validSize: pffft's sizes" {
    const testing = std.testing;
    for ([_]usize{ 32, 64, 96, 160, 480, 1024, 1 << 26 }) |n| try testing.expect(validSize(n, .real));
    for ([_]usize{ 0, 16, 48, 224, 32 * 7, 1 << 27 }) |n| try testing.expect(!validSize(n, .real));
    for ([_]usize{ 16, 48, 80, 144, 4800 }) |n| try testing.expect(validSize(n, .complex));
    for ([_]usize{ 8, 24, 16 * 7, 16 * 11 }) |n| try testing.expect(!validSize(n, .complex));
}
