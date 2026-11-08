//! Householder reflections and Givens rotations, from faer v0.24.4
//! (<https://codeberg.org/sarah-quinones/faer>, tag `faer-v0.24.4`, MIT):
//! `faer/src/linalg/householder.rs` and `jacobi.rs`. A block of reflections is
//! `I - V T⁻¹ Vᴴ`, `T` upper triangular with `tau` on its diagonal (Joffrain
//! et al., "Accumulating Householder transformations, revisited", ACM TOMS
//! 32(2), 2006); a sequence keeps its blocks' `T`s side by side.

const std = @import("std");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = ops.Float;
const Error = ops.Error;

/// What `make` found.
pub fn Info(comptime T: type) type {
    return struct {
        /// `tau`; infinite when the tail is below the safe minimum, so the
        /// reflection is the identity.
        tau: Float(T),
        /// `1 / (head + beta)`; null where `tau` is infinite.
        head_with_beta_inv: ?T,
        /// `‖x‖₂`.
        norm: Float(T),
    };
}

/// `make_householder_imp`: the reflection taking `x = [head; tail]` to
/// `beta e₀`; `beta` goes to `head` and `v`'s essential part to `out`, from
/// `input` (out of place) or from `out` itself. Strides are in elements.
pub fn make(io: Io, comptime T: type, head: *T, n: usize, out: [*]T, out_inc: usize, input: ?[*]const T, in_inc: usize) Info(T) {
    const F = Float(T);
    const src: [*]const T = input orelse out;
    const src_inc = if (input != null) in_inc else out_inc;
    const tail_norm = blas_calls.nrm2(io, T, n, @constCast(src), src_inc);
    var head_norm = ops.abs(T, head.*);
    if (head_norm < std.math.floatMin(F)) {
        head.* = ops.zero(T);
        head_norm = 0;
    }
    if (tail_norm < std.math.floatMin(F)) return .{ .tau = std.math.inf(F), .head_with_beta_inv = null, .norm = head_norm };
    const norm = std.math.hypot(head_norm, tail_norm);
    const sign = if (head_norm != 0) ops.mulReal(T, head.*, 1 / head_norm) else ops.one(T);
    const signed_norm = ops.mulReal(T, sign, norm);
    const hwbi = ops.recip(T, ops.add(T, head.*, signed_norm));
    for (0..n) |i| out[i * out_inc] = ops.mul(T, src[i * src_inc], hwbi);
    head.* = ops.neg(T, signed_norm);
    const t = tail_norm * ops.abs(T, hwbi);
    return .{ .tau = 0.5 * (1 + t * t), .head_with_beta_inv = hwbi, .norm = norm };
}

/// `upgrade_householder_factor` for single reflections into one block: `T`'s
/// strict upper triangle from `V`'s columns (`v.rows x b`, unit lower
/// trapezoidal; what lies on and above its diagonal is not read); `T`'s
/// diagonal holds the `tau`s and is kept. Faer builds it from halves
/// recursively; here it is one product, the same matrix.
pub fn blockFactor(al: Allocator, io: Io, comptime T: type, v: Mat(T), t: Mat(T)) Error!void {
    const b = v.cols;
    if (b <= 1) return;
    const taus = try al.alloc(T, b);
    defer al.free(taus);
    for (taus, 0..) |*x, i| x.* = t.get(i, i);
    // the top, explicit: unit lower
    const w = try ops.alloc(T, al, b, b);
    defer ops.free(T, al, w);
    for (0..b) |j| for (0..b) |i| {
        w.at(i, j).* = if (i == j) ops.one(T) else if (i > j) v.get(i, j) else ops.zero(T);
    };
    try blas_calls.gramUpper(al, io, T, 1, w, 0, t);
    if (v.rows > b) try blas_calls.gramUpper(al, io, T, 1, v.sub(b, 0, v.rows - b, b), 1, t);
    for (taus, 0..) |x, i| t.at(i, i).* = x;
}

/// The block factors of a sequence whose `tau`s sit in `h`'s first row.
pub fn sequenceFactors(al: Allocator, io: Io, comptime T: type, basis: Mat(T), h: Mat(T)) Error!void {
    var j: usize = 0;
    while (j < h.cols) : (j += h.rows) {
        const b = @min(h.rows, h.cols - j);
        for (0..b) |k| h.at(k, j + k).* = h.get(0, j + k);
        try blockFactor(al, io, T, basis.sub(j, j, basis.rows - j, b), h.sub(0, j, b, b));
    }
}

/// `apply_block_householder_on_the_left_in_place_generic`: `M := H M`, or
/// with `forward` `M := Hᴴ M`, where `H = I - V T⁻¹ Vᴴ` (`v.rows x b`,
/// `t` `b x b`), all conjugated with `conj`. `m` has `v.rows` rows.
pub fn applyBlock(al: Allocator, io: Io, comptime T: type, v: Mat(T), t: Mat(T), conj: bool, forward: bool, m: Mat(T)) Error!void {
    const b = v.cols;
    const k = m.cols;
    if (b == 0 or k == 0) return;
    const tmp = try ops.alloc(T, al, b, k);
    defer ops.free(T, al, tmp);
    const top = v.sub(0, 0, b, b);
    const mtop = m.sub(0, 0, b, k);
    const rest = v.rows - b;
    const one = ops.one(T);
    // tmp := Vᴴ M (conjugated: Vᵀ M)
    tmp.copyFrom(mtop);
    try blas_calls.tri(al, io, T, false, false, .of(true, !conj), true, one, top, tmp);
    if (rest > 0) try blas_calls.gemm(al, io, T, .of(true, !conj), .n, one, v.sub(b, 0, rest, b), m.sub(b, 0, rest, k), one, tmp);
    // tmp := T⁻¹ tmp, or T⁻ᴴ tmp
    try blas_calls.tri(al, io, T, true, true, if (forward) .of(true, !conj) else .of(false, conj), false, one, t, tmp);
    // M -= V tmp
    if (rest > 0) try blas_calls.gemm(al, io, T, .of(false, conj), .n, ops.neg(T, one), v.sub(b, 0, rest, b), tmp, one, m.sub(b, 0, rest, k));
    try blas_calls.tri(al, io, T, false, false, .of(false, conj), true, one, top, tmp);
    for (0..k) |j| for (0..b) |i| {
        mtop.at(i, j).* = ops.sub(T, mtop.get(i, j), tmp.get(i, j));
    };
}

/// `apply_block_householder_sequence_on_the_left_in_place_with_conj`:
/// `M := H(0) H(1) ... M` for the reflections in `basis`'s columns (unit
/// lower trapezoidal) and their factors in `factor` (`bs x size`), all
/// conjugated with `conj`; the last block first.
pub fn applySequence(al: Allocator, io: Io, comptime T: type, basis: Mat(T), factor: Mat(T), conj: bool, m: Mat(T)) Error!void {
    const size = factor.cols;
    const bs = factor.rows;
    var j = size;
    var b = size % bs;
    if (b == 0) b = bs;
    while (j > 0) {
        const j0 = j - b;
        b = bs;
        try applyBlock(al, io, T, basis.sub(j0, j0, basis.rows - j0, j - j0), factor.sub(0, j0, j - j0, j - j0), conj, false, m.sub(j0, 0, m.rows - j0, m.cols));
        j = j0;
    }
}

/// `apply_block_householder_sequence_transpose_on_the_left_in_place_with_conj`:
/// `M := ... H(1)ᵀ H(0)ᵀ M`, conjugated with `conj` (so `Hᴴ` with it); the
/// first block first.
pub fn applySequenceTransposed(al: Allocator, io: Io, comptime T: type, basis: Mat(T), factor: Mat(T), conj: bool, m: Mat(T)) Error!void {
    const size = factor.cols;
    const bs = factor.rows;
    var j: usize = 0;
    while (j < size) {
        const b = @min(bs, size - j);
        try applyBlock(al, io, T, basis.sub(j, j, basis.rows - j, b), factor.sub(0, j, b, b), !conj, true, m.sub(j, 0, m.rows - j, m.cols));
        j += b;
    }
}

/// `apply_block_householder_on_the_right_in_place_with_conj` (faer: the
/// transposed left form on `Mᵀ`): `M := M H`, or with `adjoint`
/// `M := M Hᴴ`, `H = I - V T⁻¹ Vᴴ` conjugated with `conj`. `m` has `v.rows`
/// columns.
pub fn applyBlockRight(al: Allocator, io: Io, comptime T: type, v: Mat(T), t: Mat(T), conj: bool, adjoint: bool, m: Mat(T)) Error!void {
    const b = v.cols;
    const k = m.rows;
    if (b == 0 or k == 0) return;
    const tmp = try ops.alloc(T, al, k, b);
    defer ops.free(T, al, tmp);
    const top = v.sub(0, 0, b, b);
    const mtop = m.sub(0, 0, k, b);
    const rest = v.rows - b;
    const one = ops.one(T);
    // tmp := M V
    tmp.copyFrom(mtop);
    try blas_calls.triRight(al, io, T, false, false, .of(false, conj), true, one, top, tmp);
    if (rest > 0) try blas_calls.gemm(al, io, T, .n, .of(false, conj), one, m.sub(0, b, k, rest), v.sub(b, 0, rest, b), one, tmp);
    // tmp := tmp T⁻¹, or tmp T⁻ᴴ
    try blas_calls.triRight(al, io, T, true, true, if (adjoint) .of(true, !conj) else .of(false, conj), false, one, t, tmp);
    // M -= tmp Vᴴ
    if (rest > 0) try blas_calls.gemm(al, io, T, .n, .of(true, !conj), ops.neg(T, one), tmp, v.sub(b, 0, rest, b), one, m.sub(0, b, k, rest));
    try blas_calls.triRight(al, io, T, false, false, .of(true, !conj), true, one, top, tmp);
    for (0..b) |j| for (0..k) |i| {
        mtop.at(i, j).* = ops.sub(T, mtop.get(i, j), tmp.get(i, j));
    };
}

/// `apply_block_householder_sequence_on_the_right_in_place_with_conj`:
/// `M := M H(0) H(1) ...`, or with `adjoint` `M := M ... H(1)ᴴ H(0)ᴴ`, all
/// conjugated with `conj`.
pub fn applySequenceRight(al: Allocator, io: Io, comptime T: type, basis: Mat(T), factor: Mat(T), conj: bool, adjoint: bool, m: Mat(T)) Error!void {
    const size = factor.cols;
    const bs = factor.rows;
    if (!adjoint) {
        var j: usize = 0;
        while (j < size) {
            const b = @min(bs, size - j);
            try applyBlockRight(al, io, T, basis.sub(j, j, basis.rows - j, b), factor.sub(0, j, b, b), conj, false, m.sub(0, j, m.rows, m.cols - j));
            j += b;
        }
        return;
    }
    var j = size;
    var b = size % bs;
    if (b == 0) b = bs;
    while (j > 0) {
        const j0 = j - b;
        b = bs;
        try applyBlockRight(al, io, T, basis.sub(j0, j0, basis.rows - j0, j - j0), factor.sub(0, j0, j - j0, j - j0), conj, true, m.sub(0, j0, m.rows, m.cols - j0));
        j = j0;
    }
}

/// The block factors of the reflections in `basis`'s first `factor.cols`
/// columns, `factor`'s diagonal already holding their `tau`s (faer's
/// convention), in blocks of `factor.rows`.
pub fn sequenceFactor(al: Allocator, io: Io, comptime T: type, basis: Mat(T), factor: Mat(T)) Error!void {
    const size = factor.cols;
    const bs = factor.rows;
    var j: usize = 0;
    while (j < size) : (j += bs) {
        const b = @min(bs, size - j);
        try blockFactor(al, io, T, basis.sub(j, j, basis.rows - j, b), factor.sub(0, j, b, b));
    }
}

// ---- Givens rotations (`jacobi.rs`), real -------------------------------------------

/// `[c -s; s c]`, real.
pub fn Rotation(comptime F: type) type {
    return struct {
        c: F,
        s: F,

        const Self = @This();

        /// `make_givens`: the rotation with `c p - s q = r`, `s p + c q = 0`.
        pub fn make(p: F, q: F) Self {
            if (q == 0) return .{ .c = if (p < 0) -1 else 1, .s = 0 };
            if (p == 0) return .{ .c = 0, .s = if (q < 0) 1 else -1 };
            if (@abs(p) > @abs(q)) {
                const t = q / p;
                var u = std.math.hypot(t, @as(F, 1));
                if (p < 0) u = -u;
                const c = 1 / u;
                return .{ .c = c, .s = -t * c };
            }
            const t = p / q;
            var u = std.math.hypot(t, @as(F, 1));
            if (q < 0) u = -u;
            const s = -1 / u;
            return .{ .c = -t * s, .s = s };
        }

        /// `adjoint`, `transpose` (the same for real).
        pub fn adjoint(r: Self) Self {
            return .{ .c = r.c, .s = -r.s };
        }

        /// `apply_on_the_left_in_place`: rows `x := c x - s y`,
        /// `y := s x + c y`, `n` elements at strides `incx`, `incy`. The
        /// identity is skipped, as faer's SIMD path does.
        pub fn left(r: Self, n: usize, x: [*]F, incx: usize, y: [*]F, incy: usize) void {
            if (r.c == 1 and r.s == 0) return;
            for (0..n) |i| {
                const a = x[i * incx];
                const b = y[i * incy];
                x[i * incx] = r.c * a - r.s * b;
                y[i * incy] = r.s * a + r.c * b;
            }
        }

        /// `apply_on_the_right_in_place`: columns `[x y] := [x y] J`.
        pub fn right(r: Self, n: usize, x: [*]F, y: [*]F) void {
            r.adjoint().left(n, x, 1, y, 1);
        }
    };
}

/// `[c -conj(s); s c]`, `c` real, for a complex `T` (`jacobi.rs`'s
/// `JacobiRotation` over complex fields).
pub fn ComplexRotation(comptime T: type) type {
    return struct {
        c: Float(T),
        s: T,

        const Self = @This();
        const F = Float(T);

        /// `rotg` (`xLARTG`): the rotation with `c p - conj(s) q = r`,
        /// `s p + c q = 0`, and `r`. Faer returns `r = 1` for `q = 0`; the
        /// reference's (and the rotation's) `r` is `p` (issues.md).
        pub fn rotg(a: T, b: T) struct { Self, T } {
            if (ops.isZero(T, b)) return .{ .{ .c = 1, .s = ops.zero(T) }, a };
            const eps = std.math.floatEps(F);
            const sml = std.math.floatMin(F);
            const big = std.math.floatMax(F);
            const rtmin = @sqrt(sml / eps);
            const rtmax = 1 / rtmin;
            var c: F = undefined;
            var s: T = undefined;
            var r: T = undefined;
            if (ops.isZero(T, a)) {
                c = 0;
                const g1 = @max(@abs(b.re), @abs(b.im));
                if (g1 > rtmin and g1 < rtmax) {
                    const d = @sqrt(ops.abs2(T, b));
                    s = ops.mulReal(T, ops.conj(T, b), 1 / d);
                    r = ops.scalar(T, d);
                } else {
                    const u = @min(@max(g1, sml), big);
                    const gs = ops.mulReal(T, b, 1 / u);
                    const d = @sqrt(ops.abs2(T, gs));
                    s = ops.mulReal(T, ops.conj(T, gs), 1 / d);
                    r = ops.scalar(T, d * u);
                }
            } else {
                const f1 = @max(@abs(a.re), @abs(a.im));
                const g1 = @max(@abs(b.re), @abs(b.im));
                if (f1 > rtmin and f1 < rtmax and g1 > rtmin and g1 < rtmax) {
                    const f2 = ops.abs2(T, a);
                    const g2 = ops.abs2(T, b);
                    const h2 = f2 + g2;
                    const d = if (f2 > rtmin and h2 < rtmax) @sqrt(f2 * h2) else @sqrt(f2) * @sqrt(h2);
                    const p = 1 / d;
                    c = f2 * p;
                    s = ops.mul(T, ops.conj(T, b), ops.mulReal(T, a, p));
                    r = ops.mulReal(T, a, h2 * p);
                } else {
                    const u = @min(@max(@max(f1, g1), sml), big);
                    const uu = 1 / u;
                    const gs = ops.mulReal(T, b, uu);
                    const g2 = ops.abs2(T, gs);
                    var f2: F = undefined;
                    var h2: F = undefined;
                    var w: F = undefined;
                    var fs: T = undefined;
                    if (f1 * uu < rtmin) {
                        const v = @min(@max(f1, sml), big);
                        w = v * uu;
                        fs = ops.mulReal(T, a, 1 / v);
                        f2 = ops.abs2(T, fs);
                        h2 = (f2 * w) * w + g2;
                    } else {
                        w = 1;
                        fs = ops.mulReal(T, a, uu);
                        f2 = ops.abs2(T, fs);
                        h2 = f2 + g2;
                    }
                    const d = if (f2 > rtmin and h2 < rtmax) @sqrt(f2 * h2) else @sqrt(f2) * @sqrt(h2);
                    const p = 1 / d;
                    c = (f2 * p) * w;
                    s = ops.mul(T, ops.conj(T, gs), ops.mulReal(T, fs, p));
                    r = ops.mulReal(T, ops.mulReal(T, fs, h2 * p), u);
                }
            }
            return .{ .{ .c = c, .s = ops.neg(T, ops.conj(T, s)) }, r };
        }

        pub fn adjoint(r: Self) Self {
            return .{ .c = r.c, .s = ops.neg(T, r.s) };
        }

        /// `apply_on_the_left_in_place`: `x := c x - conj(s) y`,
        /// `y := c y + s x`.
        pub fn left(r: Self, n: usize, x: [*]T, incx: usize, y: [*]T, incy: usize) void {
            const sc = ops.conj(T, r.s);
            for (0..n) |i| {
                const a = x[i * incx];
                const b = y[i * incy];
                x[i * incx] = ops.sub(T, ops.mulReal(T, a, r.c), ops.mul(T, sc, b));
                y[i * incy] = ops.add(T, ops.mulReal(T, b, r.c), ops.mul(T, r.s, a));
            }
        }

        /// `apply_on_the_right_in_place`: columns `[x y] := [x y] J`, the
        /// transposed rotation from the left (faer's portable path).
        pub fn right(r: Self, n: usize, x: [*]T, y: [*]T) void {
            const t: Self = .{ .c = r.c, .s = ops.neg(T, ops.conj(T, r.s)) };
            t.left(n, x, 1, y, 1);
        }
    };
}

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/check.zig");
const types = @import("../blas/types.zig");

test "faer householder: H x = beta e0, and blocks of reflections match one at a time" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64, f32, types.Complex32 }) |T| {
        const F = Float(T);
        const eps: F = std.math.floatEps(F);
        const m = 23;
        const b = 5;
        const k = 4;
        // V from reflections of random columns
        var a = try ops.alloc(T, al, m, b);
        defer ops.free(T, al, a);
        check.fill(T, a.ptr[0 .. m * b], 7);
        const t = try ops.alloc(T, al, b, b);
        defer ops.free(T, al, t);
        t.fill(ops.zero(T));
        for (0..b) |j| {
            const x0 = try al.dupe(T, a.col(j)[j..m]);
            defer al.free(x0);
            const info = make(check.io, T, a.at(j, j), m - j - 1, a.col(j) + j + 1, 1, null, 0);
            t.at(j, j).* = ops.scalar(T, info.tau);
            // H x0 = beta e0
            var dot = x0[0];
            for (1..m - j) |i| dot = ops.add(T, dot, ops.mul(T, ops.conj(T, a.get(j + i, j)), x0[i]));
            const kk = ops.mulReal(T, dot, -1 / info.tau);
            try testing.expectApproxEqAbs(ops.re(T, a.get(j, j)), ops.re(T, ops.add(T, x0[0], kk)), 64 * eps);
            for (1..m - j) |i| {
                const y = ops.add(T, x0[i], ops.mul(T, a.get(j + i, j), kk));
                try testing.expect(ops.abs(T, y) <= 64 * eps * info.norm);
            }
            // later columns get this reflection, so V's columns are independent
            for (j + 1..b) |c| {
                var d = a.get(j, c);
                for (1..m - j) |i| d = ops.add(T, d, ops.mul(T, ops.conj(T, a.get(j + i, j)), a.get(j + i, c)));
                const kc = ops.mulReal(T, d, -1 / info.tau);
                a.at(j, c).* = ops.add(T, a.get(j, c), kc);
                for (1..m - j) |i| a.at(j + i, c).* = ops.add(T, a.get(j + i, c), ops.mul(T, a.get(j + i, j), kc));
            }
        }
        try blockFactor(al, check.io, T, a, t);
        inline for (.{ false, true }) |cj| inline for (.{ false, true }) |fwd| {
            const x = try ops.alloc(T, al, m, k);
            defer ops.free(T, al, x);
            check.fill(T, x.ptr[0 .. m * k], 11);
            const y = try ops.alloc(T, al, m, k);
            defer ops.free(T, al, y);
            y.copyFrom(x);
            try applyBlock(al, check.io, T, a, t, cj, fwd, y);
            // one at a time: H(0) H(1) ... (last first), or Hᴴ: first first
            for (0..b) |s| {
                const j = if (fwd) s else b - 1 - s;
                for (0..k) |c| {
                    var d = x.get(j, c);
                    for (j + 1..m) |i| {
                        const vi = if (cj) a.get(i, j) else ops.conj(T, a.get(i, j));
                        d = ops.add(T, d, ops.mul(T, vi, x.get(i, c)));
                    }
                    const kc = ops.mulReal(T, d, -1 / ops.re(T, t.get(j, j)));
                    x.at(j, c).* = ops.add(T, x.get(j, c), kc);
                    for (j + 1..m) |i| {
                        const vi = if (cj) ops.conj(T, a.get(i, j)) else a.get(i, j);
                        x.at(i, c).* = ops.add(T, x.get(i, c), ops.mul(T, vi, kc));
                    }
                }
            }
            for (0..k) |c| for (0..m) |i| {
                try testing.expect(ops.abs(T, ops.sub(T, x.get(i, c), y.get(i, c))) <= 256 * eps);
            };
        };
    }
}

test "faer householder: Givens rotations zero their second entry" {
    const R = Rotation(f64);
    for ([_][2]f64{ .{ 3, 4 }, .{ -3, 4 }, .{ 3, -4 }, .{ 0, 2 }, .{ 0, -2 }, .{ 2, 0 }, .{ -2, 0 }, .{ 1e-300, 1e300 } }) |pq| {
        const r = R.make(pq[0], pq[1]);
        try testing.expectApproxEqAbs(@as(f64, 1), r.c * r.c + r.s * r.s, 4e-16);
        try testing.expect(@abs(r.s * pq[0] + r.c * pq[1]) <= 4e-16 * @abs(std.math.hypot(pq[0], pq[1])));
    }
}

test "faer householder: rotg zeroes its second entry, and keeps r = p when q = 0" {
    const T = types.Complex64;
    const R = ComplexRotation(T);
    for ([_][2]T{ .{ .init(3, 1), .init(-2, 5) }, .{ .init(0, 0), .init(1, -1) }, .{ .init(1e-300, 0), .init(0, 1e-300) }, .{ .init(1e300, 1e300), .init(1, 0) }, .{ .init(-2, 7), .init(0, 0) } }) |pq| {
        const g = R.rotg(pq[0], pq[1]);
        var x = [_]T{pq[0]};
        var y = [_]T{pq[1]};
        g[0].left(1, &x, 1, &y, 1);
        const scale = @max(ops.abs(T, pq[0]), ops.abs(T, pq[1]));
        try testing.expect(ops.abs(T, y[0]) <= 4e-16 * scale);
        try testing.expect(ops.abs(T, ops.sub(T, x[0], g[1])) <= 4e-16 * scale);
    }
}
