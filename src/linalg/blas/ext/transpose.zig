//! The tiled transposed copy, `B := alpha*A^T` for `f32` and `f64`:
//! `kernel/x86_64/domatcopy_t_avx.c` (OpenBLAS b439d63b8; BSD-3, `LICENSE`)
//! as `@Vector` code for any target. The matrix is walked in 32x32 tiles
//! of 4x4 blocks transposed in registers; where a leading dimension is a
//! multiple of the block, the first rows and columns are peeled so the
//! blocks start on vector boundaries. Each element is one product (a copy
//! for `alpha` 1), so the bits are the plain loops'. No BLAS-root imports.

const util = @import("../util.zig");

/// `TILE`: the tile's side, in elements.
pub const tile = 32;
/// The block's side: four vectors of four.
const side = 4;

/// Whether `F` has the tiled kernel: the reference's `d`, and `s` alike.
pub fn serves(comptime F: type) bool {
    return F == f32 or F == f64;
}

/// `b[j*ldb + i] = alpha*a[i*lda + j]` for `i_lo <= i < i_hi`, `j_lo <= j < j_hi`.
inline fn plain(comptime F: type, comptime scaled: bool, alpha: F, a: [*]const F, lda: usize, b: [*]F, ldb: usize, i_lo: usize, i_hi: usize, j_lo: usize, j_hi: usize) void {
    for (i_lo..i_hi) |i| for (j_lo..j_hi) |j| {
        const x = a[i * lda + j];
        b[j * ldb + i] = if (scaled) alpha * x else x;
    };
}

/// `dtrans`: `m` columns of `a`, `n` long, into the rows of `b`.
inline fn blocks(comptime F: type, comptime scaled: bool, m: usize, n: usize, alpha: F, a: [*]const F, lda: usize, b: [*]F, ldb: usize) void {
    const V = @Vector(side, F);
    const va: V = @splat(alpha);
    const m4 = m - m % side;
    const n4 = n - n % side;
    var is: usize = 0;
    while (is < m4) : (is += tile) {
        const ie = @min(is + tile, m4);
        var js: usize = 0;
        while (js < n4) : (js += tile) {
            const je = @min(js + tile, n4);
            var i = is;
            while (i < ie) : (i += side) {
                const ap = a + i * lda;
                const bp = b + i;
                var j = js;
                while (j < je) : (j += side) {
                    var r: [side]V = undefined;
                    inline for (0..side) |k| {
                        const v = util.load(V, ap + k * lda + j);
                        r[k] = if (scaled) v * va else v;
                    }
                    // `unpacklo`/`unpackhi` of each pair, then the halves
                    const t0 = @shuffle(F, r[0], r[1], [side]i32{ 0, ~@as(i32, 0), 2, ~@as(i32, 2) });
                    const t1 = @shuffle(F, r[0], r[1], [side]i32{ 1, ~@as(i32, 1), 3, ~@as(i32, 3) });
                    const t2 = @shuffle(F, r[2], r[3], [side]i32{ 0, ~@as(i32, 0), 2, ~@as(i32, 2) });
                    const t3 = @shuffle(F, r[2], r[3], [side]i32{ 1, ~@as(i32, 1), 3, ~@as(i32, 3) });
                    const lo = [side]i32{ 0, 1, ~@as(i32, 0), ~@as(i32, 1) };
                    const hi = [side]i32{ 2, 3, ~@as(i32, 2), ~@as(i32, 3) };
                    util.store(V, bp + (j + 0) * ldb, @shuffle(F, t0, t2, lo));
                    util.store(V, bp + (j + 1) * ldb, @shuffle(F, t1, t3, lo));
                    util.store(V, bp + (j + 2) * ldb, @shuffle(F, t0, t2, hi));
                    util.store(V, bp + (j + 3) * ldb, @shuffle(F, t1, t3, hi));
                }
            }
        }
    }
    plain(F, scaled, alpha, a, lda, b, ldb, 0, m4, n4, n);
    plain(F, scaled, alpha, a, lda, b, ldb, m4, m, 0, n);
}

/// Elements from `p` to the next vector boundary.
fn peel(comptime F: type, p: [*]const F) usize {
    const bytes = side * @sizeOf(F);
    return ((bytes - (@intFromPtr(p) & (bytes - 1))) & (bytes - 1)) / @sizeOf(F);
}

/// `dtrans_aligned`: with a leading dimension a multiple of the block
/// every column is aligned alike, so the first few rows and columns go
/// one by one and the blocks start aligned.
fn aligned(comptime F: type, comptime scaled: bool, m: usize, n: usize, alpha: F, a: [*]const F, lda: usize, b: [*]F, ldb: usize) void {
    const pi = if (ldb % side == 0) @min(peel(F, b), m) else 0;
    const pj = if (lda % side == 0) @min(peel(F, a), n) else 0;
    plain(F, scaled, alpha, a, lda, b, ldb, 0, pi, 0, n);
    plain(F, scaled, alpha, a, lda, b, ldb, pi, m, 0, pj);
    blocks(F, scaled, m - pi, n - pj, alpha, a + pi * lda + pj, lda, b + pj * ldb + pi, ldb);
}

/// `B := alpha*A^T`, column-major: `A` `rows x cols`, `B` `cols x rows`,
/// both non-empty, `alpha` not zero. `alpha` 1 copies the bits.
pub fn run(comptime F: type, rows: usize, cols: usize, alpha: F, a: [*]const F, lda: usize, b: [*]F, ldb: usize) void {
    if (alpha == 1) return aligned(F, false, cols, rows, alpha, a, lda, b, ldb);
    aligned(F, true, cols, rows, alpha, a, lda, b, ldb);
}
