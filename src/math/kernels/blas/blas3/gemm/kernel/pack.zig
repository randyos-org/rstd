//! The gemm kernels' packing copies: panels of an unroll's width, their remainders, and the side kernels.

const std = @import("std");
const types = @import("../../../types.zig");
const kroot = @import("../kernel.zig");

const BlasLong = types.BlasLong;
const panelWidth = kroot.panelWidth;

/// Packs `n` lines of depth `m` of elements `E` in panels of `u`. `across`:
/// line `j` is element `j` of each of `m` rows `lda` apart (`tcopy`); else it
/// is column `j`, `lda` apart (`ncopy`).
pub fn pack(comptime E: type, comptime u: usize, comptime across: bool, m: BlasLong, n: BlasLong, a: [*]const E, lda: BlasLong, b: [*]E) void {
    const depth: usize = @intCast(m);
    const lines: usize = @intCast(n);
    const ld: usize = @intCast(lda);
    var c0: usize = 0;
    while (c0 + u <= lines) : (c0 += u) panel(E, u, across, depth, a, ld, c0, b + depth * c0);
    if (c0 < lines) packRest(E, across, depth, a, ld, c0, lines, b);
}

/// The panels past the whole ones, widths as values: they are few, and one
/// function for every width keeps debug builds small.
fn packRest(comptime E: type, comptime across: bool, depth: usize, a: [*]const E, ld: usize, from: usize, lines: usize, b: [*]E) void {
    var c0 = from;
    while (c0 < lines) {
        const w = std.math.floorPowerOfTwo(usize, lines - c0);
        const out = b + depth * c0;
        for (0..depth) |l| for (0..w) |c| {
            out[l * w + c] = if (across) a[l * ld + c0 + c] else a[(c0 + c) * ld + l];
        };
        c0 += w;
    }
}

inline fn panel(comptime E: type, comptime w: usize, comptime across: bool, depth: usize, a: [*]const E, ld: usize, c0: usize, dst: [*]E) void {
    if (across) {
        for (0..depth) |l| dst[l * w ..][0..w].* = a[l * ld + c0 ..][0..w].*;
    } else {
        var cols: [w][*]const E = undefined;
        inline for (0..w) |c| cols[c] = a + (c0 + c) * ld;
        for (0..depth) |l| {
            inline for (0..w) |c| dst[l * w + c] = cols[c][l];
        }
    }
}

/// Packs `lines` lines of `depth` as `pack` does, in panels of `u`, element
/// `(j, l)` being `src.at(j, l)`; where that is an optional, null leaves it
/// unwritten. For the other routines' copies, which transform what they pack;
/// the widths are run-time values, so each copy is one small function.
pub fn packWith(comptime E: type, u: usize, depth: usize, lines: usize, dst: [*]E, src: anytype) void {
    var c0: usize = 0;
    while (c0 < lines) {
        const w = panelWidth(lines - c0, u);
        const out = dst + depth * c0;
        for (0..depth) |l| for (0..w) |c| {
            const v = src.at(c0 + c, l);
            if (@TypeOf(v) == ?E) {
                if (v) |x| out[l * w + c] = x;
            } else out[l * w + c] = v;
        };
        c0 += w;
    }
}

/// Packs as `packWith` the triangle the trmm and trsm copies take: line `j`,
/// depth `i` is `a[j + i ld]` (`across`) or `a[i + j ld]`, the diagonal where
/// `i == j + off`, kept on `upper`'s side of it and 1 there for `unit`; the
/// other side is left unwritten. Each depth's kept lines are one run, copied
/// without a test per element.
pub fn packTriangle(comptime E: type, comptime upper: bool, comptime across: bool, comptime unit: bool, u: usize, depth: usize, lines: usize, dst: [*]E, a: [*]const E, ld: usize, off: isize) void {
    const one: E = if (@typeInfo(E) == .array) .{ 1, 0 } else 1;
    // the kept lines are `j <= i - off`, else `j >= i - off`
    const below = across == upper;
    var c0: usize = 0;
    while (c0 < lines) {
        const w = panelWidth(lines - c0, u);
        const wi: isize = @intCast(w);
        const out = dst + depth * c0;
        for (0..depth) |i| {
            // the panel's diagonal line
            const t = @as(isize, @intCast(i)) - off - @as(isize, @intCast(c0));
            const lo: usize = if (below) 0 else @intCast(std.math.clamp(t, 0, wi));
            const hi: usize = if (below) @intCast(std.math.clamp(t + 1, 0, wi)) else w;
            const row = out + i * w;
            for (lo..hi) |c| row[c] = if (across) a[c0 + c + i * ld] else a[i + (c0 + c) * ld];
            if (unit and t >= 0 and t < wi) row[@intCast(t)] = one;
        }
        c0 += w;
    }
}

/// A copy in the reference's ABI on `packWith`: `m` deep, `n` lines into
/// `b`, from `src`; the trmm, trsm, symm and 3M copies.
pub fn copyFrom(comptime E: type, u: usize, m: BlasLong, n: BlasLong, b: anytype, src: anytype) i32 {
    packWith(E, u, @intCast(m), @intCast(n), @ptrCast(b), src);
    return 0;
}

/// `C += alpha * A * B` on packed panels through `G`'s gemm kernel
/// (`Kernels(T, geo)`), conjugating the triangle's side, `left` or right,
/// for `conj`: trmm's and trsm's kernels. `alpha` is `[re, im]`, or `[re]`.
pub fn sideKernel(comptime G: type, comptime F: type, comptime complex: bool, left: bool, conj: bool, m: usize, n: usize, k: usize, alpha: [*]const F, a: anytype, b: anytype, c: anytype, ldc: usize) void {
    const bm: BlasLong = @intCast(m);
    const bn: BlasLong = @intCast(n);
    const bk: BlasLong = @intCast(k);
    if (comptime complex) {
        const f = if (!conj) &G.gemm_kernel_n else if (left) &G.gemm_kernel_l else &G.gemm_kernel_r;
        _ = f(bm, bn, bk, alpha[0], alpha[1], @ptrCast(a), @ptrCast(b), @ptrCast(c), @intCast(ldc));
    } else {
        _ = G.gemm_kernel(bm, bn, bk, alpha[0], @ptrCast(a), @ptrCast(b), @ptrCast(c), @intCast(ldc));
    }
}
