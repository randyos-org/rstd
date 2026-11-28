//! One column of the bidiagonalization's fused pass (`bidiag.zig`), from
//! faer v0.24.4 (<https://codeberg.org/sarah-quinones/faer>, tag
//! `faer-v0.24.4`, MIT): `bidiag_fused_op_simd` in
//! `faer/src/linalg/svd/bidiag.rs`, `a += up yj + z vj` and then `uᴴ a` in
//! four vector partial sums (its `simd_iter!` over `[indices; 4]`), by
//! `@Vector`, multiply-adds fused where the build fuses them (`fma.zig`).
//! The remainder past whole vectors runs scalar, where faer masks it.

const std = @import("std");
const types = @import("../blas/types.zig");
const fma = @import("../blas/fma.zig");
const vector = @import("../blas/vector/util.zig");
const ops = @import("ops.zig");

const Float = types.Float;
const isComplex = types.isComplex;

/// `a := a + up yj + z vj` over `m` entries, then returns `uᴴ a`.
pub fn column(comptime T: type, m: usize, a: [*]T, u: [*]const T, up: [*]const T, z: [*]const T, yj: T, vj: T) T {
    const F = Float(T);
    const width = vector.lanes(F);
    if (comptime width < 2) return scalar(T, m, a, u, up, z, yj, vj);
    return if (comptime isComplex(T)) complexColumn(T, width, m, a, u, up, z, yj, vj) else realColumn(F, width, m, a, u, up, z, yj, vj);
}

/// The column without vectors, four partial sums by index.
fn scalar(comptime T: type, m: usize, a: [*]T, u: [*]const T, up: [*]const T, z: [*]const T, yj: T, vj: T) T {
    var acc: [4]T = @splat(ops.zero(T));
    for (0..m) |i| {
        var x = ops.add(T, a[i], ops.mul(T, up[i], yj));
        x = ops.add(T, x, ops.mul(T, z[i], vj));
        a[i] = x;
        acc[i % 4] = ops.add(T, acc[i % 4], ops.mul(T, ops.conj(T, u[i]), x));
    }
    return ops.add(T, ops.add(T, acc[0], acc[1]), ops.add(T, acc[2], acc[3]));
}

fn realColumn(comptime F: type, comptime w: usize, m: usize, a: [*]F, u: [*]const F, up: [*]const F, z: [*]const F, yj: F, vj: F) F {
    const V = @Vector(w, F);
    const ys: V = @splat(yj);
    const vs: V = @splat(vj);
    var acc: [4]V = @splat(@splat(0));
    var i: usize = 0;
    // four vectors a step, then the whole vectors left, one accumulator each
    while (i + 4 * w <= m) : (i += 4 * w) {
        inline for (0..4) |k| step(F, w, i + k * w, a, u, up, z, ys, vs, &acc[k]);
    }
    var k: usize = 0;
    while (i + w <= m) : ({
        i += w;
        k += 1;
    }) step(F, w, i, a, u, up, z, ys, vs, &acc[k]);
    var tail: F = 0;
    while (i < m) : (i += 1) {
        const x = fma.madd(true, F, z[i], vj, fma.madd(true, F, up[i], yj, a[i]));
        a[i] = x;
        tail = fma.madd(true, F, u[i], x, tail);
    }
    return @reduce(.Add, (acc[0] + acc[1]) + (acc[2] + acc[3])) + tail;
}

/// One vector of `realColumn`: `a += up ys + z vs`, `acc += u a`.
inline fn step(comptime F: type, comptime w: usize, o: usize, a: [*]F, u: [*]const F, up: [*]const F, z: [*]const F, ys: @Vector(w, F), vs: @Vector(w, F), acc: *@Vector(w, F)) void {
    const V = @Vector(w, F);
    const pa: *[w]F = a[o..][0..w];
    var x: V = pa.*;
    x = fma.madd(true, V, up[o..][0..w].*, ys, x);
    x = fma.madd(true, V, z[o..][0..w].*, vs, x);
    pa.* = x;
    acc.* = fma.madd(true, V, u[o..][0..w].*, x, acc.*);
}

/// The complex column on interleaved parts, `w / 2` elements a vector:
/// `p · s = P sr + swap(P) [-si, si, ...]`, and `conj(u) x` summed as `U X`
/// (both lanes into the real part) and `U swap(X)` (even lanes less odd
/// into the imaginary).
fn complexColumn(comptime T: type, comptime w: usize, m: usize, a: [*]T, u: [*]const T, up: [*]const T, z: [*]const T, yj: T, vj: T) T {
    const F = Float(T);
    const V = @Vector(w, F);
    const half = w / 2;
    const fa: [*]F = @ptrCast(a);
    const fu: [*]const F = @ptrCast(u);
    const fup: [*]const F = @ptrCast(up);
    const fz: [*]const F = @ptrCast(z);
    const yr: V = @splat(yj.re);
    const yi: V = alternate(F, w, yj.im);
    const vr: V = @splat(vj.re);
    const vi: V = alternate(F, w, vj.im);
    var re: [4]V = @splat(@splat(0));
    var im: [4]V = @splat(@splat(0));
    var i: usize = 0;
    while (i + 4 * half <= m) : (i += 4 * half) {
        inline for (0..4) |k| complexStep(F, w, 2 * (i + k * half), fa, fu, fup, fz, .{ yr, yi, vr, vi }, &re[k], &im[k]);
    }
    var k: usize = 0;
    while (i + half <= m) : ({
        i += half;
        k += 1;
    }) complexStep(F, w, 2 * i, fa, fu, fup, fz, .{ yr, yi, vr, vi }, &re[k], &im[k]);
    var tail = ops.zero(T);
    while (i < m) : (i += 1) {
        var x = ops.add(T, a[i], ops.mul(T, up[i], yj));
        x = ops.add(T, x, ops.mul(T, z[i], vj));
        a[i] = x;
        tail = ops.add(T, tail, ops.mul(T, ops.conj(T, u[i]), x));
    }
    const sign = alternate(F, w, -1);
    const sr = @reduce(.Add, (re[0] + re[1]) + (re[2] + re[3]));
    const si = @reduce(.Add, ((im[0] + im[1]) + (im[2] + im[3])) * sign);
    return ops.add(T, .init(sr, si), tail);
}

/// `[-v, v, -v, v, ...]`.
fn alternate(comptime F: type, comptime w: usize, v: F) @Vector(w, F) {
    var out: [w]F = undefined;
    for (&out, 0..) |*e, l| e.* = if (l % 2 == 0) -v else v;
    return out;
}

/// The pairs of `x` swapped: `[x1, x0, x3, x2, ...]`.
inline fn swapPairs(comptime F: type, comptime w: usize, x: @Vector(w, F)) @Vector(w, F) {
    const mask = comptime blk: {
        var s: [w]i32 = undefined;
        for (&s, 0..) |*e, l| e.* = @intCast(l ^ 1);
        break :blk s;
    };
    return @shuffle(F, x, undefined, mask);
}

/// One vector of `complexColumn` at float offset `o`; `s` is `yr`, `yi`,
/// `vr`, `vi`.
inline fn complexStep(comptime F: type, comptime w: usize, o: usize, a: [*]F, u: [*]const F, up: [*]const F, z: [*]const F, s: [4]@Vector(w, F), re: *@Vector(w, F), im: *@Vector(w, F)) void {
    const V = @Vector(w, F);
    const pa: *[w]F = a[o..][0..w];
    const vu: V = u[o..][0..w].*;
    const vup: V = up[o..][0..w].*;
    const vz: V = z[o..][0..w].*;
    var x: V = pa.*;
    x = fma.madd(true, V, vup, s[0], x);
    x = fma.madd(true, V, swapPairs(F, w, vup), s[1], x);
    x = fma.madd(true, V, vz, s[2], x);
    x = fma.madd(true, V, swapPairs(F, w, vz), s[3], x);
    pa.* = x;
    re.* = fma.madd(true, V, vu, x, re.*);
    im.* = fma.madd(true, V, vu, swapPairs(F, w, x), im.*);
}
