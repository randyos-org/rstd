//! Elementwise (entrywise) binary operations on same-shaped matrices - the
//! REPL's `+`, `-`, and "dot" operators `.*`, `./`, `.\` (as opposed to
//! `generators.kron`'s Kronecker-block versions `.*.`, `./.`, `.\.`, whose
//! `KronOp` enum this reuses - same three cases, just applied entrywise
//! instead of block-expanded). Kept as its own file rather than folded
//! into `CxMat.zig` to keep files under this project's line-count
//! guideline.

const blas1 = @import("blas1.zig");
const CxMat = @import("CxMat.zig");
const generators = @import("generators.zig");
const Cx = blas1.Cx;

/// ADDITION (`X + Y`): `out(i,j) = a(i,j) + b(i,j)`. `a`, `b`, and `out`
/// must have the same shape; `out` may alias `a` or `b`.
pub fn add(a: CxMat, b: CxMat, out: CxMat) void {
    var j: usize = 1;
    while (j <= a.n) : (j += 1) {
        var i: usize = 1;
        while (i <= a.m) : (i += 1) {
            const av = a.get(i, j);
            const bv = b.get(i, j);
            out.set(i, j, Cx.init(av.re + bv.re, av.im + bv.im));
        }
    }
}

/// SUBTRACTION (`X - Y`): `out(i,j) = a(i,j) - b(i,j)`. Same shape and
/// aliasing rules as `add`.
pub fn sub(a: CxMat, b: CxMat, out: CxMat) void {
    var j: usize = 1;
    while (j <= a.n) : (j += 1) {
        var i: usize = 1;
        while (i <= a.m) : (i += 1) {
            const av = a.get(i, j);
            const bv = b.get(i, j);
            out.set(i, j, Cx.init(av.re - bv.re, av.im - bv.im));
        }
    }
}

fn dotOp(a: CxMat, b: CxMat, out: CxMat, op: generators.KronOp) void {
    var j: usize = 1;
    while (j <= a.n) : (j += 1) {
        var i: usize = 1;
        while (i <= a.m) : (i += 1) {
            const av = a.get(i, j);
            const bv = b.get(i, j);
            const v = switch (op) {
                .mul => blas1.wmul(av, bv),
                .div => blas1.wdiv(av, bv),
                .rdiv => blas1.wdiv(bv, av),
            };
            out.set(i, j, v);
        }
    }
}

/// `.*`: entrywise multiplication (Hadamard product), `out(i,j) =
/// a(i,j)*b(i,j)`. Same shape and aliasing rules as `add`.
pub fn dotMul(a: CxMat, b: CxMat, out: CxMat) void {
    dotOp(a, b, out, .mul);
}

/// `./`: entrywise division, `out(i,j) = a(i,j)/b(i,j)`.
pub fn dotDiv(a: CxMat, b: CxMat, out: CxMat) void {
    dotOp(a, b, out, .div);
}

/// `.\`: entrywise left division, `out(i,j) = b(i,j)/a(i,j)` (REPL docs:
/// "`A .\ B` has elements `a(i,j)\b(i,j)`", and `x\y == y/x` for scalars).
pub fn dotRdiv(a: CxMat, b: CxMat, out: CxMat) void {
    dotOp(a, b, out, .rdiv);
}

const std = @import("std");

const zigstd = @import("../../zigstd.zig");

test "add and sub match hand computation" {
    var ar = [_]f64{ 1.0, 2.0, 3.0, 4.0 }; // column-major 2x2
    var ai: [4]f64 = @splat(0.0);
    var br = [_]f64{ 5.0, 6.0, 7.0, 8.0 };
    var bi: [4]f64 = @splat(0.0);
    var out_r: [4]f64 = undefined;
    var out_i: [4]f64 = undefined;
    const a = CxMat.initPacked(&ar, &ai, 2, 2);
    const b = CxMat.initPacked(&br, &bi, 2, 2);
    const out = CxMat.initPacked(&out_r, &out_i, 2, 2);

    add(a, b, out);
    try zigstd.testExpectEqualSlices(f64, &[_]f64{ 6.0, 8.0, 10.0, 12.0 }, &out_r);

    sub(a, b, out);
    try zigstd.testExpectEqualSlices(f64, &[_]f64{ -4.0, -4.0, -4.0, -4.0 }, &out_r);
}

test "dotMul, dotDiv and dotRdiv match hand computation" {
    var ar = [_]f64{ 2.0, 3.0, 4.0, 5.0 }; // column-major 2x2
    var ai: [4]f64 = @splat(0.0);
    var br = [_]f64{ 10.0, 12.0, 8.0, 20.0 };
    var bi: [4]f64 = @splat(0.0);
    var out_r: [4]f64 = undefined;
    var out_i: [4]f64 = undefined;
    const a = CxMat.initPacked(&ar, &ai, 2, 2);
    const b = CxMat.initPacked(&br, &bi, 2, 2);
    const out = CxMat.initPacked(&out_r, &out_i, 2, 2);

    dotMul(a, b, out);
    try zigstd.testExpectEqualSlices(f64, &[_]f64{ 20.0, 36.0, 32.0, 100.0 }, &out_r);

    dotDiv(b, a, out);
    try zigstd.testExpectEqualSlices(f64, &[_]f64{ 5.0, 4.0, 2.0, 4.0 }, &out_r);

    dotRdiv(a, b, out);
    try zigstd.testExpectEqualSlices(f64, &[_]f64{ 5.0, 4.0, 2.0, 4.0 }, &out_r);
}
