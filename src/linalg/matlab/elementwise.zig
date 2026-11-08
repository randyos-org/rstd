//! The REPL's entrywise operators: `+`, `-`, `.*`, `./` and `.\`.

const scalar = @import("scalar.zig");
const CxMat = @import("CxMat.zig");
const Cx = scalar.Cx;

/// ADDITION (`X + Y`): `out(i,j) = a(i,j) + b(i,j)`. `a`, `b`, and `out`
/// must have the same shape; `out` may alias `a` or `b`.
pub fn add(a: CxMat, b: CxMat, out: CxMat) void {
    binary(a, b, out, .add);
}

/// SUBTRACTION (`X - Y`): `out(i,j) = a(i,j) - b(i,j)`. Same shape and
/// aliasing rules as `add`.
pub fn sub(a: CxMat, b: CxMat, out: CxMat) void {
    binary(a, b, out, .sub);
}

/// `out(i,j) = a(i,j) op b(i,j)`; `.rdiv` is `b(i,j) / a(i,j)`.
fn binary(a: CxMat, b: CxMat, out: CxMat, comptime op: enum { add, sub, mul, div, rdiv }) void {
    var j: usize = 1;
    while (j <= a.n) : (j += 1) {
        var i: usize = 1;
        while (i <= a.m) : (i += 1) {
            const av = a.get(i, j);
            const bv = b.get(i, j);
            const v = switch (op) {
                .add => Cx.init(av.re + bv.re, av.im + bv.im),
                .sub => Cx.init(av.re - bv.re, av.im - bv.im),
                .mul => av.mul(bv),
                .div => scalar.wdiv(av, bv),
                .rdiv => scalar.wdiv(bv, av),
            };
            out.set(i, j, v);
        }
    }
}

/// `.*`: entrywise multiplication (Hadamard product), `out(i,j) =
/// a(i,j)*b(i,j)`. Same shape and aliasing rules as `add`.
pub fn dotMul(a: CxMat, b: CxMat, out: CxMat) void {
    binary(a, b, out, .mul);
}

/// `./`: entrywise division, `out(i,j) = a(i,j)/b(i,j)`.
pub fn dotDiv(a: CxMat, b: CxMat, out: CxMat) void {
    binary(a, b, out, .div);
}

/// `.\`: entrywise left division, `out(i,j) = b(i,j)/a(i,j)` (REPL docs:
/// "`A .\ B` has elements `a(i,j)\b(i,j)`", and `x\y == y/x` for scalars).
pub fn dotRdiv(a: CxMat, b: CxMat, out: CxMat) void {
    binary(a, b, out, .rdiv);
}

// ---- tests ------------------------------------------------------------------

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
