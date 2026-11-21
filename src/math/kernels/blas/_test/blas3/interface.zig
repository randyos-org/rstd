//! `blas3/interface.zig`'s tests.

const std = @import("std");
const blas = @import("../../__root__.zig");
const types = @import("../../types.zig");
const accel = @import("../../apple_accel.zig");
const interface = @import("../../blas3/interface.zig");

const isComplex = types.isComplex;
const gemm = interface.gemm;
const gemm3m = interface.gemm3m;
const gemmt = interface.gemmt;
const symm = interface.symm;
const symm3m = interface.symm3m;
const hemm = interface.hemm;
const hemm3m = interface.hemm3m;
const syrk = interface.syrk;
const herk = interface.herk;
const syr2k = interface.syr2k;
const her2k = interface.her2k;
const trmm = interface.trmm;
const trsm = interface.trsm;

const testing = std.testing;
const check3 = @import("check.zig");
const check = @import("../check.zig");

test "interface: gemm reports the first bad argument and does nothing for empty problems" {
    const fa = testing.failing_allocator;
    var a: [16]f64 = @splat(1);
    var c: [16]f64 = @splat(std.math.nan(f64));
    const nan = std.math.nan(f64);
    try testing.expectError(error.InvalidTrans, gemm(fa, check3.io, f64, .n, .o, -1, -1, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidM, gemm(fa, check3.io, f64, .n, .n, -1, -1, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidN, gemm(fa, check3.io, f64, .n, .n, 0, -1, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidK, gemm(fa, check3.io, f64, .n, .n, 0, 0, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidLda, gemm(fa, check3.io, f64, .n, .n, 3, 2, 4, 1, &a, 2, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidLda, gemm(fa, check3.io, f64, .t, .n, 3, 2, 4, 1, &a, 3, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidLdb, gemm(fa, check3.io, f64, .n, .n, 3, 2, 4, 1, &a, 3, &a, 3, 0, &c, 1));
    try testing.expectError(error.InvalidLdb, gemm(fa, check3.io, f64, .n, .c, 3, 2, 4, 1, &a, 3, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidLdc, gemm(fa, check3.io, f64, .n, .n, 3, 2, 4, 1, &a, 3, &a, 4, 0, &c, 2));
    try testing.expectError(error.InvalidLdc, gemm(fa, check3.io, f64, .n, .n, 0, 2, 0, 1, &a, 1, &a, 1, 0, &c, 0));
    // empty: nothing read or written, nothing allocated
    try gemm(fa, check3.io, f64, .n, .n, 0, 3, 2, nan, &a, 1, &a, 2, 0, &c, 1);
    try gemm(fa, check3.io, f64, .n, .n, 3, 0, 2, nan, &a, 3, &a, 2, 0, &c, 3);
    try testing.expect(std.math.isNan(c[0]));
    // k == 0 only scales C, and allocates nothing
    try gemm(fa, check3.io, f64, .n, .n, 2, 2, 0, nan, &a, 2, &a, 1, 0, &c, 2);
    try testing.expectEqualSlices(f64, &.{ 0, 0, 0, 0 }, c[0..4]);
}

const all = check.sd_types;
const unslotted = .{ f16, f128, std.math.Complex(f16), std.math.Complex(f128) } ++ check3.f80_unslotted;

test "interface: gemm matches the wide reference through the global vtable" {
    inline for (all) |T| try check3.interfaceCases(T, .gemm);
}

test "interface: gemm over types the vtable has no slots for runs the portable kernels" {
    inline for (unslotted) |T| try check3.interfaceCases(T, .gemm);
}

test "interface: a matrix past the address space is its leading dimension's error" {
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(1);
    const big = std.math.maxInt(isize) / 2;
    try testing.expectError(error.InvalidLda, gemm(fa, check3.io, f64, .n, .n, 2, 2, 3, 1, &a, big, &a, 3, 0, &a, 2));
    try testing.expectError(error.InvalidLdb, gemm(fa, check3.io, f64, .n, .n, 2, 2, 3, 1, &a, 2, &a, big, 0, &a, 2));
    try testing.expectError(error.InvalidLdc, gemm(fa, check3.io, f64, .n, .n, 2, 2, 3, 1, &a, 2, &a, 3, 0, &a, big));
}

test "interface: gemmt matches the wide reference" {
    inline for (all ++ unslotted) |T| try check3.interfaceCases(T, .gemmt);
}

test "interface: gemmt reports the first bad argument, and returns for empty problems" {
    const fa = testing.failing_allocator;
    const nan = std.math.nan(f64);
    var a: [16]f64 = @splat(nan);
    var c: [16]f64 = @splat(nan);
    try testing.expectError(error.InvalidUplo, gemmt(fa, check3.io, f64, .v, .o, .n, -1, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidTrans, gemmt(fa, check3.io, f64, .u, .o, .n, -1, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidTrans, gemmt(fa, check3.io, f64, .u, .n, .d, -1, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidM, gemmt(fa, check3.io, f64, .u, .n, .n, -1, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidK, gemmt(fa, check3.io, f64, .l, .n, .n, 0, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidLda, gemmt(fa, check3.io, f64, .l, .n, .n, 3, 2, 1, &a, 2, &a, 2, 0, &c, 3));
    try testing.expectError(error.InvalidLda, gemmt(fa, check3.io, f64, .l, .t, .n, 3, 4, 1, &a, 3, &a, 4, 0, &c, 3));
    try testing.expectError(error.InvalidLdb, gemmt(fa, check3.io, f64, .u, .n, .n, 3, 4, 1, &a, 3, &a, 3, 0, &c, 3));
    try testing.expectError(error.InvalidLdb, gemmt(fa, check3.io, f64, .u, .n, .c, 3, 4, 1, &a, 3, &a, 2, 0, &c, 3));
    try testing.expectError(error.InvalidLdc, gemmt(fa, check3.io, f64, .u, .n, .n, 3, 4, 1, &a, 3, &a, 4, 0, &c, 2));
    try testing.expectError(error.InvalidLdc, gemmt(fa, check3.io, f64, .u, .n, .n, 0, 4, 1, &a, 1, &a, 4, 0, &c, 0));
    try gemmt(fa, check3.io, f64, .u, .n, .n, 0, 3, nan, &a, 1, &a, 3, nan, &c, 1);
    try testing.expect(std.math.isNan(c[0]));
    // `k` 0 scales the triangle alone, and allocates nothing
    try gemmt(fa, check3.io, f64, .l, .n, .n, 2, 0, nan, &a, 2, &a, 1, 0, &c, 2);
    try testing.expectEqual(0, c[0]);
    try testing.expectEqual(0, c[1]);
    try testing.expect(std.math.isNan(c[2]));
    try testing.expectEqual(0, c[3]);
}

test "interface: symm and hemm match the wide reference" {
    inline for (all ++ unslotted) |T| {
        try check3.interfaceCases(T, .symm);
        if (comptime isComplex(T)) try check3.interfaceCases(T, .hemm);
    }
}

test "interface: symm and hemm report the first bad argument, and return for empty problems" {
    const fa = testing.failing_allocator;
    var a: [16]f64 = @splat(std.math.nan(f64));
    var c: [16]f64 = @splat(std.math.nan(f64));
    try testing.expectError(error.InvalidUplo, symm(fa, check3.io, f64, .l, .v, -1, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidM, symm(fa, check3.io, f64, .l, .u, -1, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidN, symm(fa, check3.io, f64, .l, .u, 0, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    // the order of `S` is `m` on the left, `n` on the right
    try testing.expectError(error.InvalidLda, symm(fa, check3.io, f64, .l, .u, 3, 0, 1, &a, 2, &a, 3, 0, &c, 3));
    try testing.expectError(error.InvalidLda, symm(fa, check3.io, f64, .r, .l, 1, 3, 1, &a, 2, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidLdb, symm(fa, check3.io, f64, .r, .u, 2, 1, 1, &a, 1, &a, 1, 0, &c, 2));
    try testing.expectError(error.InvalidLdc, symm(fa, check3.io, f64, .l, .u, 2, 1, 1, &a, 2, &a, 2, 0, &c, 1));
    try symm(fa, check3.io, f64, .l, .u, 0, 3, std.math.nan(f64), &a, 1, &a, 1, 0, &c, 1);
    try symm(fa, check3.io, f64, .r, .l, 3, 0, std.math.nan(f64), &a, 1, &a, 3, 0, &c, 3);
    try testing.expect(std.math.isNan(c[0]));
    var z: [4]types.Complex64 = @splat(.init(0, 0));
    try testing.expectError(error.InvalidLda, hemm(fa, check3.io, types.Complex64, .r, .u, 1, 2, .init(1, 0), &z, 1, &z, 1, .init(0, 0), &z, 1));
}

test "interface: gemm3m, symm3m and hemm3m match the wide reference" {
    inline for (all ++ unslotted) |T| {
        if (comptime isComplex(T)) {
            try check3.interfaceCases(T, .gemm3m);
            try check3.interfaceCases(T, .symm3m);
            try check3.interfaceCases(T, .hemm3m);
        }
    }
}

test "interface: the 3M routines check their arguments as gemm, symm and hemm" {
    const fa = testing.failing_allocator;
    const Z = types.Complex64;
    const one: Z = .init(1, 0);
    var z: [16]Z = @splat(.init(std.math.nan(f64), 0));
    try testing.expectError(error.InvalidTrans, gemm3m(fa, check3.io, Z, .n, .o, -1, -1, -1, one, &z, 0, &z, 0, one, &z, 0));
    try testing.expectError(error.InvalidK, gemm3m(fa, check3.io, Z, .c, .r, 0, 0, -1, one, &z, 0, &z, 0, one, &z, 0));
    try testing.expectError(error.InvalidLdb, gemm3m(fa, check3.io, Z, .n, .t, 3, 2, 4, one, &z, 3, &z, 1, one, &z, 3));
    try testing.expectError(error.InvalidLda, symm3m(fa, check3.io, Z, .r, .l, 1, 3, one, &z, 2, &z, 1, one, &z, 1));
    try testing.expectError(error.InvalidUplo, hemm3m(fa, check3.io, Z, .l, .v, 1, 1, one, &z, 1, &z, 1, one, &z, 1));
    // empty: nothing read, written or allocated
    try gemm3m(fa, check3.io, Z, .n, .n, 0, 3, 2, one, &z, 1, &z, 2, one, &z, 1);
    try hemm3m(fa, check3.io, Z, .r, .u, 3, 0, one, &z, 3, &z, 3, one, &z, 3);
    try testing.expect(std.math.isNan(z[0].re));
    // `k` 0 only scales `C`, and allocates nothing
    try gemm3m(fa, check3.io, Z, .n, .n, 2, 2, 0, one, &z, 2, &z, 1, .init(0, 0), &z, 2);
    try testing.expectEqual(Z.init(0, 0), z[0]);
}

test "interface: syrk, herk, syr2k and her2k match the wide reference" {
    inline for (all ++ unslotted) |T| {
        try check3.interfaceCases(T, .syrk);
        try check3.interfaceCases(T, .syr2k);
        if (comptime isComplex(T)) {
            try check3.interfaceCases(T, .herk);
            try check3.interfaceCases(T, .her2k);
        }
    }
}

test "interface: the rank-k updates report the first bad argument, and return for empty problems" {
    const fa = testing.failing_allocator;
    const Z = types.Complex64;
    var a: [16]f64 = @splat(std.math.nan(f64));
    var c: [16]f64 = @splat(std.math.nan(f64));
    var z: [4]Z = @splat(.init(0, 0));
    try testing.expectError(error.InvalidUplo, syrk(fa, check3.io, f64, .m, .r, -1, 0, 1, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidTrans, syrk(fa, check3.io, f64, .u, .r, -1, 0, 1, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidTrans, syr2k(fa, check3.io, f64, .u, .u, -1, 0, 1, &a, 1, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidTrans, syrk(fa, check3.io, Z, .u, .c, -1, 0, .init(1, 0), &z, 1, .init(0, 0), &z, 1));
    try testing.expectError(error.InvalidTrans, herk(fa, check3.io, Z, .l, .t, -1, 0, 1, &z, 1, 0, &z, 1));
    try testing.expectError(error.InvalidN, syrk(fa, check3.io, f64, .u, .c, -1, -1, 1, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidK, syrk(fa, check3.io, f64, .u, .n, 0, -1, 1, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidLda, syrk(fa, check3.io, f64, .u, .t, 1, 3, 1, &a, 2, 0, &c, 1));
    try testing.expectError(error.InvalidLdb, syr2k(fa, check3.io, f64, .l, .n, 3, 1, 1, &a, 3, &a, 2, 0, &c, 3));
    try testing.expectError(error.InvalidLdc, syr2k(fa, check3.io, f64, .l, .n, 3, 1, 1, &a, 3, &a, 3, 0, &c, 2));
    try testing.expectError(error.InvalidLdc, her2k(fa, check3.io, Z, .u, .c, 2, 0, .init(1, 0), &z, 1, &z, 1, 0, &z, 1));
    try syrk(fa, check3.io, f64, .u, .n, 0, 3, std.math.nan(f64), &a, 1, 0, &c, 1);
    try testing.expect(std.math.isNan(c[0]));
    // nothing to add and `beta` 1: a Hermitian diagonal keeps its imaginary part
    z[0] = .init(2, 5);
    try herk(fa, check3.io, Z, .u, .n, 1, 0, 1, &z, 1, 1, &z, 1);
    try testing.expectEqual(Z.init(2, 5), z[0]);
    try herk(fa, check3.io, Z, .u, .n, 1, 0, 1, &z, 1, 2, &z, 1);
    try testing.expectEqual(Z.init(4, 0), z[0]);
}

test "interface: trmm and trsm match the wide reference" {
    inline for (all) |T| {
        try check3.interfaceCases(T, .trmm);
        try check3.interfaceCases(T, .trsm);
    }
}

test "interface: unslotted trmm and trsm match the wide reference" {
    inline for (unslotted) |T| {
        try check3.interfaceCases(T, .trmm);
        try check3.interfaceCases(T, .trsm);
    }
}

test "interface: trmm and trsm report the first bad argument, and return for empty problems" {
    const fa = testing.failing_allocator;
    var a: [16]f64 = @splat(std.math.nan(f64));
    var b: [16]f64 = @splat(std.math.nan(f64));
    inline for (.{ trmm, trsm }) |f| {
        try testing.expectError(error.InvalidUplo, f(fa, check3.io, f64, .l, .v, .s, .n, -1, -1, 1, &a, 0, &b, 0));
        try testing.expectError(error.InvalidTrans, f(fa, check3.io, f64, .l, .u, .s, .n, -1, -1, 1, &a, 0, &b, 0));
        try testing.expectError(error.InvalidM, f(fa, check3.io, f64, .l, .u, .n, .n, -1, -1, 1, &a, 0, &b, 0));
        try testing.expectError(error.InvalidN, f(fa, check3.io, f64, .r, .l, .t, .u, 0, -1, 1, &a, 0, &b, 0));
        // the order of `A` is `m` on the left, `n` on the right
        try testing.expectError(error.InvalidLda, f(fa, check3.io, f64, .l, .u, .n, .n, 2, 0, 1, &a, 1, &b, 2));
        try testing.expectError(error.InvalidLda, f(fa, check3.io, f64, .r, .u, .c, .n, 0, 2, 1, &a, 1, &b, 1));
        try testing.expectError(error.InvalidLdb, f(fa, check3.io, f64, .r, .l, .n, .n, 2, 1, 1, &a, 1, &b, 1));
        try f(fa, check3.io, f64, .l, .u, .n, .n, 0, 3, std.math.nan(f64), &a, 1, &b, 1);
        try f(fa, check3.io, f64, .r, .u, .n, .n, 3, 0, std.math.nan(f64), &a, 1, &b, 3);
        try testing.expect(std.math.isNan(b[0]));
        // a zero `alpha` zeroes `B` without reading either
        try f(fa, check3.io, f64, .l, .l, .t, .u, 2, 2, 0, &a, 2, &b, 3);
        try testing.expectEqualSlices(f64, &.{ 0, 0 }, b[0..2]);
        try testing.expectEqualSlices(f64, &.{ 0, 0 }, b[3..5]);
        try testing.expect(std.math.isNan(b[2]));
        @memset(&b, std.math.nan(f64));
    }
}

test "the routines through Accelerate: special scalars, and special values in each operand" {
    if (comptime !accel.covers(f64)) return error.SkipZigTest;
    _ = blas.vtable();
    const was = accel.routing();
    defer _ = accel.route(was);
    if (!accel.route(true)) return error.SkipZigTest;
    var failed: ?anyerror = null;
    inline for (check.sd_types) |T| {
        inline for (.{ .gemm, .symm, .hemm, .syrk, .herk, .syr2k, .her2k, .trmm, .trsm }) |r| {
            if (comptime isComplex(T) or !(r == .hemm or r == .herk or r == .her2k)) {
                check3.accelerateCases(T, r) catch |e| {
                    failed = e;
                };
            }
        }
    }
    if (failed) |e| return e;
}
