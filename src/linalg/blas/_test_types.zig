//! Tests for `types.zig`. They live here because `types.zig` is part of the
//! per-core kernel objects, which are compiled with only `linalg/blas/`
//! importable, and `zigstd` is outside that.

const std = @import("std");
const zigstd = @import("../../zigstd.zig");
const types = @import("types.zig");

const BF16 = types.BF16;
const isBlasType = types.isBlasType;
const HF16 = types.HF16;
const XDouble = types.XDouble;
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;
const ComplexXD = types.ComplexXD;

test "BF16 is a distinct 16-bit type, not an arithmetic one" {
    try zigstd.testExpectEqual(@as(usize, 2), @sizeOf(BF16));
    try zigstd.testExpectTrue(BF16 != u16);
    try zigstd.testExpectTrue(isBlasType(BF16));
}

test "isBlasType covers the full family and rejects everything else" {
    inline for (.{ BF16, HF16, f32, f64, XDouble, Complex32, Complex64, ComplexXD }) |T| {
        try zigstd.testExpectTrue(isBlasType(T));
    }
    inline for (.{ u8, i32, bool, usize, []const u8 }) |T| {
        try zigstd.testExpectTrue(!isBlasType(T));
    }
}
