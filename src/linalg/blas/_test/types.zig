//! Tests for `types.zig`. They live here because `types.zig` is part of the
//! per-core kernel objects, which are compiled with only `linalg/blas/`
//! importable, and `zigstd` is outside that.

const std = @import("std");
const zigstd = @import("../../../zigstd.zig");
const types = @import("../types.zig");

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

test "XDouble is the x87 long double exactly where the reference enables EXPRECISION" {
    const sysinfo = @import("builtin");
    const x86 = sysinfo.cpu.arch == .x86 or sysinfo.cpu.arch == .x86_64;
    if (types.exprecision) {
        try zigstd.testExpectTrue(x86 and sysinfo.os.tag != .windows);
        try zigstd.testExpectEqual(@as(u16, 80), @typeInfo(XDouble).float.bits);
        try zigstd.testExpectEqual(@sizeOf(c_longdouble), @sizeOf(XDouble));
        try zigstd.testExpectTrue(XDouble != f64 and ComplexXD != Complex64);
    } else if (!types.quad_precision) {
        try zigstd.testExpectTrue(XDouble == f64 and ComplexXD == Complex64);
    }
    if (sysinfo.os.tag == .windows or sysinfo.os.tag == .uefi or !x86) try zigstd.testExpectTrue(!types.exprecision);
    try zigstd.testExpectTrue(types.isExtended(f80) and types.isExtended(f128) and !types.isExtended(f64));
}

test "EXPRECISION by target: x87 long double, never a Microsoft ABI" {
    const on = types.exprecisionOn;
    try zigstd.testExpectTrue(on(.x86_64, .linux, 80, true));
    try zigstd.testExpectTrue(on(.x86, .linux, 80, true));
    try zigstd.testExpectTrue(on(.x86_64, .macos, 80, true));
    try zigstd.testExpectTrue(on(.x86_64, .freestanding, 80, true));
    // Microsoft ABIs, even with an 80-bit long double (windows-gnu, uefi-gnu)
    inline for (.{ std.Target.Cpu.Arch.x86, .x86_64 }) |arch| {
        try zigstd.testExpectTrue(!on(arch, .windows, 80, true));
        try zigstd.testExpectTrue(!on(arch, .uefi, 80, true));
    }
    // no x87 unit, or a long double that is not x87's
    try zigstd.testExpectTrue(!on(.x86_64, .linux, 80, false));
    try zigstd.testExpectTrue(!on(.x86_64, .linux, 64, true));
    try zigstd.testExpectTrue(!on(.aarch64, .linux, 128, false));
    try zigstd.testExpectTrue(!on(.aarch64, .macos, 64, false));
}
