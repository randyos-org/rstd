//! Element and index types.
//!
//! Sizes, strides and leading dimensions are `isize` at the interface, so a
//! call can address any array that fits in memory; this is not the reference's
//! Fortran `INTEGER` (`blasint`), since there is no Fortran or C ABI here.
//! `BlasLong` (`BLASLONG`) is the kernel index type, pointer-sized as `long`
//! is on LP64 and as the reference's `long long` is on 64-bit Windows.

const std = @import("std");

const SrcInfo = std.builtin.SourceLocation;

/// Kernel-level index/stride type: `long`, so pointer-sized where the pointer
/// is 64 bits (LP64 and Win64's `long long`) and 32 bits on 32-bit targets.
pub const BlasLong = if (@bitSizeOf(usize) >= 64) i64 else i32;

/// Whether extended-precision routines use `f80` rather than `f64`.
/// `EXPRECISION` in the reference; off by default, where `XDouble` is just
/// `f64` and the `q`/`x` families are aliases of `d`/`z`.
pub const exprecision = false;

/// Whether extended-precision routines use `f128`. `QUAD_PRECISION` in the
/// reference, which builds it as a software type and rules out
/// `EXPRECISION`; Zig's `f128` (compiler-rt arithmetic) plays that role.
pub const quad_precision = false;

comptime {
    if (quad_precision and exprecision) @compileError("quad_precision and exprecision are mutually exclusive");
}

// ---- element types ---------------------------------------------------------
// The full set the reference supports, not just the classic four. Which of
// these are actually built is a per-build choice there (`BUILD_SINGLE`,
// `BUILD_BFLOAT16`, `BUILD_HFLOAT16`, `QUAD_PRECISION`, ...).

/// `BF16`. A distinct type over a 16-bit payload rather than an
/// arithmetic type, because that is what it is in the reference too: a
/// `uint16_t` bit pattern with explicit conversion kernels
/// (`sbstobf16_k`, `sbf16tos_k`, ...). Nothing computes in bf16 directly --
/// it is a storage format that gets widened to `f32` to do arithmetic, so
/// making it non-arithmetic here keeps that honest.
pub const BF16 = enum(u16) { _ };

/// IEEE half. `_Float16` in the reference, under `BUILD_HFLOAT16`.
pub const HF16 = f16;

/// `XDouble`: `f64` in a default build, `f80` with `exprecision`, `f128`
/// with `quad_precision`.
pub const XDouble = if (quad_precision) f128 else if (exprecision) f80 else f64;

/// Complex counterparts, matching the `c`/`z`/`x` families.
pub const Complex32 = std.math.Complex(f32);
pub const Complex64 = std.math.Complex(f64);
pub const ComplexXD = std.math.Complex(XDouble);

/// C-ABI complex values, for kernel signatures that return one (the `dotu`/
/// `dotc` kernels): `std.math.Complex` has automatic layout, which a
/// `callconv(.c)` return type cannot have.
pub fn CComplex(comptime F: type) type {
    return extern struct {
        re: F,
        im: F,

        pub fn toComplex(c: @This()) std.math.Complex(F) {
            return .{ .re = c.re, .im = c.im };
        }
    };
}
pub const CComplex32 = CComplex(f32);
pub const CComplex64 = CComplex(f64);
pub const CComplexXD = CComplex(XDouble);

/// Every element type a BLAS routine can be defined over.
///
/// Note these are not all distinct: with both extended modes off, `XDouble` *is*
/// `f64` and `ComplexXD` *is* `Complex64`, so the `q`/`x` families alias
/// `d`/`z`. That is the reference's behavior too, and it is why this is a
/// membership test rather than a `switch` -- a switch cannot have the same
/// type in two prongs.
pub const blas_types = .{
    BF16,    HF16,      f32,       f64,
    XDouble, Complex32, Complex64, ComplexXD,
};

/// True for any type a BLAS routine can have as its element type.
pub fn isBlasType(comptime T: type) bool {
    inline for (blas_types) |Candidate| {
        if (T == Candidate) return true;
    }
    return false;
}

/// Asserts `T` is an element type BLAS routines are defined over -- the
/// full set the reference supports, not just the classic four: `BF16` and
/// `HF16` alongside `f32`/`f64`/`XDouble` and their complex counterparts.
/// See `runtime/types.zig`.
///
/// This says nothing about whether a given routine implements `T`; that is
/// each routine's own business, and should be its own compile error naming
/// what is missing. Call with `comptime assertBlasType(T, @src())` so the
/// message names the function that got the bad type.
pub fn assertBlasType(comptime T: type, comptime src: SrcInfo) void {
    if (!isBlasType(T)) {
        @compileError(src.fn_name ++ ": " ++ @typeName(T) ++ " is not a BLAS element type" ++
            " (expected BF16, HF16, f32, f64, XDouble, or a Complex of one)");
    }
}
