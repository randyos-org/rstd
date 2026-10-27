//! Element and index types.
//!
//! Sizes, strides and leading dimensions are `isize` at the interface, so a
//! call can address any array that fits in memory; this is not the reference's
//! Fortran `INTEGER` (`blasint`), since there is no Fortran or C ABI here.
//! `BlasLong` (`BLASLONG`) is the kernel index type, pointer-sized as `long`
//! is on LP64 and as the reference's `long long` is on 64-bit Windows.

const std = @import("std");
const sysinfo = @import("builtin");

const SrcInfo = std.lang.SourceLocation;

/// Kernel-level index/stride type: `long`, pointer-sized (LP64, and Win64's
/// `long long`), the same as the interface's `isize`.
pub const BlasLong = isize;

/// Whether extended-precision routines use `f128`. `QUAD_PRECISION` in the
/// reference, which builds it as a software type and rules out
/// `EXPRECISION`; Zig's `f128` (compiler-rt arithmetic) plays that role.
/// Opt-in: off unless set here.
pub const quad_precision = false;

/// `EXPRECISION`: the `q`/`x` families compute in x87 `long double`. The
/// reference sets it for x86 and x86_64 built with gcc or clang, never for a
/// Microsoft-ABI target (`Makefile.system`); here likewise (`exprecisionOn`).
/// Elsewhere `q`/`x` alias `d`/`z`.
pub const exprecision = !quad_precision and exprecisionOn(
    sysinfo.target.cpu.arch,
    sysinfo.target.os.tag,
    @typeInfo(c_longdouble).float.bits,
    sysinfo.target.cpu.arch.isX86() and std.Target.x86.featureSetHas(sysinfo.target.cpu.features, .x87),
);

/// Whether a target gets `EXPRECISION`: x86 or x86_64 with an x87 unit, whose
/// `long double` is the 80-bit x87 type, and not a Microsoft-ABI OS (Windows,
/// UEFI), which the reference never builds it for. Assumes the x87 control
/// word's extended precision (0x037F), as Linux and x64 UEFI's handoff leave
/// it; Windows sets double (0x027F). An OS without either must set it.
pub fn exprecisionOn(arch: std.Target.Cpu.Arch, os: std.Target.Os.Tag, longdouble_bits: u16, x87: bool) bool {
    if (!arch.isX86() or !x87 or longdouble_bits != 80) return false;
    return switch (os) {
        .windows, .uefi => false,
        else => true,
    };
}

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

/// `XDouble` (`xdouble`): `c_longdouble` (x87 `f80`) with `exprecision`,
/// `f128` with `quad_precision`, else `f64`. `f80` keeps Zig's native size: 16
/// bytes on x86_64, 12 on x86, where the reference builds with
/// `-m128bit-long-double` (`SIZE` 16). Not an ABI drop-in, so nothing depends
/// on 16; buffer formulas count bytes and round up to whole elements. The gemm
/// blocking (`blas3/gemm/params.zig`) keeps the reference's 16 and 32, which
/// only make its blocks a little smaller than they could be at 12 bytes.
pub const XDouble = if (quad_precision) f128 else if (exprecision) c_longdouble else f64;

/// Whether float `F` is wider than `f64` (`f80`, `f128`): such types never
/// fuse (`fma.zig`), and the portable `q`/`x` reductions and gemv take the
/// x87 kernels' accumulator order.
pub fn isExtended(comptime F: type) bool {
    return @typeInfo(F).float.bits > 64;
}

/// Complex counterparts, matching the `c`/`z`/`x` families.
pub const Complex32 = std.math.Complex(f32);
pub const Complex64 = std.math.Complex(f64);
pub const ComplexXD = std.math.Complex(XDouble);

/// True for `std.math.Complex` of any float type.
pub fn isComplex(comptime T: type) bool {
    return @typeInfo(T) == .@"struct" and @hasField(T, "re") and T == std.math.Complex(@FieldType(T, "re"));
}

/// The float type of `T`'s parts: `T` itself for a real type.
pub fn Float(comptime T: type) type {
    return if (isComplex(T)) @FieldType(T, "re") else T;
}

/// Whether the vtable has slots over `F` (`f32`, `f64`, `XDouble`).
pub fn hasSlots(comptime F: type) bool {
    return F == f32 or F == f64 or F == XDouble;
}

/// Calling convention of the portable kernels over `F`: C for the vtable's
/// types, Zig's own for the rest, which only Zig calls.
pub fn kernelCc(comptime F: type) std.lang.CallingConvention {
    return if (hasSlots(F)) .c else .auto;
}

/// The type the portable kernels accumulate `F` in: `f32` below it, so a
/// long sum of `f16` neither stalls nor overflows; `F` itself otherwise.
pub fn Acc(comptime F: type) type {
    return if (@bitSizeOf(F) < 32) f32 else F;
}

/// C-ABI complex values, for kernel signatures that return one (the `dotu`/
/// `dotc` kernels): `std.math.Complex` has automatic layout, which a
/// `callconv(.c)` return type cannot have. An extern struct, not C's
/// `_Complex`: `CComplexXD` at `f80` returns through a hidden pointer on
/// SysV, not in `st(0)`/`st(1)` as the x87 `xdotu`/`xdotc` do, so only Zig
/// kernels fill those slots.
pub fn CComplex(comptime F: type) type {
    return extern struct {
        re: F,
        im: F,

        pub fn toComplex(c: @This()) std.math.Complex(F) {
            return .{ .re = c.re, .im = c.im };
        }
    };
}

/// A `CComplex` or `std.math.Complex` value as the latter.
pub fn complexOf(v: anytype) std.math.Complex(@TypeOf(v.re)) {
    return .init(v.re, v.im);
}

pub const CComplex32 = CComplex(f32);
pub const CComplex64 = CComplex(f64);
pub const CComplexXD = CComplex(XDouble);

/// Every element type a BLAS routine can be defined over.
///
/// Note these are not always distinct: without `exprecision` or
/// `quad_precision`, `XDouble` *is* `f64` and `ComplexXD` *is* `Complex64`,
/// so the `q`/`x` families alias `d`/`z`. That is the reference's behavior too, and it is why this is a
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

// ---- argument enums ----------------------------------------------------------

/// `TRANS`: `n` is `A`, `t` `A^T`, `r` `conj(A)`, `c` `A^H`; real `r` and
/// `c` are `n` and `t`. `o`, `u`, `s`, `d` are those four with `x`
/// conjugated too (complex `gemv` and `gbmv` only). A routine given a letter
/// it does not take returns `error.InvalidTrans`.
pub const Trans = enum {
    n,
    t,
    r,
    c,
    o,
    u,
    s,
    d,

    pub fn transposed(t: Trans) bool {
        return switch (t) {
            .t, .c, .u, .d => true,
            .n, .r, .o, .s => false,
        };
    }

    /// Whether `A` is conjugated; only complex types have conjugates.
    pub fn conjA(t: Trans) bool {
        return switch (t) {
            .r, .c, .s, .d => true,
            .n, .t, .o, .u => false,
        };
    }

    /// Whether `x` is conjugated (`o`, `u`, `s`, `d`).
    pub fn conjX(t: Trans) bool {
        return switch (t) {
            .o, .u, .s, .d => true,
            .n, .t, .r, .c => false,
        };
    }
};

/// `UPLO`: the stored triangle, upper or lower. `v` and `m` are the upper and
/// lower triangles of the transposed Hermitian matrix (OpenBLAS's
/// `HEMVREV`), taken by the Hermitian level-2 routines only; the others
/// return `error.InvalidUplo` for them.
pub const Uplo = enum {
    u,
    l,
    v,
    m,

    pub fn isLower(uplo: Uplo) bool {
        return uplo == .l or uplo == .m;
    }

    pub fn isRev(uplo: Uplo) bool {
        return uplo == .v or uplo == .m;
    }
};

/// `DIAG`: a unit (`u`) diagonal is taken as ones and never read.
pub const Diag = enum { u, n };

/// `SIDE`: whether the symmetric or triangular matrix multiplies from the
/// left (`l`) or the right (`r`).
pub const Side = enum { l, r };
