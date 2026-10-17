//! Element and index types.
//!
//! Sizes, strides and leading dimensions are `isize` at the interface, so a
//! call can address any array that fits in memory; this is not the reference's
//! Fortran `INTEGER` (`blasint`), since there is no Fortran or C ABI here.
//! `BlasLong` (`BLASLONG`) is the kernel index type, pointer-sized as `long`
//! is on LP64 and as the reference's `long long` is on 64-bit Windows.

const std = @import("std");
const sysinfo = @import("builtin");

const SrcInfo = std.builtin.SourceLocation;

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
    sysinfo.cpu.arch,
    sysinfo.os.tag,
    @typeInfo(c_longdouble).float.bits,
    sysinfo.cpu.arch.isX86() and std.Target.x86.featureSetHas(sysinfo.cpu.features, .x87),
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
pub fn kernelCc(comptime F: type) std.builtin.CallingConvention {
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
