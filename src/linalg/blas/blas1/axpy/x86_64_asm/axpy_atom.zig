//! Copyright (c) 2014, The OpenBLAS Project, under the BSD 3-clause license
//! reproduced in the OpenBLAS section of the repository `LICENSE` file.
//!
//! `daxpy` for the Atom core (`kernel/x86_64/axpy_atom.S`), as a naked function
//! whose body is that source with its macros as named fields, filled in by
//! `std.fmt.comptimePrint`.
//! `%` is doubled and local labels carry the kernel's name, for Zig's asm
//! syntax and so several kernels can share an object. `PREFETCHW` is not
//! defined for ATOM, so its blocks are absent.

const std = @import("std");
const sysinfo = @import("builtin");

const windows = sysinfo.os.tag == .windows or sysinfo.os.tag == .uefi;

/// `PROLOGUE`'s `.align`: 512 bytes for ELF, 16 for Win64 (common_x86_64.h).
const alignment = if (windows) 16 else 512;

/// The `#define`s of the source, `common.h` and `l1param.h` for ATOM in double
/// precision, filled into the `{[NAME]s}` fields below. `M`..`INCY` are
/// `ARG1`.. of the ABI (`#ifndef WINDOWS_ABI`).
const macros = .{
    .M = if (windows) "%%rcx" else "%%rdi",
    .X = if (windows) "%%rdx" else "%%rcx",
    .INCX = "%%r8",
    .Y = "%%r9",
    .INCY = if (windows) "%%r10" else "%%rsi",
    .YY = "%%r11",
    .ALPHA = "%%xmm15",
    .SIZE = "8",
    .PREFETCH = "prefetcht0",
    .PREFETCHSIZE = "( 64 * 12 + 32)",
    .PREOFFSET = "128",
    .ALIGN_3 = ".p2align 3",
    .SAVEREGISTERS = if (windows) save_registers else "",
    .RESTOREREGISTERS = if (windows) restore_registers else "",
};

/// Win64 makes `xmm6`-`xmm15` callee-saved (`SAVEREGISTERS` in common_x86_64.h).
const save_registers =
    "subq $256, %%rsp;" ++
    "movups %%xmm6, 0(%%rsp); movups %%xmm7, 16(%%rsp); movups %%xmm8, 32(%%rsp); movups %%xmm9, 48(%%rsp);" ++
    "movups %%xmm10, 64(%%rsp); movups %%xmm11, 80(%%rsp); movups %%xmm12, 96(%%rsp); movups %%xmm13, 112(%%rsp);" ++
    "movups %%xmm14, 128(%%rsp); movups %%xmm15, 144(%%rsp)";

const restore_registers =
    "movups 0(%%rsp), %%xmm6; movups 16(%%rsp), %%xmm7; movups 32(%%rsp), %%xmm8; movups 48(%%rsp), %%xmm9;" ++
    "movups 64(%%rsp), %%xmm10; movups 80(%%rsp), %%xmm11; movups 96(%%rsp), %%xmm12; movups 112(%%rsp), %%xmm13;" ++
    "movups 128(%%rsp), %%xmm14; movups 144(%%rsp), %%xmm15; addq $256, %%rsp";

/// The `#ifndef WINDOWS_ABI` split at the top of the source (`XDOUBLE` is not
/// defined).
const sysv_args =
    \\        movq     8(%%rsp), {[INCY]s}
    \\        movaps  %%xmm0,  {[ALPHA]s}
;

const win64_args =
    \\        movaps  %%xmm3,  {[ALPHA]s}
    \\
    \\        movq    40(%%rsp), {[X]s}
    \\        movq    48(%%rsp), {[INCX]s}
    \\        movq    56(%%rsp), {[Y]s}
    \\        movq    64(%%rsp), {[INCY]s}
;

const body =
    \\        {[SAVEREGISTERS]s}
    \\
    \\        leaq    (, {[INCX]s}, {[SIZE]s}), {[INCX]s}
    \\        leaq    (, {[INCY]s}, {[SIZE]s}), {[INCY]s}
    \\
    \\        testq   {[M]s}, {[M]s}
    \\        jle     .Ldaxpy_atom_29
    \\
    \\        cmpq    ${[SIZE]s}, {[INCX]s}
    \\        jne     .Ldaxpy_atom_20
    \\        cmpq    ${[SIZE]s}, {[INCY]s}
    \\        jne     .Ldaxpy_atom_20
    \\
    \\        movq    {[M]s},  %%rax
    \\        sarq    $3, %%rax
    \\        jle     .Ldaxpy_atom_13
    \\
    \\        movsd   0 * {[SIZE]s}({[X]s}), %%xmm0
    \\        movsd   1 * {[SIZE]s}({[X]s}), %%xmm1
    \\        movsd   2 * {[SIZE]s}({[X]s}), %%xmm2
    \\        movsd   3 * {[SIZE]s}({[X]s}), %%xmm3
    \\
    \\        movsd   0 * {[SIZE]s}({[Y]s}), %%xmm4
    \\        movsd   1 * {[SIZE]s}({[Y]s}), %%xmm5
    \\        movsd   2 * {[SIZE]s}({[Y]s}), %%xmm6
    \\        movsd   3 * {[SIZE]s}({[Y]s}), %%xmm7
    \\
    \\        movsd   4 * {[SIZE]s}({[X]s}), %%xmm8
    \\        mulsd   {[ALPHA]s}, %%xmm0
    \\        movsd   5 * {[SIZE]s}({[X]s}), %%xmm9
    \\        mulsd   {[ALPHA]s}, %%xmm1
    \\        movsd   6 * {[SIZE]s}({[X]s}), %%xmm10
    \\        mulsd   {[ALPHA]s}, %%xmm2
    \\        movsd   7 * {[SIZE]s}({[X]s}), %%xmm11
    \\        mulsd   {[ALPHA]s}, %%xmm3
    \\
    \\        decq    %%rax
    \\        jle .Ldaxpy_atom_12
    \\        {[ALIGN_3]s}
    \\
    \\.Ldaxpy_atom_11:
    \\        {[PREFETCH]s} ({[PREFETCHSIZE]s} +  0) - {[PREOFFSET]s}({[X]s})
    \\
    \\        addsd   %%xmm4, %%xmm0
    \\        movsd   4 * {[SIZE]s}({[Y]s}), %%xmm4
    \\        addsd   %%xmm5, %%xmm1
    \\        movsd   5 * {[SIZE]s}({[Y]s}), %%xmm5
    \\        addsd   %%xmm6, %%xmm2
    \\        movsd   6 * {[SIZE]s}({[Y]s}), %%xmm6
    \\        addsd   %%xmm7, %%xmm3
    \\        movsd   7 * {[SIZE]s}({[Y]s}), %%xmm7
    \\
    \\        movsd   %%xmm0, 0 * {[SIZE]s}({[Y]s})
    \\        mulsd   {[ALPHA]s}, %%xmm8
    \\        movsd    8 * {[SIZE]s}({[X]s}), %%xmm0
    \\
    \\        movsd   %%xmm1, 1 * {[SIZE]s}({[Y]s})
    \\        mulsd   {[ALPHA]s}, %%xmm9
    \\        movsd    9 * {[SIZE]s}({[X]s}), %%xmm1
    \\
    \\        movsd   %%xmm2, 2 * {[SIZE]s}({[Y]s})
    \\        mulsd   {[ALPHA]s}, %%xmm10
    \\        movsd   10 * {[SIZE]s}({[X]s}), %%xmm2
    \\
    \\        movsd   %%xmm3, 3 * {[SIZE]s}({[Y]s})
    \\        mulsd   {[ALPHA]s}, %%xmm11
    \\        movsd   11 * {[SIZE]s}({[X]s}), %%xmm3
    \\
    \\
    \\        addsd   %%xmm4, %%xmm8
    \\        movsd    8 * {[SIZE]s}({[Y]s}), %%xmm4
    \\        addsd   %%xmm5, %%xmm9
    \\        movsd    9 * {[SIZE]s}({[Y]s}), %%xmm5
    \\        addsd   %%xmm6, %%xmm10
    \\        movsd   10 * {[SIZE]s}({[Y]s}), %%xmm6
    \\        addsd   %%xmm7, %%xmm11
    \\        movsd   11 * {[SIZE]s}({[Y]s}), %%xmm7
    \\
    \\        movsd   %%xmm8,  4 * {[SIZE]s}({[Y]s})
    \\        mulsd   {[ALPHA]s}, %%xmm0
    \\        movsd   12 * {[SIZE]s}({[X]s}), %%xmm8
    \\
    \\        movsd   %%xmm9,  5 * {[SIZE]s}({[Y]s})
    \\        mulsd   {[ALPHA]s}, %%xmm1
    \\        movsd   13 * {[SIZE]s}({[X]s}), %%xmm9
    \\
    \\        movsd   %%xmm10, 6 * {[SIZE]s}({[Y]s})
    \\        mulsd   {[ALPHA]s}, %%xmm2
    \\        movsd   14 * {[SIZE]s}({[X]s}), %%xmm10
    \\
    \\        movsd   %%xmm11, 7 * {[SIZE]s}({[Y]s})
    \\        mulsd   {[ALPHA]s}, %%xmm3
    \\        movsd   15 * {[SIZE]s}({[X]s}), %%xmm11
    \\
    \\        addq    $8 * {[SIZE]s}, {[Y]s}
    \\        addq    $8 * {[SIZE]s}, {[X]s}
    \\        decq    %%rax
    \\        jg      .Ldaxpy_atom_11
    \\        {[ALIGN_3]s}
    \\
    \\.Ldaxpy_atom_12:
    \\        addsd   %%xmm4, %%xmm0
    \\        movsd   4 * {[SIZE]s}({[Y]s}), %%xmm4
    \\        addsd   %%xmm5, %%xmm1
    \\        movsd   5 * {[SIZE]s}({[Y]s}), %%xmm5
    \\        addsd   %%xmm6, %%xmm2
    \\        movsd   6 * {[SIZE]s}({[Y]s}), %%xmm6
    \\        addsd   %%xmm7, %%xmm3
    \\        movsd   7 * {[SIZE]s}({[Y]s}), %%xmm7
    \\
    \\        movsd   %%xmm0, 0 * {[SIZE]s}({[Y]s})
    \\        mulsd   {[ALPHA]s}, %%xmm8
    \\        movsd   %%xmm1, 1 * {[SIZE]s}({[Y]s})
    \\        mulsd   {[ALPHA]s}, %%xmm9
    \\        movsd   %%xmm2, 2 * {[SIZE]s}({[Y]s})
    \\        mulsd   {[ALPHA]s}, %%xmm10
    \\        movsd   %%xmm3, 3 * {[SIZE]s}({[Y]s})
    \\        mulsd   {[ALPHA]s}, %%xmm11
    \\
    \\        addsd   %%xmm4, %%xmm8
    \\        addsd   %%xmm5, %%xmm9
    \\        addsd   %%xmm6, %%xmm10
    \\        addsd   %%xmm7, %%xmm11
    \\
    \\        movsd   %%xmm8,  4 * {[SIZE]s}({[Y]s})
    \\        movsd   %%xmm9,  5 * {[SIZE]s}({[Y]s})
    \\        movsd   %%xmm10, 6 * {[SIZE]s}({[Y]s})
    \\        movsd   %%xmm11, 7 * {[SIZE]s}({[Y]s})
    \\
    \\        addq    $8 * {[SIZE]s}, {[X]s}
    \\        addq    $8 * {[SIZE]s}, {[Y]s}
    \\        {[ALIGN_3]s}
    \\
    \\.Ldaxpy_atom_13:
    \\        movq    {[M]s},  %%rax
    \\        andq    $4, %%rax
    \\        jle     .Ldaxpy_atom_15
    \\        {[ALIGN_3]s}
    \\
    \\        movsd   0 * {[SIZE]s}({[X]s}), %%xmm0
    \\        movsd   1 * {[SIZE]s}({[X]s}), %%xmm1
    \\        movsd   2 * {[SIZE]s}({[X]s}), %%xmm2
    \\        movsd   3 * {[SIZE]s}({[X]s}), %%xmm3
    \\
    \\        movsd   0 * {[SIZE]s}({[Y]s}), %%xmm4
    \\        mulsd   {[ALPHA]s}, %%xmm0
    \\        movsd   1 * {[SIZE]s}({[Y]s}), %%xmm5
    \\        mulsd   {[ALPHA]s}, %%xmm1
    \\        movsd   2 * {[SIZE]s}({[Y]s}), %%xmm6
    \\        mulsd   {[ALPHA]s}, %%xmm2
    \\        movsd   3 * {[SIZE]s}({[Y]s}), %%xmm7
    \\        mulsd   {[ALPHA]s}, %%xmm3
    \\
    \\        addsd   %%xmm4, %%xmm0
    \\        addsd   %%xmm5, %%xmm1
    \\        addsd   %%xmm6, %%xmm2
    \\        addsd   %%xmm7, %%xmm3
    \\
    \\        movsd   %%xmm0, 0 * {[SIZE]s}({[Y]s})
    \\        movsd   %%xmm1, 1 * {[SIZE]s}({[Y]s})
    \\        movsd   %%xmm2, 2 * {[SIZE]s}({[Y]s})
    \\        movsd   %%xmm3, 3 * {[SIZE]s}({[Y]s})
    \\
    \\        addq    $4 * {[SIZE]s}, {[X]s}
    \\        addq    $4 * {[SIZE]s}, {[Y]s}
    \\        {[ALIGN_3]s}
    \\
    \\.Ldaxpy_atom_15:
    \\        movq    {[M]s},  %%rax
    \\        andq    $2, %%rax
    \\        jle     .Ldaxpy_atom_16
    \\        {[ALIGN_3]s}
    \\
    \\        movsd   0 * {[SIZE]s}({[X]s}), %%xmm0
    \\        movsd   0 * {[SIZE]s}({[Y]s}), %%xmm4
    \\        movsd   1 * {[SIZE]s}({[X]s}), %%xmm1
    \\        movsd   1 * {[SIZE]s}({[Y]s}), %%xmm5
    \\
    \\        mulsd   {[ALPHA]s}, %%xmm0
    \\        mulsd   {[ALPHA]s}, %%xmm1
    \\        addsd   %%xmm4, %%xmm0
    \\        addsd   %%xmm5, %%xmm1
    \\
    \\        movsd   %%xmm0, 0 * {[SIZE]s}({[Y]s})
    \\        movsd   %%xmm1, 1 * {[SIZE]s}({[Y]s})
    \\
    \\        addq    $2 * {[SIZE]s}, {[X]s}
    \\        addq    $2 * {[SIZE]s}, {[Y]s}
    \\        {[ALIGN_3]s}
    \\
    \\.Ldaxpy_atom_16:
    \\        movq    {[M]s},  %%rax
    \\        andq    $1, %%rax
    \\        jle     .Ldaxpy_atom_19
    \\        {[ALIGN_3]s}
    \\
    \\        movsd   0 * {[SIZE]s}({[X]s}), %%xmm0
    \\        mulsd   {[ALPHA]s}, %%xmm0
    \\        addsd   0 * {[SIZE]s}({[Y]s}), %%xmm0
    \\
    \\        movsd   %%xmm0,  0 * {[SIZE]s}({[Y]s})
    \\        addq    ${[SIZE]s}, {[Y]s}
    \\        {[ALIGN_3]s}
    \\
    \\.Ldaxpy_atom_19:
    \\        xorq    %%rax,%%rax
    \\
    \\        {[RESTOREREGISTERS]s}
    \\
    \\        ret
    \\        {[ALIGN_3]s}
    \\
    \\.Ldaxpy_atom_20:
    \\        movq    {[Y]s}, {[YY]s}
    \\
    \\        movq    {[M]s},  %%rax
    \\        sarq    $3, %%rax
    \\        jle     .Ldaxpy_atom_23
    \\
    \\        movsd   ({[X]s}), %%xmm0
    \\        addq    {[INCX]s}, {[X]s}
    \\        movsd   ({[X]s}), %%xmm1
    \\        addq    {[INCX]s}, {[X]s}
    \\        movsd   ({[X]s}), %%xmm2
    \\        addq    {[INCX]s}, {[X]s}
    \\        movsd   ({[X]s}), %%xmm3
    \\        addq    {[INCX]s}, {[X]s}
    \\
    \\        movsd   ({[Y]s}), %%xmm4
    \\        addq    {[INCY]s}, {[Y]s}
    \\        movsd   ({[Y]s}), %%xmm5
    \\        addq    {[INCY]s}, {[Y]s}
    \\        movsd   ({[Y]s}), %%xmm6
    \\        addq    {[INCY]s}, {[Y]s}
    \\        movsd   ({[Y]s}), %%xmm7
    \\        addq    {[INCY]s}, {[Y]s}
    \\
    \\        movsd   ({[X]s}), %%xmm8
    \\        addq    {[INCX]s}, {[X]s}
    \\        mulsd   {[ALPHA]s}, %%xmm0
    \\        movsd   ({[X]s}), %%xmm9
    \\        addq    {[INCX]s}, {[X]s}
    \\        mulsd   {[ALPHA]s}, %%xmm1
    \\        movsd   ({[X]s}), %%xmm10
    \\        addq    {[INCX]s}, {[X]s}
    \\        mulsd   {[ALPHA]s}, %%xmm2
    \\        movsd   ({[X]s}), %%xmm11
    \\        addq    {[INCX]s}, {[X]s}
    \\        mulsd   {[ALPHA]s}, %%xmm3
    \\
    \\        decq    %%rax
    \\        jle .Ldaxpy_atom_22
    \\        {[ALIGN_3]s}
    \\
    \\.Ldaxpy_atom_21:
    \\        addsd   %%xmm4, %%xmm0
    \\        movsd   ({[Y]s}), %%xmm4
    \\        addq    {[INCY]s}, {[Y]s}
    \\        addsd   %%xmm5, %%xmm1
    \\        movsd   ({[Y]s}), %%xmm5
    \\        addq    {[INCY]s}, {[Y]s}
    \\
    \\        addsd   %%xmm6, %%xmm2
    \\        movsd   ({[Y]s}), %%xmm6
    \\        addq    {[INCY]s}, {[Y]s}
    \\        addsd   %%xmm7, %%xmm3
    \\        movsd   ({[Y]s}), %%xmm7
    \\        addq    {[INCY]s}, {[Y]s}
    \\
    \\        movsd   %%xmm0, ({[YY]s})
    \\        addq    {[INCY]s}, {[YY]s}
    \\        movsd   ({[X]s}), %%xmm0
    \\        addq    {[INCX]s}, {[X]s}
    \\        mulsd   {[ALPHA]s}, %%xmm8
    \\
    \\        movsd   %%xmm1, ({[YY]s})
    \\        addq    {[INCY]s}, {[YY]s}
    \\        movsd   ({[X]s}), %%xmm1
    \\        addq    {[INCX]s}, {[X]s}
    \\        mulsd   {[ALPHA]s}, %%xmm9
    \\
    \\        movsd   %%xmm2, ({[YY]s})
    \\        addq    {[INCY]s}, {[YY]s}
    \\        movsd   ({[X]s}), %%xmm2
    \\        addq    {[INCX]s}, {[X]s}
    \\        mulsd   {[ALPHA]s}, %%xmm10
    \\
    \\        movsd   %%xmm3, ({[YY]s})
    \\        addq    {[INCY]s}, {[YY]s}
    \\        movsd   ({[X]s}), %%xmm3
    \\        addq    {[INCX]s}, {[X]s}
    \\        mulsd   {[ALPHA]s}, %%xmm11
    \\
    \\        addsd   %%xmm4, %%xmm8
    \\        movsd   ({[Y]s}), %%xmm4
    \\        addq    {[INCY]s}, {[Y]s}
    \\        addsd   %%xmm5, %%xmm9
    \\        movsd   ({[Y]s}), %%xmm5
    \\        addq    {[INCY]s}, {[Y]s}
    \\
    \\        addsd   %%xmm6, %%xmm10
    \\        movsd   ({[Y]s}), %%xmm6
    \\        addq    {[INCY]s}, {[Y]s}
    \\        addsd   %%xmm7, %%xmm11
    \\        movsd   ({[Y]s}), %%xmm7
    \\        addq    {[INCY]s}, {[Y]s}
    \\
    \\        movsd   %%xmm8,  ({[YY]s})
    \\        addq    {[INCY]s}, {[YY]s}
    \\        movsd   ({[X]s}), %%xmm8
    \\        addq    {[INCX]s}, {[X]s}
    \\        mulsd   {[ALPHA]s}, %%xmm0
    \\
    \\        movsd   %%xmm9,  ({[YY]s})
    \\        addq    {[INCY]s}, {[YY]s}
    \\        movsd   ({[X]s}), %%xmm9
    \\        addq    {[INCX]s}, {[X]s}
    \\        mulsd   {[ALPHA]s}, %%xmm1
    \\
    \\        movsd   %%xmm10, ({[YY]s})
    \\        addq    {[INCY]s}, {[YY]s}
    \\        movsd   ({[X]s}), %%xmm10
    \\        addq    {[INCX]s}, {[X]s}
    \\        mulsd   {[ALPHA]s}, %%xmm2
    \\
    \\        movsd   %%xmm11, ({[YY]s})
    \\        addq    {[INCY]s}, {[YY]s}
    \\        movsd   ({[X]s}), %%xmm11
    \\        addq    {[INCX]s}, {[X]s}
    \\        mulsd   {[ALPHA]s}, %%xmm3
    \\
    \\        decq    %%rax
    \\        jg      .Ldaxpy_atom_21
    \\        {[ALIGN_3]s}
    \\
    \\.Ldaxpy_atom_22:
    \\        addsd   %%xmm4, %%xmm0
    \\        movsd   ({[Y]s}), %%xmm4
    \\        addq    {[INCY]s}, {[Y]s}
    \\        addsd   %%xmm5, %%xmm1
    \\        movsd   ({[Y]s}), %%xmm5
    \\        addq    {[INCY]s}, {[Y]s}
    \\        addsd   %%xmm6, %%xmm2
    \\        movsd   ({[Y]s}), %%xmm6
    \\        addq    {[INCY]s}, {[Y]s}
    \\        addsd   %%xmm7, %%xmm3
    \\        movsd   ({[Y]s}), %%xmm7
    \\        addq    {[INCY]s}, {[Y]s}
    \\
    \\        movsd   %%xmm0, ({[YY]s})
    \\        addq    {[INCY]s}, {[YY]s}
    \\        mulsd   {[ALPHA]s}, %%xmm8
    \\
    \\        movsd   %%xmm1, ({[YY]s})
    \\        addq    {[INCY]s}, {[YY]s}
    \\        mulsd   {[ALPHA]s}, %%xmm9
    \\
    \\        movsd   %%xmm2, ({[YY]s})
    \\        addq    {[INCY]s}, {[YY]s}
    \\        mulsd   {[ALPHA]s}, %%xmm10
    \\
    \\        movsd   %%xmm3, ({[YY]s})
    \\        addq    {[INCY]s}, {[YY]s}
    \\        mulsd   {[ALPHA]s}, %%xmm11
    \\
    \\        addsd   %%xmm4, %%xmm8
    \\        addsd   %%xmm5, %%xmm9
    \\        addsd   %%xmm6, %%xmm10
    \\        addsd   %%xmm7, %%xmm11
    \\
    \\        movsd   %%xmm8,  ({[YY]s})
    \\        addq    {[INCY]s}, {[YY]s}
    \\        movsd   %%xmm9,  ({[YY]s})
    \\        addq    {[INCY]s}, {[YY]s}
    \\        movsd   %%xmm10, ({[YY]s})
    \\        addq    {[INCY]s}, {[YY]s}
    \\        movsd   %%xmm11, ({[YY]s})
    \\        addq    {[INCY]s}, {[YY]s}
    \\        {[ALIGN_3]s}
    \\
    \\.Ldaxpy_atom_23:
    \\        movq    {[M]s},  %%rax
    \\        andq    $4, %%rax
    \\        jle     .Ldaxpy_atom_25
    \\        {[ALIGN_3]s}
    \\
    \\        movsd   ({[X]s}), %%xmm0
    \\        addq    {[INCX]s}, {[X]s}
    \\        movsd   ({[Y]s}), %%xmm4
    \\        addq    {[INCY]s}, {[Y]s}
    \\        movsd   ({[X]s}), %%xmm1
    \\        addq    {[INCX]s}, {[X]s}
    \\        movsd   ({[Y]s}), %%xmm5
    \\        addq    {[INCY]s}, {[Y]s}
    \\
    \\        movsd   ({[X]s}), %%xmm2
    \\        addq    {[INCX]s}, {[X]s}
    \\        mulsd   {[ALPHA]s}, %%xmm0
    \\        movsd   ({[Y]s}), %%xmm6
    \\        addq    {[INCY]s}, {[Y]s}
    \\        mulsd   {[ALPHA]s}, %%xmm1
    \\        movsd   ({[X]s}), %%xmm3
    \\        addq    {[INCX]s}, {[X]s}
    \\        mulsd   {[ALPHA]s}, %%xmm2
    \\        movsd   ({[Y]s}), %%xmm7
    \\        addq    {[INCY]s}, {[Y]s}
    \\        mulsd   {[ALPHA]s}, %%xmm3
    \\
    \\        addsd   %%xmm4, %%xmm0
    \\        addsd   %%xmm5, %%xmm1
    \\        addsd   %%xmm6, %%xmm2
    \\        addsd   %%xmm7, %%xmm3
    \\
    \\        movsd   %%xmm0, ({[YY]s})
    \\        addq    {[INCY]s}, {[YY]s}
    \\        movsd   %%xmm1, ({[YY]s})
    \\        addq    {[INCY]s}, {[YY]s}
    \\        movsd   %%xmm2, ({[YY]s})
    \\        addq    {[INCY]s}, {[YY]s}
    \\        movsd   %%xmm3, ({[YY]s})
    \\        addq    {[INCY]s}, {[YY]s}
    \\        {[ALIGN_3]s}
    \\
    \\.Ldaxpy_atom_25:
    \\        movq    {[M]s},  %%rax
    \\        andq    $2, %%rax
    \\        jle     .Ldaxpy_atom_26
    \\        {[ALIGN_3]s}
    \\
    \\        movsd   ({[X]s}), %%xmm0
    \\        addq    {[INCX]s}, {[X]s}
    \\        movsd   ({[Y]s}), %%xmm4
    \\        addq    {[INCY]s}, {[Y]s}
    \\        movsd   ({[X]s}), %%xmm1
    \\        addq    {[INCX]s}, {[X]s}
    \\        movsd   ({[Y]s}), %%xmm5
    \\        addq    {[INCY]s}, {[Y]s}
    \\
    \\        mulsd   {[ALPHA]s}, %%xmm0
    \\        mulsd   {[ALPHA]s}, %%xmm1
    \\        addsd   %%xmm4, %%xmm0
    \\        addsd   %%xmm5, %%xmm1
    \\
    \\        movsd   %%xmm0, ({[YY]s})
    \\        addq    {[INCY]s}, {[YY]s}
    \\        movsd   %%xmm1, ({[YY]s})
    \\        addq    {[INCY]s}, {[YY]s}
    \\        {[ALIGN_3]s}
    \\
    \\.Ldaxpy_atom_26:
    \\        movq    {[M]s},  %%rax
    \\        andq    $1, %%rax
    \\        jle     .Ldaxpy_atom_29
    \\        {[ALIGN_3]s}
    \\
    \\        movsd   ({[X]s}), %%xmm0
    \\        mulsd   {[ALPHA]s}, %%xmm0
    \\        addsd   ({[Y]s}), %%xmm0
    \\
    \\        movsd   %%xmm0,  ({[YY]s})
    \\        addq    ${[SIZE]s}, {[Y]s}
    \\        {[ALIGN_3]s}
    \\
    \\.Ldaxpy_atom_29:
    \\        xorq    %%rax, %%rax
    \\
    \\        {[RESTOREREGISTERS]s}
    \\
    \\        ret
;

const code = blk: {
    @setEvalBranchQuota(1 << 24);
    break :blk std.fmt.comptimePrint((if (windows) win64_args else sysv_args) ++ "\n\n" ++ body, macros);
};

/// The kernel; the export name is given where it is exported.
pub fn kernel() align(alignment) callconv(.naked) void {
    asm volatile (code);
}
