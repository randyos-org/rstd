//! Where rstd's kernels run instead of Accelerate's for speed alone, by
//! architecture (design.md, Apple Accelerate, "Kept for speed"): a slot is
//! kept where routing ran under 0.97 of rstd's own past noise. arm64: the
//! flyoff of 2026-10-29 for level 1, the M6 sweep of 2026-11-25 (issue #726)
//! for levels 2 and 3, by x86_64's rule. x86_64: the Intel Mac's sweep of
//! 2026-11-26 (issue #772), each slot's median of `rstd/ob` routed
//! (`--threads=6`) over `--rstd-accelerate=0`, the one-thread runs deciding
//! the slots within noise; the slots it left, the `accroute` bench section's
//! runs of 2026-11-27 (issue #797), routed against rstd's own in one process.
//! Rerun timings revise the tables, not the code.

const std = @import("std");
const portable = @import("../portable.zig");
const binding = @import("binding.zig");

const x86_64 = binding.x86_64;

/// The level-1 routines `sends` covers; `dot` over `f32` is `sdot`.
pub const L1 = enum { copy, swap, scal, dot, dotc, dsdot, nrm2, asum, rot, rotm, axpy, axpby };

/// Which calls go to Accelerate: all, those whose every stride is 1, those
/// with another stride, none.
const Strides = enum { all, unit, strided, none };

/// `Strides` by letter.
const Letters = struct { s: Strides = .all, d: Strides = .all, c: Strides = .all, z: Strides = .all };

/// The level-1 table; on arm64 the slots the flyoff did not time (`rotm`,
/// `dsdot`, `axpby`) stay routed.
const sends: std.enums.EnumFieldStruct(L1, Letters, .{}) = if (x86_64) .{
    .asum = .{ .c = .strided, .z = .unit },
    .axpy = .{ .s = .unit, .c = .none, .z = .strided },
    .axpby = .{ .s = .unit },
    .copy = .{ .s = .unit, .d = .unit, .c = .unit },
    .dot = .{ .z = .none },
    .dotc = .{ .z = .none },
    .nrm2 = .{ .s = .none, .d = .none, .c = .none, .z = .none },
    .rot = .{ .s = .strided, .d = .unit },
    .scal = .{ .s = .none, .d = .unit, .c = .strided, .z = .strided },
    .swap = .{ .d = .unit, .c = .unit },
} else .{
    .nrm2 = .{ .s = .none, .d = .none, .c = .none, .z = .none },
    .dot = .{ .s = .none, .c = .unit, .z = .unit },
    .dotc = .{ .c = .unit, .z = .unit },
    .copy = .{ .s = .unit, .d = .unit },
    .swap = .{ .s = .unit, .d = .unit, .c = .unit, .z = .unit },
    .axpy = .{ .c = .unit, .z = .unit },
    .rot = .{ .c = .unit },
};

/// Whether `sends` lets `r` over `T` go to Accelerate; `unit` if every
/// stride is 1.
pub inline fn sent(comptime r: L1, comptime T: type, unit: bool) bool {
    return switch (comptime @field(@field(sends, @tagName(r)), portable.letter(T).?)) {
        .all => true,
        .unit => unit,
        .strided => !unit,
        .none => false,
    };
}

/// The level-2 and level-3 routines by their reference names, the Hermitian
/// ones apart: `gemv` both ways, `ger` (`gerc` follows it), `syrk` and `herk`
/// without `B`, `syr2k` and `her2k` with it.
pub const L23 = enum { gemv, gbmv, ger, symv, hemv, sbmv, hbmv, spmv, hpmv, syr, her, spr, hpr, syr2, her2, spr2, hpr2, trmv, trsv, tbmv, tbsv, tpmv, tpsv, gemm, symm, hemm, syrk, herk, syr2k, her2k, trmm, trsm };

/// Whether rstd's kernels run, by letter.
const Kept = struct { s: bool = false, d: bool = false, c: bool = false, z: bool = false };

const all: Kept = .{ .s = true, .d = true, .c = true, .z = true };
const real: Kept = .{ .s = true, .d = true };
const complex: Kept = .{ .c = true, .z = true };

/// The level-2 and level-3 table.
const kept: std.enums.EnumFieldStruct(L23, Kept, .{}) = if (x86_64) .{
    .gemv = all,
    .gbmv = all,
    .ger = all,
    .symv = all,
    .hemv = complex,
    .sbmv = real,
    .hbmv = complex,
    .spmv = real,
    .hpmv = complex,
    .syr = real,
    .her = complex,
    .spr = .{ .s = true },
    .hpr = complex,
    .syr2 = real,
    .her2 = complex,
    .spr2 = .{ .s = true },
    .hpr2 = complex,
    .trmv = .{ .z = true },
    .trsv = complex,
    .tbmv = all,
    .tbsv = complex,
    .tpmv = all,
    .tpsv = complex,
    .gemm = all,
    .symm = .{ .s = true, .c = true, .z = true },
    .hemm = complex,
    .syrk = all,
    .herk = complex,
    .syr2k = all,
    .her2k = complex,
    .trmm = complex,
    .trsm = complex,
} else .{
    .gbmv = all,
    .ger = complex,
    .symv = all,
};

/// Whether `kept` lets `r` over `T` go to Accelerate.
pub inline fn routes(comptime r: L23, comptime T: type) bool {
    return comptime !@field(@field(kept, @tagName(r)), portable.letter(T).?);
}
