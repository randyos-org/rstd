//! The table-driven races of `lib2`, `lib3` and `lapack`.

const std = @import("std");
const rstd = @import("rstd");
const bb = @import("blas_bench.zig");
const kernels = @import("kernels.zig");
const libs = @import("libs.zig");
const rc = @import("ext/race.zig");
const e2 = @import("ext/level2.zig");
const e3 = @import("ext/level3.zig");
const el = @import("ext/lapack.zig");

const Io = std.Io;
const lapack = rstd.math.kernels.lapack;
const clock = bb.clock;
const l1_timing = bb.l1_timing;
const l3_timing = bb.l3_timing;
const precs = bb.precs;
const gemm_sizes = kernels.gemm_sizes;
const Col = libs.Col;
const lib_cols = libs.lib_cols;
const printHeader = libs.printHeader;
const sectionMark = libs.sectionMark;
const raceExt = rc.raceExt;
const Gbmv = e2.Gbmv;
const Ger = e2.Ger;
const Symv = e2.Symv;
const Trsv = e2.Trsv;
const Gemm3 = e3.Gemm3;
const Symm = e3.Symm;
const Syrk = e3.Syrk;
const Trmm = e3.Trmm;
const Trsm = e3.Trsm;
const Gbtrf = el.Gbtrf;
const Gehrd = el.Gehrd;
const Geqrf = el.Geqrf;
const Gesdd = el.Gesdd;
const Gesv = el.Gesv;
const Gesvd = el.Gesvd;
const Getrf = el.Getrf;
const Getrs = el.Getrs;
const Potrf = el.Potrf;
const Potri = el.Potri;
const Syevd = el.Syevd;
const SyevdN = el.SyevdN;
const Sytrf = el.Sytrf;
const Trtri = el.Trtri;

// ---- table-driven races ------------------------------------------------------
//
// The routines beyond `gemv` and `gemm` (`ext_l2`, `ext_l3`) and LAPACK
// (`lapack_ops`), one entry each, a struct of declarations:
//
// - `name`; `shape` (`Shape`); `fill` (`A`'s `Fill`, default `general`);
//   `prep` (`Prep`, default `none`); `restore` (what a call overwrites and
//   must find fresh again; default none, it accumulates); `out` (the `Buf`
//   checked, compared and hashed); `cols` (`out`'s columns where it is `b`
//   or `c` of an `mm` shape, default n); `level` (default its table's);
//   `band_rows` (a `band` shape's, default `2 band_k + 1`);
// - `flops(n)`, real; a complex call counts 4 times as many;
// - `sym(p)` (default `p.name ++ name`), `Sig(p)` and `lib(p, f, w)`: the
//   Fortran-ABI call (`lookupF`); a library without it shows `-`;
// - `rstd(p, w)`, absent until rstd has the routine; `has_rstd` gates one
//   written ahead of it;
// - `residual(p, w, i, j)`: `|want - got|` at output entry (i, j) and the
//   size of its terms, null off the stored part; without it only the
//   agreement with rstd is shown.
//
// Adding a routine is adding an entry to its table.

/// Where an entry's orders come from: `l2_bytes`, `gemm_sizes`,
/// `lapack_sizes`.
pub const Level = enum { l2, l3, lapack };

/// LAPACK's orders: `gemm_sizes` less the largest, where reference LAPACK
/// and `prepare` take many seconds a call.
pub const lapack_sizes = gemm_sizes[0 .. gemm_sizes.len - 1].*;

pub const Shape = enum {
    /// `A` n by n, `x` in `b` and `y` in `c`, n long; n from `l2_bytes`
    mv,
    /// `A` in band storage, `band_rows` by n; `x` and `y` as `mv`; n from
    /// `l2_bytes` over the rows
    band,
    /// `A`, `B` and `C` n by n
    mm,
};

/// `gbmv`'s and `gbtrf`'s `kl` and `ku`.
pub const band_k = 32;

pub const Fill = enum {
    /// every float in [-1, 1]
    general,
    /// `general` with 2n added to the diagonal's real parts: strictly
    /// diagonally dominant by rows and columns, so LU with partial pivoting
    /// exchanges no rows and every triangle is well conditioned
    dominant,
    /// `dominant` made Hermitian (symmetric if real), its diagonal real:
    /// positive definite
    spd,
};

pub const Prep = enum {
    none,
    /// `A` replaced by its Cholesky factor `U`, `A = Uᴴ U` (`prepare`)
    potrf,
    /// `A` replaced by `L` and `U` of `A = L U`, pivots 1 to n (`dominant`
    /// takes no row exchange; `prepare`)
    getrf,
};

/// `orig` is `A` before `prep`; `a0`, `b0` and `c0` what each call starts
/// from; `s` singular values, eigenvalues or `tau`.
pub const Buf = enum { a, a0, orig, b, b0, c, c0, s, u, vt };

pub const Restore = struct { a: bool = false, b: bool = false, c: bool = false };

pub const Cols = enum { one, n };

/// The table-driven races' contestants: the libraries' Fortran-ABI routines
/// and rstd's public ones. `max_n` holds in level 3 and LAPACK only.
pub const ext_cols = lib_cols[0..4].* ++ [_]Col{
    .{ .label = "reference", .short = "ref", .lib = .reference, .max_n = 1024 },
    .{ .label = "eigen", .short = "eigen", .lib = .eigen },
    .{ .label = "aocl-libflame", .short = "aocl", .lib = .aoclflame },
};

pub const ext_l2 = .{ Ger, Symv, Trsv, Gbmv };
pub const ext_l3 = .{ Gemm3, Symm, Syrk, Trmm, Trsm };
/// rstd's LAPACK. An entry gates rstd's call on `@hasDecl`, so it races
/// once the routine lands; `Geqrf` and `Syevd` race the libraries alone until
/// their rstd call is added (issue #120).
pub const lapack_ops = .{ Potrf, Getrf, Getrs, Gesv, Trtri, Potri, Gesvd, Gesdd, Geqrf, Gehrd, Syevd, SyevdN, Sytrf, Gbtrf };

pub fn libLapack(io: Io, gpa: std.mem.Allocator) !void {
    sectionMark("lapack");
    try clock(io, "lapack start");
    try extTable(io, gpa, .lapack, lapack_ops);
}

pub fn extTable(io: Io, gpa: std.mem.Allocator, comptime lv: Level, comptime ops: anytype) !void {
    const t = if (lv == .l2) l1_timing else l3_timing;
    std.debug.print("\n{s} GFLOP/s, best of {d} batches of at least {d} ms; ratios are over OpenBLAS. {s} {s}\n\n", .{ ext_leads[@backingInt(lv)], t.batches, @divTrunc(t.min_ns, std.time.ns_per_ms), libs.io_note, ext_note });
    printHeader("routine | p | n", "------- | - | -:", &ext_cols);
    inline for (ops) |E| {
        if (libs.wanted(E.name)) inline for (precs) |p| try raceExt(io, gpa, E, lv, p);
    }
}

pub const ext_leads = [_][]const u8{
    "Level 2 beyond `gemv`, each library's Fortran-ABI routine, unit strides, n as `gemv`'s (`A` of 32 KiB, 4 MiB, 64 MiB): `ger` `A += 1e-6 x yᵀ` (`geru` complex), `symv` `y += 1e-3 S x` (upper), `trsv` upper, no transpose, non-unit on a diagonally dominant `A`, `x` restored before each call (counted), `gbmv` `y += 1e-3 A x` with kl = ku = 32 (`A`'s band of the same bytes); rates by 2n², 2n², n², 2n (kl + ku + 1), four times as many complex.",
    "Level 3 beyond `gemm`, each library's Fortran-ABI routine, square n, `alpha` and `beta` 1: `gemm` NN, `symm` left upper, `syrk` upper NN, `trmm` and `trsm` left upper NN non-unit (`trsm` on a diagonally dominant `A`), `B` restored before each `trmm` and `trsm` call (an n² copy, counted); rates by 2n³, 2n³, n²(n + 1), n³, n³, four times as many complex; `reference` to n = 1024.",
    "LAPACK, square n, upper where it matters: rates by LAPACK's counts, `potrf` n³/3, `getrf` 2n³/3, `getrs` 2n³ (n right-hand sides), `gesv` 2n³/3 + 2n² (one), `trtri` n³/3, `potri` 2n³/3, and nominal ones: `gesvd` and `gesdd` 22n³ (thin `U` and `Vᵀ`, Golub and Van Loan's Golub-Reinsch count), `geqrf` 4n³/3, `gehrd` 10n³/3 (the whole of `A`), `syevd` (`heevd` complex) 9n³ with vectors and `syevdN` 4n³/3 without (the reduction's), `sytrf` n³/3, `gbtrf` 4nk² (kl = ku = k = 32, at `gemv`'s band sizes); four times as many complex. `A` is diagonally dominant for `getrf`, `getrs`, `gesv` and `trtri` (no row exchanges, so pivots agree), Hermitian positive definite for `potrf`, `potri`, `syevd` and `sytrf`; `getrs` and `potri` start from a factorization made here in f64. Each call restores what it overwrites first (an n² copy, counted). `gesvd`, `gesdd`, `geqrf`, `syevd`, `sytrf` and `gbtrf` race the libraries until rstd has them.",
};

pub const ext_note = "`rstd err` is rstd's largest residual at the checked entries (all, or 256 sampled), in eps of the type relative to the size of its terms; `rstd-<lib>` is rstd's largest difference from that library's output in eps relative to the output's largest entry, `!` past 2 (terms + 8). `rstd bits` hashes rstd's output. Eigen runs one thread; `-` where a library lacks the routine";
