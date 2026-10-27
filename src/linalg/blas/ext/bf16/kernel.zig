//! The bfloat16 and half-precision kernels: the conversions, `bscal`, `sbdot`
//! and the `bgemv`/`sbgemv`/`shgemv` family. Nothing here depends on the BLAS
//! root, so the per-core kernel objects (`_kernel.zig`) can use it too.
//!
//! Every bfloat16 conversion, on every core, is `kernel/x86_64/{tobf16,bf16to}.c`
//! (`Makefile.L1` selects them for all architectures; only the AVX-512 BF16
//! microkernels of COOPERLAKE and SAPPHIRERAPIDS differ, and those give the
//! same bits, so there is one portable kernel). Denormals flush to a signed
//! zero, NaNs come out quiet, the rest round to nearest even.
//!
//! These are portable kernels, the same on every core (see the readme,
//! "bfloat16 and half precision"). The loops are plain, for the compiler to
//! vectorize.

const std = @import("std");
const types = @import("../../types.zig");
const fma = @import("../../fma.zig");
const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");

const BlasLong = types.BlasLong;
const HF16 = types.HF16;
const BF16 = types.BF16;

// ---- conversions -----------------------------------------------------------

/// `kernel/x86_64/tobf16.c`.
pub fn toBf16(v: f32) BF16 {
    var bits: u32 = @bitCast(v);
    const top: u16 = switch (bits & 0xff800000) {
        0x00000000 => 0x0000,
        0x80000000 => 0x8000,
        0x7f800000, 0xff800000 => @as(u16, @truncate(bits >> 16)) | @as(u16, if (bits & 0x007fffff != 0) 0x0040 else 0),
        else => blk: {
            bits +%= ((bits >> 16) & 1) + 0x7fff;
            break :blk @truncate(bits >> 16);
        },
    };
    return @fromBackingInt(@intCast(top));
}

/// `kernel/x86_64/bf16to.c`.
pub fn fromBf16(b: BF16) f32 {
    const h: u16 = @backingInt(b);
    const top: u16 = switch (h & 0xff80) {
        0x0000 => 0x0000,
        0x8000 => 0x8000,
        0x7f80, 0xff80 => h | @as(u16, if (h & 0x007f != 0) 0x0040 else 0),
        else => h,
    };
    return @bitCast(@as(u32, top) << 16);
}

/// To bfloat16 from `F`; a double rounds to float first, as upstream does.
fn To(comptime F: type) type {
    return struct {
        fn k(n: BlasLong, in: [*]F, inc_in: BlasLong, out: [*]BF16, inc_out: BlasLong) callconv(.c) void {
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                util.at(BF16, out, i * inc_out).* = toBf16(@floatCast(util.at(F, in, i * inc_in).*));
            }
        }
    };
}

/// From bfloat16 to `F`.
fn From(comptime F: type) type {
    return struct {
        fn k(n: BlasLong, in: [*]BF16, inc_in: BlasLong, out: [*]F, inc_out: BlasLong) callconv(.c) void {
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                util.at(F, out, i * inc_out).* = fromBf16(util.at(BF16, in, i * inc_in).*);
            }
        }
    };
}

pub const sbstobf16 = &To(f32).k;
pub const sbdtobf16 = &To(f64).k;
pub const sbf16tos = &From(f32).k;
pub const dbf16tod = &From(f64).k;

// ---- bscal -----------------------------------------------------------------

/// `kernel/generic/scal.c` built `BFLOAT16`+`BGEMM` (the CMake build; the
/// Makefile one forgets the flags and compiles the `f32` kernel). `flag` is
/// `dummy2`: with 0 a zero `da` zeroes `x`; with 1, which the interface
/// passes, it turns NaN and infinity into NaN. As in the C, the flag-1 test
/// is on the bits of `da` (`bfloat16` is a `uint16_t`, so only +0 counts).
fn bscal(n: BlasLong, _: BlasLong, _: BlasLong, da: BF16, x: [*]BF16, inc_x: BlasLong, _: [*]BF16, _: BlasLong, _: [*]BF16, flag: BlasLong) callconv(.c) i32 {
    if (n <= 0 or inc_x <= 0) return 0;
    const da_f = fromBf16(da);
    var j: BlasLong = 0;
    while (j < n) : (j += 1) {
        const p = util.at(BF16, x, j * inc_x);
        var v: f32 = undefined;
        if (flag == 0) {
            v = if (da_f == 0) 0 else da_f * fromBf16(p.*);
        } else {
            v = fromBf16(p.*);
            v = if (@as(u16, @backingInt(da)) == 0)
                (if (std.math.isFinite(v)) 0 else std.math.nan(f32))
            else
                da_f * v;
        }
        p.* = toBf16(v);
    }
    return 0;
}

pub const bscal_k = &bscal;

// ---- sbdot -----------------------------------------------------------------

/// `x . y` with both operands widened to `f32`, the products and the sum in
/// `f64`, returned as `f32`. (`kernel/x86_64/sbdot.c` widens into `malloc`ed
/// buffers and calls `sdot`; this needs neither.)
fn sbdot(n: BlasLong, x: [*]BF16, inc_x: BlasLong, y: [*]BF16, inc_y: BlasLong) callconv(.c) f32 {
    var sum: f64 = 0;
    var i: BlasLong = 0;
    while (i < n) : (i += 1) {
        sum += @as(f64, fromBf16(util.at(BF16, x, i * inc_x).*)) * fromBf16(util.at(BF16, y, i * inc_y).*);
    }
    return @floatCast(sum);
}

pub const sbdot_portable = &sbdot;

// ---- half-precision gemv ---------------------------------------------------

pub fn wide(comptime T: type, v: T) f32 {
    return if (T == BF16) fromBf16(v) else @floatCast(v);
}

pub fn narrow(comptime T: type, v: f32) T {
    return if (T == BF16) toBf16(v) else @floatCast(v);
}

/// `y = alpha*op(A)*x + beta*y`, `A` column-major `m` by `n`. `E` is the type
/// of `A` and `x`, `Y` that of `y`, `alpha` and `beta`. Each result sums
/// `a*x` in `f32` in index order and rounds once on the store; `beta == 0`
/// never reads `y`; with `m` or `n` 0 the result is `beta*y` (the interface
/// returns before that). `fuse` says whether `acc + a*x` and
/// `alpha*acc + beta*y` contract; the vtable's kernels do not.
pub fn Gemv(comptime E: type, comptime Y: type, comptime trans: bool, comptime fuse: bool) type {
    return struct {
        pub fn k(m: BlasLong, n: BlasLong, alpha: Y, a: [*]E, lda: BlasLong, x: [*]E, inc_x: BlasLong, beta: Y, y: [*]Y, inc_y: BlasLong) callconv(.c) i32 {
            const alpha_f = wide(Y, alpha);
            const beta_f = wide(Y, beta);
            if (trans) {
                // each result a dot down a column of `A`
                var r: BlasLong = 0;
                while (r < n) : (r += 1) {
                    var acc: f32 = 0;
                    var q: BlasLong = 0;
                    while (q < m) : (q += 1) {
                        acc = fma.madd(fuse, f32, wide(E, util.at(E, a, q + r * lda).*), wide(E, util.at(E, x, q * inc_x).*), acc);
                    }
                    store(y, r * inc_y, alpha_f, beta_f, acc);
                }
                return 0;
            }
            // columns outermost over a block of rows, so `A` is read down its
            // columns; each row still sums in column order
            var r0: BlasLong = 0;
            while (r0 < m) : (r0 += row_block) {
                const rows: usize = @intCast(@min(m - r0, row_block));
                var acc: [row_block]f32 = @splat(0);
                var q: BlasLong = 0;
                while (q < n) : (q += 1) {
                    const xv = wide(E, util.at(E, x, q * inc_x).*);
                    const col = util.at(E, a, r0 + q * lda);
                    const c: [*]E = @ptrCast(col);
                    for (acc[0..rows], c[0..rows]) |*s, e| s.* = fma.madd(fuse, f32, wide(E, e), xv, s.*);
                }
                for (acc[0..rows], 0..) |s, i| store(y, (r0 + @as(BlasLong, @intCast(i))) * inc_y, alpha_f, beta_f, s);
            }
            return 0;
        }

        /// `y[at] = alpha*acc + beta*y[at]`, rounded once; `beta == 0` never
        /// reads `y`.
        inline fn store(y: [*]Y, at: BlasLong, alpha_f: f32, beta_f: f32, acc: f32) void {
            const yp = util.at(Y, y, at);
            yp.* = narrow(Y, util.scalar(f32, if (beta_f == 0) alpha_f * acc else fma.madd(fuse, f32, alpha_f, acc, beta_f * wide(Y, yp.*))));
        }
    };
}

/// Rows of the non-transposed gemv summed together: 1 KiB of `f32`
/// accumulators on the stack.
const row_block = 256;

/// The vtable roles of the family: slot prefix, type of `A` and `x`, of `y`.
pub const families = .{
    .{ "shgemv", HF16, f32 },
    .{ "bgemv", BF16, BF16 },
    .{ "sbgemv", BF16, f32 },
};

// ---- installing ------------------------------------------------------------

/// `kernel/setparam-ref.c` fills the bfloat16-named slots of `gotoblas_t`
/// that have no kernel of their own with the `f32` kernels: these level-1
/// slots (`dsbdot` with `dsdot`), which the level-1 interfaces reach.
pub const aliases = [_][2][]const u8{
    .{ "sbamax_k", "samax_k" },   .{ "sbamin_k", "samin_k" },
    .{ "sbmax_k", "smax_k" },     .{ "sbmin_k", "smin_k" },
    .{ "isbamax_k", "isamax_k" }, .{ "isbamin_k", "isamin_k" },
    .{ "isbmax_k", "ismax_k" },   .{ "isbmin_k", "ismin_k" },
    .{ "sbnrm2_k", "snrm2_k" },   .{ "sbasum_k", "sasum_k" },
    .{ "sbsum_k", "ssum_k" },     .{ "sbcopy_k", "scopy_k" },
    .{ "dsbdot_k", "dsdot_k" },   .{ "sbrot_k", "srot_k" },
    .{ "sbrotm_k", "srotm_k" },   .{ "sbaxpy_k", "saxpy_k" },
    .{ "sbscal_k", "sscal_k" },   .{ "sbswap_k", "sswap_k" },
};

/// The same for levels 2 and 3: `sbger_k`, `sbsymv_*`, and every `sbtrsm_*`,
/// `sbtrmm_*` and `sbsymm_*` slot, `sbneg_tcopy` and `sblaswp_ncopy`. No
/// interface calls them; they are null wherever the `f32` slot still is.
pub const aliases_l23: []const [2][]const u8 = blk: {
    @setEvalBranchQuota(100_000);
    var out: []const [2][]const u8 = &.{
        .{ "sbger_k", "sger_k" },             .{ "sbsymv_L", "ssymv_L" },
        .{ "sbsymv_U", "ssymv_U" },           .{ "sbneg_tcopy", "sneg_tcopy" },
        .{ "sblaswp_ncopy", "slaswp_ncopy" },
    };
    for (@typeInfo(Vtable).@"struct".field_names) |f| {
        for ([_][]const u8{ "sbtrsm_", "sbtrmm_", "sbsymm_" }) |p| {
            if (std.mem.startsWith(u8, f, p)) out = out ++ [_][2][]const u8{.{ f, "s" ++ f[2..] }};
        }
    }
    break :blk out;
};

/// Points every alias at what `vt` holds for its `f32` slot. Run once, after
/// all of a vtable's kernels are installed.
pub fn alias(vt: *Vtable) void {
    inline for (aliases) |a| @field(vt, a[0]) = @field(vt, a[1]);
    inline for (aliases_l23) |a| @field(vt, a[0]) = @field(vt, a[1]);
}

/// Points the conversion, `bscal` and `sbdot` slots at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    vt.sbstobf16_k = sbstobf16;
    vt.sbdtobf16_k = sbdtobf16;
    vt.sbf16tos_k = sbf16tos;
    vt.dbf16tod_k = dbf16tod;
    vt.bscal_k = bscal_k;
    vt.sbdot_k = sbdot_portable;
}

/// Points the half-precision gemv slots at the portable kernels.
pub fn installGenericGemv(vt: *Vtable) void {
    inline for (families) |f| {
        @field(vt, f[0] ++ "_n") = &Gemv(f[1], f[2], false, false).k;
        @field(vt, f[0] ++ "_t") = &Gemv(f[1], f[2], true, false).k;
    }
}
