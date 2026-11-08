//! The bfloat16 and half-precision kernels, portable and free of the BLAS root
//! so the per-core objects use them too: the conversions (every core's are
//! `kernel/x86_64/{tobf16,bf16to}.c`, by `Makefile.L1`; denormals flush to a
//! signed zero, NaNs come out quiet), `bscal`, `sbdot`, `?bgemv` and `?bgemm`.

const std = @import("std");
const types = @import("../../types.zig");
const fma = @import("../../fma.zig");
const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");
const gemm_kernels = @import("../../blas3/gemm/kernel.zig");
const vector = @import("../../vector.zig");

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
/// `alpha*acc + beta*y` contract; the vtable's kernels do not. On the x87
/// each step is rounded to `f32` (`madd`).
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
                        acc = madd(fuse, wide(E, util.at(E, a, q + r * lda).*), wide(E, util.at(E, x, q * inc_x).*), acc);
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
                    for (acc[0..rows], c[0..rows]) |*s, e| s.* = madd(fuse, wide(E, e), xv, s.*);
                }
                for (acc[0..rows], 0..) |s, i| store(y, (r0 + @as(BlasLong, @intCast(i))) * inc_y, alpha_f, beta_f, s);
            }
            return 0;
        }

        /// `y[at] = alpha*acc + beta*y[at]`, rounded once; `beta == 0` never
        /// reads `y`.
        inline fn store(y: [*]Y, at: BlasLong, alpha_f: f32, beta_f: f32, acc: f32) void {
            const yp = util.at(Y, y, at);
            yp.* = narrow(Y, util.scalar(f32, if (beta_f == 0) alpha_f * acc else madd(fuse, alpha_f, acc, util.rounded(f32, beta_f * wide(Y, yp.*)))));
        }
    };
}

/// `fma.madd` in `f32`; on the x87 the product and the sum are each rounded
/// to `f32` (`util.rounded`), so the bits do not hang on where values spill.
inline fn madd(comptime fused: bool, a: f32, b: f32, c: f32) f32 {
    if (comptime !fused and util.onX87(f32)) return util.rounded(f32, util.rounded(f32, a * b) + c);
    return fma.madd(fused, f32, a, b, c);
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

// ---- half-precision gemm ---------------------------------------------------

/// The gemm slot prefix of a `families` entry: `sbgemv` is `sbgemm`.
pub fn gemmPrefix(comptime gemv: []const u8) []const u8 {
    return gemv[0 .. gemv.len - 1] ++ "m";
}

/// `kernel/generic/gemmkernel_2x2.c`, `gemm_beta.c`, `gemm_ncopy_2.c` and
/// `gemm_tcopy_2.c` as `Makefile.L3:190-243` builds them for `SHGEMM`,
/// `SBGEMM` and `BGEMM` (bfloat16 through `conversion_macros.h`), over the
/// `f32` gemm tile (`blas3/gemm/kernel.zig`'s layout and `pack`). `I` is the
/// type of `A` and `B`, `O` that of `C` and the scalars. Each element sums
/// `a*b` in `f32` in depth order, unfused, then `C := C + alpha*sum`
/// rounded once into `O`: the same operations in every lane and tile shape,
/// so where a tile falls never changes the bits. On the x87 each step is
/// rounded to `f32` (`madd`).
pub fn GemmKernels(comptime I: type, comptime O: type) type {
    const geo = gemm_kernels.geometry(f32);
    const vl: usize = geo.vl;
    const mr: usize = geo.mr();
    const nr: usize = geo.nr;
    return struct {
        pub const unroll_m: i32 = mr;
        pub const unroll_n: i32 = nr;

        /// `C := beta*C`; a zero `beta` writes zeros without reading `C`.
        /// Upstream's bfloat16 test is on `beta`'s bits (only +0), and its
        /// bfloat16 product truncates; here the value counts, and the
        /// product rounds as `toBf16` does (`issues.md`).
        pub fn beta(m: BlasLong, n: BlasLong, _: BlasLong, b: O, _: [*]I, _: BlasLong, _: [*]I, _: BlasLong, c: [*]O, ldc: BlasLong) callconv(.c) i32 {
            const bf = wide(O, b);
            const h: usize = @intCast(m);
            const ld: usize = @intCast(ldc);
            for (0..@intCast(n)) |j| {
                const col = c[j * ld ..][0..h];
                if (bf == 0) {
                    @memset(col, narrow(O, 0));
                } else {
                    for (col) |*v| v.* = narrow(O, util.scalar(f32, bf * wide(O, v.*)));
                }
            }
            return 0;
        }

        /// Rows of `op(A) = A^T`, `mr` a panel.
        pub fn incopy(m: BlasLong, n: BlasLong, a: [*]I, lda: BlasLong, b: [*]I) callconv(.c) i32 {
            gemm_kernels.pack(I, mr, false, m, n, a, lda, b);
            return 0;
        }
        /// Rows of `op(A) = A`.
        pub fn itcopy(m: BlasLong, n: BlasLong, a: [*]I, lda: BlasLong, b: [*]I) callconv(.c) i32 {
            gemm_kernels.pack(I, mr, true, m, n, a, lda, b);
            return 0;
        }
        /// Columns of `op(B) = B`, `nr` a panel.
        pub fn oncopy(m: BlasLong, n: BlasLong, a: [*]I, lda: BlasLong, b: [*]I) callconv(.c) i32 {
            gemm_kernels.pack(I, nr, false, m, n, a, lda, b);
            return 0;
        }
        /// Columns of `op(B) = B^T`.
        pub fn otcopy(m: BlasLong, n: BlasLong, a: [*]I, lda: BlasLong, b: [*]I) callconv(.c) i32 {
            gemm_kernels.pack(I, nr, true, m, n, a, lda, b);
            return 0;
        }

        /// `C += alpha * A * B` on packed `A` (`bm x bk`) and `B` (`bk x bn`).
        pub fn kernel(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: O, ba: [*]I, bb: [*]I, c: [*]O, ldc: BlasLong) callconv(.c) i32 {
            const m: usize = @intCast(bm);
            const n: usize = @intCast(bn);
            const k: usize = @intCast(bk);
            const ld: usize = @intCast(ldc);
            const al = wide(O, alpha);
            if (vl == 1) {
                scalarPanels(m, n, k, al, ba, bb, c, ld);
                return 0;
            }
            var q0: usize = 0;
            while (q0 + nr <= n) : (q0 += nr) rows(nr, m, k, al, ba, bb + q0 * k, c + q0 * ld, ld);
            inline for (comptime gemm_kernels.tails(nr)) |w| {
                if ((n - q0) & w != 0) {
                    rows(w, m, k, al, ba, bb + q0 * k, c + q0 * ld, ld);
                    q0 += w;
                }
            }
            return 0;
        }

        /// Without vectors: whole tiles unrolled, partial ones shaped by
        /// values, so debug builds of wide scalar tiles stay small.
        fn scalarPanels(m: usize, n: usize, k: usize, al: f32, ba: [*]const I, bb: [*]const I, c: [*]O, ld: usize) void {
            var q0: usize = 0;
            while (q0 < n) {
                const wn = gemm_kernels.panelWidth(n - q0, nr);
                var r0: usize = 0;
                while (r0 < m) {
                    const wm = gemm_kernels.panelWidth(m - r0, mr);
                    if (wm == mr and wn == nr) {
                        tile(mr, nr, k, al, ba + r0 * k, bb + q0 * k, c + r0 + q0 * ld, ld);
                    } else {
                        partTile(wm, wn, k, al, ba + r0 * k, bb + q0 * k, c + r0 + q0 * ld, ld);
                    }
                    r0 += wm;
                }
                q0 += wn;
            }
        }

        fn partTile(wm: usize, wn: usize, k: usize, al: f32, pa: [*]const I, pb: [*]const I, c: [*]O, ld: usize) void {
            var acc: [nr][mr]f32 = @splat(@splat(0));
            var a = pa;
            var b = pb;
            for (0..k) |_| {
                for (0..wn) |j| {
                    const y = wide(I, b[j]);
                    for (0..wm) |i| acc[j][i] = madd(false, wide(I, a[i]), y, acc[j][i]);
                }
                a += wm;
                b += wn;
            }
            for (0..wn) |j| for (0..wm) |i| put(c + i + j * ld, acc[j][i], al);
        }

        fn rows(comptime wn: usize, m: usize, k: usize, al: f32, ba: [*]const I, pb: [*]const I, c: [*]O, ld: usize) void {
            var r0: usize = 0;
            while (r0 + mr <= m) : (r0 += mr) @call(.always_inline, tile, .{ mr, wn, k, al, ba + r0 * k, pb, c + r0, ld });
            inline for (comptime gemm_kernels.tails(mr)) |w| {
                if ((m - r0) & w != 0) {
                    @call(.never_inline, tile, .{ w, wn, k, al, ba + r0 * k, pb, c + r0, ld });
                    r0 += w;
                }
            }
        }

        /// One `wm x wn` tile: `wm / lanes` vectors per column over all of `k`.
        fn tile(comptime wm: usize, comptime wn: usize, k: usize, al: f32, pa: [*]const I, pb: [*]const I, c: [*]O, ld: usize) void {
            @setEvalBranchQuota(100_000);
            const lanes = @min(vl, wm);
            const nv = wm / lanes;
            const V = Lanes(lanes);
            var acc: [wn][nv]V = undefined;
            inline for (0..wn) |j| inline for (0..nv) |v| {
                acc[j][v] = vector.splatTo(V, @as(f32, 0));
            };
            var a = pa;
            var b = pb;
            for (0..k) |_| {
                var av: [nv]V = undefined;
                inline for (0..nv) |v| av[v] = widen(I, lanes, a + v * lanes);
                inline for (0..wn) |j| {
                    const bv = vector.splatTo(V, wide(I, b[j]));
                    inline for (0..nv) |v| acc[j][v] = if (lanes == 1) madd(false, av[v], bv, acc[j][v]) else fma.madd(false, V, av[v], bv, acc[j][v]);
                }
                a += wm;
                b += wn;
            }
            inline for (0..wn) |j| inline for (0..nv) |v| {
                const p = c + j * ld + v * lanes;
                if (lanes == 1) {
                    put(p, acc[j][v], al);
                } else {
                    const s: [lanes]f32 = acc[j][v];
                    inline for (0..lanes) |l| put(p + l, s[l], al);
                }
            };
        }

        /// `C := C + alpha*sum` (`res = res*ALPHA; C = TO_OUTPUT(C + res)`).
        inline fn put(p: [*]O, sum: f32, al: f32) void {
            p[0] = narrow(O, madd(false, sum, al, wide(O, p[0])));
        }
    };
}

fn Lanes(comptime n: usize) type {
    return if (n == 1) f32 else @Vector(n, f32);
}

/// `n` elements at `p` as `f32`, converted as `wide` does each.
pub inline fn widen(comptime T: type, comptime n: usize, p: [*]const T) Lanes(n) {
    if (n == 1) return wide(T, p[0]);
    if (T != BF16) return @floatCast(@as(@Vector(n, T), p[0..n].*));
    const U = @Vector(n, u16);
    const q: [*]const u16 = @ptrCast(p);
    const h: U = q[0..n].*;
    // `fromBf16` by lanes: a zero exponent keeps the sign alone, a NaN quiets
    const e = h & @as(U, @splat(0x7f80));
    const signed = @select(u16, e == @as(U, @splat(0)), h & @as(U, @splat(0x8000)), h);
    const nan = @select(bool, e == @as(U, @splat(0x7f80)), (h & @as(U, @splat(0x007f))) != @as(U, @splat(0)), @as(@Vector(n, bool), @splat(false)));
    const top = @select(u16, nan, signed | @as(U, @splat(0x0040)), signed);
    const w: @Vector(n, u32) = top;
    return @bitCast(w << @as(@Vector(n, u5), @splat(16)));
}

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

/// The same for levels 2 and 3: every `sbtrsm_*`, `sbtrmm_*` and `sbsymm_*`
/// slot, `sbneg_tcopy` and `sblaswp_ncopy`. No interface calls them; they are
/// null wherever the `f32` slot still is.
pub const aliases_l23: []const [2][]const u8 = blk: {
    @setEvalBranchQuota(100_000);
    var out: []const [2][]const u8 = &.{
        .{ "sbneg_tcopy", "sneg_tcopy" }, .{ "sblaswp_ncopy", "slaswp_ncopy" },
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

/// The half-precision gemm slots a kernel set fills, less their prefix.
pub const gemm_slots = [_][]const u8{ "_kernel", "_beta", "_incopy", "_itcopy", "_oncopy", "_otcopy" };

/// Points the half-precision gemm slots at the portable kernels, and their
/// unrolls at the tile.
pub fn installGenericGemm(vt: *Vtable) void {
    inline for (families) |f| {
        const pre = comptime gemmPrefix(f[0]);
        const K = GemmKernels(f[1], f[2]);
        inline for (gemm_slots) |s| @field(vt, pre ++ s) = &@field(K, s[1..]);
    }
    settleGemm(vt);
}

/// Puts back the portable kernels' unrolls where `vt` still holds them: a
/// core's blocking (`params.installCore`) writes `param.h`'s, which describe
/// that core's own kernels, not ported (`blas3/kernels.zig`'s `settle`).
pub fn settleGemm(vt: *Vtable) void {
    inline for (families) |f| {
        const pre = comptime gemmPrefix(f[0]);
        const K = GemmKernels(f[1], f[2]);
        if (@field(vt, pre ++ "_kernel") == &K.kernel) {
            @field(vt, pre ++ "_unroll_m") = K.unroll_m;
            @field(vt, pre ++ "_unroll_n") = K.unroll_n;
            @field(vt, pre ++ "_unroll_mn") = @max(K.unroll_m, K.unroll_n);
        }
    }
}
