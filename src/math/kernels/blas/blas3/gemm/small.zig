//! The small-matrix gemm path (`SMALL_MATRIX_OPT`, `interface/gemm.c`): where a
//! core's `gemm_small_matrix_permit` allows (`permitFor`; the generic one
//! allows nothing), `gemm` multiplies straight from the caller's matrices.

const types = @import("../../types.zig");
const Vtable = @import("../../Vtable.zig");
const portable = @import("../../portable.zig");
const params = @import("params.zig");
const stest = @import("../../_test/blas3/small.zig");

const BlasLong = types.BlasLong;
const isComplex = types.isComplex;
const skylakex = stest.skylakex;

// ---- the permits ---------------------------------------------------------------

/// The permit rules, by the reference's source: `generic` allows nothing
/// (`kernel/generic/gemm_small_matrix_permit.c`).
pub const Permit = enum {
    generic,
    /// `kernel/x86_64/sgemm_small_kernel_permit_skylakex.c`, `dgemm_...`
    skylakex,
    /// `kernel/arm64/gemm_small_kernel_permit_sve.c`
    sve,
    /// `kernel/power/gemm_small_kernel_permit_power10.c`
    power10,
    /// `kernel/loongarch64/[sd]gemm_small_matrix_permit.c` under `LA264`
    la264,
    /// the same under `LA464`
    la464,
};

/// The rule `core` takes: its `KERNEL.<core>`'s `[SD]GEMM_SMALL_M_PERMIT`,
/// through the files it includes (`KERNEL.COOPERLAKE` includes
/// `KERNEL.SKYLAKEX`, `SAPPHIRERAPIDS` `COOPERLAKE`; `A64FX`, `ARMV9SME`,
/// `CORTEXA510`, `CORTEXA710`, `CORTEXX2` and `NEOVERSEV1` `ARMV8SVE`).
pub fn permitFor(comptime family: params.Family, core: params.CoreOf(family)) Permit {
    return switch (family) {
        .x86_64 => switch (core) {
            .skylakex, .cooperlake, .sapphirerapids => .skylakex,
            else => .generic,
        },
        .arm64 => switch (core) {
            .armv8sve, .a64fx, .armv9sme, .cortexa510, .cortexa710, .cortexx2, .neoversev1 => .sve,
            else => .generic,
        },
        .power64, .power64be => switch (core) {
            .power10 => .power10,
            else => .generic,
        },
        .loongarch64 => switch (core) {
            .la264 => .la264,
            .la464 => .la464,
            else => .generic,
        },
        else => .generic,
    };
}

/// The `sbgemm` permit rules: `generic` allows nothing.
pub const HalfPermit = enum {
    generic,
    /// `kernel/x86_64/sbgemm_small_kernel_permit_cooperlake.c`
    cooperlake,
};

/// The rule `core` takes for `sbgemm`: `KERNEL.COOPERLAKE`'s
/// `SBGEMM_SMALL_M_PERMIT`; `KERNEL.SAPPHIRERAPIDS` empties it, so the
/// generic one, as every other core.
pub fn halfPermitFor(comptime family: params.Family, core: params.CoreOf(family)) HalfPermit {
    return switch (family) {
        .x86_64 => switch (core) {
            .cooperlake => .cooperlake,
            else => .generic,
        },
        else => .generic,
    };
}

/// Whether `rule` lets an `sbgemm` of `m x n` by `k` take the small path:
/// Cooper Lake's up to `256^3`, for 8, 16 or 32 columns ("small matrix
/// kernel works well for N = 8, 16, 32").
pub fn allowsHalf(rule: HalfPermit, m: BlasLong, n: BlasLong, k: BlasLong) bool {
    const mnk = @as(f64, @floatFromInt(m)) * @as(f64, @floatFromInt(n)) * @as(f64, @floatFromInt(k));
    return switch (rule) {
        .generic => false,
        .cooperlake => mnk <= cube(256) and (n == 8 or n == 16 or n == 32),
    };
}

fn HalfPermits(comptime rule: HalfPermit) type {
    return struct {
        fn f(_: i32, _: i32, m: BlasLong, n: BlasLong, k: BlasLong, _: f32, _: f32) callconv(.c) i32 {
            return @intFromBool(allowsHalf(rule, m, n, k));
        }
    };
}

/// `100.0*100.0*100.0` and the like: the cubes the permits compare
/// `M * N * K` against.
fn cube(comptime e: comptime_int) f64 {
    return e * e * e;
}

/// Whether `rule` lets a real `gemm` (`double` for `d`) of `m x n` by `k`
/// take the small path; `ta`, `tb` whether `op` transposes.
pub fn allows(rule: Permit, double: bool, ta: bool, tb: bool, m: BlasLong, n: BlasLong, k: BlasLong) bool {
    const mnk = @as(f64, @floatFromInt(m)) * @as(f64, @floatFromInt(n)) * @as(f64, @floatFromInt(k));
    // the integer products, wide enough not to wrap as upstream's may
    const mn = @as(i128, m) * n;
    return switch (rule) {
        .generic => false,
        // "TN kernel perform not good when C matrix is too big, K is too
        // small"; `s` also "TT ... when K is too small"
        .skylakex => mnk <= cube(100) and if (!ta) true else if (tb) double or k >= 4 else !(mn > 1200 or k < 32),
        // `BLASLONG MNK = M * N * K`, an integer product
        .sve => mn * k <= 64 * 64 * 64,
        // `d` past gcc 11.2's `assemble_pairs` issue; `s` as gcc builds it,
        // where more than one thread caps it at 64^3 and one thread at
        // 100^3: always 64^3, so the path, and the bits, never follow the
        // thread count (`comparisons.md`, 8)
        .power10 => mnk <= if (double) cube(54) else cube(64),
        .la264 => mnk <= if (ta) (if (double) cube(24) else cube(80)) else (if (double) cube(80) else cube(100)),
        .la464 => mnk <= if (ta) (if (double) cube(24) else cube(14)) else (if (double) cube(64) else cube(80)),
    };
}

/// The permit slots' functions for `rule`; complex follows the generic
/// permit on every core, as upstream.
fn Permits(comptime rule: Permit) type {
    return struct {
        fn real(comptime F: type) fn (i32, i32, BlasLong, BlasLong, BlasLong, F, F) callconv(.c) i32 {
            return struct {
                fn f(transa: i32, transb: i32, m: BlasLong, n: BlasLong, k: BlasLong, _: F, _: F) callconv(.c) i32 {
                    return @intFromBool(allows(rule, F == f64, transa != 0, transb != 0, m, n, k));
                }
            }.f;
        }
        fn complex(comptime F: type) fn (i32, i32, BlasLong, BlasLong, BlasLong, F, F, F, F) callconv(.c) i32 {
            return struct {
                fn f(_: i32, _: i32, _: BlasLong, _: BlasLong, _: BlasLong, _: F, _: F, _: F, _: F) callconv(.c) i32 {
                    return 0;
                }
            }.f;
        }
    };
}

/// Points `vt`'s permits at `rule`'s, and `sbgemm`'s at `half`'s.
pub fn installPermits(vt: *Vtable, rule: Permit, half: HalfPermit) void {
    switch (rule) {
        inline else => |r| {
            const P = Permits(r);
            vt.sgemm_small_matrix_permit = &P.real(f32);
            vt.dgemm_small_matrix_permit = &P.real(f64);
            vt.cgemm_small_matrix_permit = &P.complex(f32);
            vt.zgemm_small_matrix_permit = &P.complex(f64);
        },
    }
    switch (half) {
        inline else => |h| vt.sbgemm_small_matrix_permit = &HalfPermits(h).f,
    }
}

// ---- the path --------------------------------------------------------------------

/// The reference's `transa`, `transb` codes: `N` 0, `T` 1, `R` 2, `C` 3.
fn code(trans: bool, conj: bool) i32 {
    return @as(i32, @intFromBool(trans)) | @as(i32, if (conj) 2 else 0);
}

/// `gemm` on the small path if `vt`'s permit allows it: true if done. After
/// `driver.gemm`'s quick returns, so a zero `alpha` or depth reads nothing
/// here either.
pub fn gemm(comptime T: type, vt: *const Vtable, ta: bool, tb: bool, ca: bool, cb: bool, m: usize, n: usize, k: usize, alpha: T, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize) bool {
    if (comptime portable.letter(T) == null) return false;
    const g = comptime portable.Gemm.of(T, false);
    if (comptime !@hasField(Vtable, portable.gemmSlotName(g, .small_matrix_permit))) return false;
    const permit = portable.gemmSlotPtr(vt, g, .small_matrix_permit).* orelse return false;
    const M: BlasLong = @intCast(m);
    const N: BlasLong = @intCast(n);
    const K: BlasLong = @intCast(k);
    const ap: [*]types.Float(T) = @ptrCast(@constCast(a));
    const bp: [*]types.Float(T) = @ptrCast(@constCast(b));
    const cp: [*]types.Float(T) = @ptrCast(c);
    const xa = code(ta, ca);
    const xb = code(tb, cb);
    if (comptime isComplex(T)) {
        if (permit(xa, xb, M, N, K, alpha.re, alpha.im, beta.re, beta.im) == 0) return false;
        const b0 = beta.re == 0 and beta.im == 0;
        switch (xa) {
            inline 0...3 => |ia| switch (xb) {
                inline 0...3 => |ib| {
                    if (b0) {
                        _ = portable.gemmSlotPtr(vt, g, comptime .small(true, ia, ib)).*.?(M, N, K, ap, @intCast(lda), alpha.re, alpha.im, bp, @intCast(ldb), cp, @intCast(ldc));
                    } else {
                        _ = portable.gemmSlotPtr(vt, g, comptime .small(false, ia, ib)).*.?(M, N, K, ap, @intCast(lda), alpha.re, alpha.im, bp, @intCast(ldb), beta.re, beta.im, cp, @intCast(ldc));
                    }
                },
                else => unreachable,
            },
            else => unreachable,
        }
    } else {
        if (permit(xa, xb, M, N, K, alpha, beta) == 0) return false;
        switch (xa) {
            inline 0, 1 => |ia| switch (xb) {
                inline 0, 1 => |ib| {
                    if (beta == 0) {
                        _ = portable.gemmSlotPtr(vt, g, comptime .small(true, ia, ib)).*.?(M, N, K, ap, @intCast(lda), alpha, bp, @intCast(ldb), cp, @intCast(ldc));
                    } else {
                        _ = portable.gemmSlotPtr(vt, g, comptime .small(false, ia, ib)).*.?(M, N, K, ap, @intCast(lda), alpha, bp, @intCast(ldb), beta, cp, @intCast(ldc));
                    }
                },
                else => unreachable,
            },
            else => unreachable,
        }
    }
    return true;
}

/// `sbgemm` on the small path if `vt`'s permit allows it and the core has
/// the kernel: true if done. After `half.gemm`'s quick returns, as `gemm`.
pub fn gemmHalf(vt: *const Vtable, ta: bool, tb: bool, m: usize, n: usize, k: usize, alpha: f32, a: [*]const types.BFloat16, lda: usize, b: [*]const types.BFloat16, ldb: usize, beta: f32, c: [*]f32, ldc: usize) bool {
    const permit = vt.sbgemm_small_matrix_permit orelse return false;
    const M: BlasLong = @intCast(m);
    const N: BlasLong = @intCast(n);
    const K: BlasLong = @intCast(k);
    const ap: [*]types.BFloat16 = @constCast(a);
    const bp: [*]types.BFloat16 = @constCast(b);
    switch (ta) {
        inline else => |xa| switch (tb) {
            inline else => |xb| {
                if (beta == 0) {
                    const f = portable.gemmSlotPtr(vt, .sb, comptime .small(true, @intFromBool(xa), @intFromBool(xb))).* orelse return false;
                    if (permit(code(ta, false), code(tb, false), M, N, K, alpha, beta) == 0) return false;
                    _ = f(M, N, K, ap, @intCast(lda), alpha, bp, @intCast(ldb), c, @intCast(ldc));
                } else {
                    const f = portable.gemmSlotPtr(vt, .sb, comptime .small(false, @intFromBool(xa), @intFromBool(xb))).* orelse return false;
                    if (permit(code(ta, false), code(tb, false), M, N, K, alpha, beta) == 0) return false;
                    _ = f(M, N, K, ap, @intCast(lda), alpha, bp, @intCast(ldb), beta, c, @intCast(ldc));
                }
            },
        },
    }
    return true;
}
