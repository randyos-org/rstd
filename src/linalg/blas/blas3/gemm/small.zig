//! The small-matrix gemm path (`SMALL_MATRIX_OPT`, `interface/gemm.c`): where a
//! core's `gemm_small_matrix_permit` allows (`permitFor`; the generic one
//! allows nothing), `gemm` multiplies straight from the caller's matrices.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../../types.zig");
const fma = @import("../../fma.zig");
const Vtable = @import("../../Vtable.zig");
const portable = @import("../../portable.zig");
const gk = @import("kernel.zig");
const params = @import("params.zig");
const small_kernel = @import("small_kernel.zig");

const BlasLong = types.BlasLong;
const isComplex = types.isComplex;

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
pub fn gemmHalf(vt: *const Vtable, ta: bool, tb: bool, m: usize, n: usize, k: usize, alpha: f32, a: [*]const types.BF16, lda: usize, b: [*]const types.BF16, ldb: usize, beta: f32, c: [*]f32, ldc: usize) bool {
    const permit = vt.sbgemm_small_matrix_permit orelse return false;
    const M: BlasLong = @intCast(m);
    const N: BlasLong = @intCast(n);
    const K: BlasLong = @intCast(k);
    const ap: [*]types.BF16 = @constCast(a);
    const bp: [*]types.BF16 = @constCast(b);
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

// ---- tests -----------------------------------------------------------------

const check3 = @import("../../_test/blas3/check.zig");
const check = @import("../../_test/check.zig");

const fillSeq = check.fillSeq;

/// Every kernel against plain arithmetic, around every vector and column
/// group edge; exact on small integers. `C` padding must survive.
fn checkKernels(comptime T: type, comptime geo: gk.Geometry) !void {
    const K = small_kernel.Kernels(T, geo);
    const complex = comptime isComplex(T);
    const F = types.Float(T);
    const cs: usize = if (complex) 2 else 1;
    const al = std.testing.allocator;
    const mr: usize = geo.mr();
    const sizes = [_]usize{ 1, 2, geo.vl + 1, mr, mr + geo.vl + 1, geo.nr, geo.nr + 2 };
    const ld = std.mem.max(usize, &sizes) + 2;
    const buf = try al.alloc(F, 5 * cs * ld * ld);
    defer al.free(buf);
    const a = buf[0 .. cs * ld * ld];
    const b = buf[cs * ld * ld ..][0 .. cs * ld * ld];
    const c = buf[2 * cs * ld * ld ..][0 .. cs * ld * ld];
    const want = buf[3 * cs * ld * ld ..][0 .. cs * ld * ld];
    const got = buf[4 * cs * ld * ld ..][0 .. cs * ld * ld];
    fillSeq(F, a, 1);
    fillSeq(F, b, 2);
    fillSeq(F, c, 3);
    const ops: []const u8 = if (complex) "ntrc" else "nt";
    const lm: BlasLong = @intCast(ld);
    inline for (comptime ops) |oa| inline for (comptime ops) |ob| inline for (.{ false, true }) |b0| {
        const ta = oa == 't' or oa == 'c';
        const tb = ob == 't' or ob == 'c';
        const sa: F = if (oa == 'r' or oa == 'c') -1 else 1;
        const sb: F = if (ob == 'r' or ob == 'c') -1 else 1;
        const name = "gemm_small_kernel_" ++ (if (b0) "b0_" else "") ++ &[_]u8{ oa, ob };
        const f = @field(K, name);
        for (sizes) |m| for (sizes) |n| for ([_]usize{ 1, 5 }) |k| {
            @memcpy(want, c);
            @memcpy(got, c);
            for (0..n) |j| for (0..m) |i| {
                var sr: F = 0;
                var si: F = 0;
                for (0..k) |l| {
                    const ia = cs * (if (ta) l + i * ld else i + l * ld);
                    const ib = cs * (if (tb) j + l * ld else l + j * ld);
                    const x = a[ia];
                    const y = if (complex) sa * a[ia + 1] else 0;
                    const u = b[ib];
                    const v = if (complex) sb * b[ib + 1] else 0;
                    sr += x * u - y * v;
                    si += x * v + y * u;
                }
                // alpha 2 (2 - i), beta 0 or 3 (3 + 2i)
                const p = cs * (i + j * ld);
                if (complex) {
                    const cr = if (b0) 0 else 3 * c[p] - 2 * c[p + 1];
                    const ci = if (b0) 0 else 3 * c[p + 1] + 2 * c[p];
                    want[p] = cr + 2 * sr + si;
                    want[p + 1] = ci + 2 * si - sr;
                } else want[p] = (if (b0) 0 else 3 * c[p]) + 2 * sr;
            };
            if (complex) {
                if (b0) _ = f(@intCast(m), @intCast(n), @intCast(k), a.ptr, lm, 2, -1, b.ptr, lm, got.ptr, lm) else _ = f(@intCast(m), @intCast(n), @intCast(k), a.ptr, lm, 2, -1, b.ptr, lm, 3, 2, got.ptr, lm);
            } else {
                if (b0) _ = f(@intCast(m), @intCast(n), @intCast(k), a.ptr, lm, 2, b.ptr, lm, got.ptr, lm) else _ = f(@intCast(m), @intCast(n), @intCast(k), a.ptr, lm, 2, b.ptr, lm, 3, got.ptr, lm);
            }
            if (!std.mem.eql(F, want, got)) {
                std.debug.print("{s} {s} m={d} n={d} k={d}\n", .{ @typeName(T), name, m, n, k });
                return error.TestUnexpectedResult;
            }
        };
    };
}

test "small kernels match plain arithmetic around every vector and column edge" {
    inline for (.{ f32, f64 }) |T| try checkKernels(T, gk.geometry(T));
    try checkKernels(f64, .{ .vl = 1, .mv = 3, .nr = 5 });
    try checkKernels(f32, .{ .vl = 4, .mv = 2, .nr = 3 });
    inline for (.{ types.Complex32, types.Complex64 }) |T| try checkKernels(T, gk.geometry(T));
}

test "small permits: the reference's thresholds" {
    // Skylake-X: past 100^3 never; `TN` not for a big `C` or a short depth
    try std.testing.expect(allows(.skylakex, true, false, false, 100, 100, 100));
    try std.testing.expect(!allows(.skylakex, true, false, false, 100, 100, 101));
    try std.testing.expect(!allows(.skylakex, true, true, false, 40, 31, 64));
    try std.testing.expect(allows(.skylakex, true, true, false, 40, 30, 32));
    try std.testing.expect(!allows(.skylakex, false, true, true, 4, 4, 3));
    try std.testing.expect(allows(.skylakex, true, true, true, 4, 4, 3));
    try std.testing.expect(allows(.sve, false, true, true, 64, 64, 64));
    try std.testing.expect(!allows(.sve, true, false, false, 64, 64, 65));
    try std.testing.expect(!allows(.power10, true, false, false, 55, 54, 54));
    try std.testing.expect(allows(.power10, false, false, false, 64, 64, 64));
    try std.testing.expect(!allows(.la464, false, true, false, 15, 14, 14));
    try std.testing.expect(allows(.la264, false, true, false, 80, 80, 80));
    try std.testing.expect(!allows(.generic, false, false, false, 1, 1, 1));
    try std.testing.expectEqual(Permit.skylakex, permitFor(.x86_64, .sapphirerapids));
    try std.testing.expectEqual(Permit.generic, permitFor(.x86, .skylakex));
    try std.testing.expectEqual(Permit.sve, permitFor(.arm64, .neoversev1));
    try std.testing.expectEqual(Permit.generic, permitFor(.arm64, .neoversev2));
    // Cooper Lake's `sbgemm`: 8, 16 or 32 columns, up to 256^3
    try std.testing.expect(allowsHalf(.cooperlake, 256, 32, 2048));
    try std.testing.expect(!allowsHalf(.cooperlake, 256, 32, 2049));
    try std.testing.expect(!allowsHalf(.cooperlake, 4, 9, 4));
    try std.testing.expectEqual(HalfPermit.cooperlake, halfPermitFor(.x86_64, .cooperlake));
    try std.testing.expectEqual(HalfPermit.generic, halfPermitFor(.x86_64, .sapphirerapids));
    // a core's blocking brings its permit
    var vt: Vtable = .{};
    params.installCore(&vt, .x86_64, .skylakex, .{});
    try std.testing.expectEqual(@as(i32, 1), vt.dgemm_small_matrix_permit.?(0, 0, 10, 10, 10, 1, 0));
    try std.testing.expectEqual(@as(i32, 0), vt.zgemm_small_matrix_permit.?(0, 0, 1, 1, 1, 1, 0, 0, 0));
    try std.testing.expectEqual(@as(i32, 0), vt.sbgemm_small_matrix_permit.?(0, 0, 8, 8, 8, 1, 0));
    params.installCore(&vt, .x86_64, .haswell, .{});
    try std.testing.expectEqual(@as(i32, 0), vt.dgemm_small_matrix_permit.?(0, 0, 1, 1, 1, 1, 0));
}

fn always(comptime F: type) fn (i32, i32, BlasLong, BlasLong, BlasLong, F, F) callconv(.c) i32 {
    return struct {
        fn f(_: i32, _: i32, _: BlasLong, _: BlasLong, _: BlasLong, _: F, _: F) callconv(.c) i32 {
            return 1;
        }
    }.f;
}

fn alwaysComplex(comptime F: type) fn (i32, i32, BlasLong, BlasLong, BlasLong, F, F, F, F) callconv(.c) i32 {
    return struct {
        fn f(_: i32, _: i32, _: BlasLong, _: BlasLong, _: BlasLong, _: F, _: F, _: F, _: F) callconv(.c) i32 {
            return 1;
        }
    }.f;
}

/// `cases`' vtable: static, as a vtable is large for a rig's stack.
var permissive: Vtable = undefined;

/// `gemm`'s cases with every problem on the small path.
fn cases(vt: *const Vtable) !void {
    permissive = vt.*;
    permissive.sgemm_small_matrix_permit = &always(f32);
    permissive.dgemm_small_matrix_permit = &always(f64);
    permissive.cgemm_small_matrix_permit = &alwaysComplex(f32);
    permissive.zgemm_small_matrix_permit = &alwaysComplex(f64);
    try check3.routineAll(&permissive, .gemm);
}

test "portable small kernels match the wide reference through gemm" {
    var vt: Vtable = .{};
    portable.install(&vt);
    try cases(&vt);
}

test "every linked core's small kernels match the wide reference, where this machine runs them" {
    try check.eachCore(cases);
}

// ---- Skylake-X's kernels, on any x86_64 --------------------------------------------

/// `small_x86_64.zig`'s kernels, which LLVM splits into narrower vectors
/// on a machine without AVX-512: the same lanes, so the same bits.
const skylakex = @import("small_x86_64.zig");

/// Skylake-X's element `(i, j)` of an `m x n x k` problem of `op` (`"nn"`
/// ...), as upstream computes it: the depth in order, fused, except in the
/// dot blocks (`L` lanes, the masked tail, `REDUCE_4` or
/// `_mm512_reduce_add`); then the store. Fused as the kernel's `fma.madd`
/// does: not without an FMA instruction (a baseline x86_64 build, T154).
fn skylakexModel(comptime T: type, comptime op: []const u8, b0: bool, m: usize, n: usize, k: usize, a: []const T, ld: usize, b: []const T, alpha: T, beta: T, c: T, i: usize, j: usize) T {
    const L = skylakex.Skylakex(T).lanes;
    const ta = op[0] == 't';
    const tb = op[1] == 't';
    const at = struct {
        fn f(x: []const T, t: bool, ldx: usize, r: usize, q: usize) T {
            return if (t) x[q + r * ldx] else x[r + q * ldx];
        }
    }.f;
    // the dot rows: all of `tn`'s; `nn`'s past its vectors, for at most
    // half a vector of them over a depth of two vectors or more
    const mm = m % L;
    const start: ?usize = if (ta and !tb) 0 else if (!ta and !tb and mm > 0 and mm <= L / 2 and k >= 2 * L) m - mm else null;
    const fused = struct {
        fn f(x: T, y: T, z: T) T {
            return fma.madd(true, T, x, y, z);
        }
    }.f;
    if (start == null or i < start.?) {
        var acc: T = 0;
        for (0..k) |l| acc = fused(at(a, ta, ld, i, l), at(b, tb, ld, l, j), acc);
        const r = acc * alpha;
        return if (b0) r else fused(c, beta, r);
    }
    var lanes: [64]T = @splat(0);
    var l0: usize = 0;
    while (l0 < k) : (l0 += L) for (0..L) |x| {
        const l = l0 + x;
        lanes[x] = if (l < k) fused(at(a, ta, ld, i, l), at(b, tb, ld, l, j), lanes[x]) else fused(0, 0, lanes[x]);
    };
    const rows4 = (m - start.?) & ~@as(usize, 3);
    const grouped = i - start.? < rows4 or j < (n & ~@as(usize, 3));
    var w: usize = L;
    while (w > 1) : (w /= 2) for (0..w / 2) |x| {
        lanes[x] = if (grouped) lanes[2 * x] + lanes[2 * x + 1] else lanes[x] + lanes[x + w / 2];
    };
    if (grouped) {
        const s = lanes[0] * alpha;
        return if (b0) s else fused(c, beta, s);
    }
    return if (b0) alpha * lanes[0] else fma.madd(fma.fuse, T, alpha, lanes[0], beta * c);
}

/// Every Skylake-X kernel over sizes that reach each of its blocks: equal
/// bits to `skylakexModel` on uneven data, and to plain arithmetic on small
/// integers; `C`'s padding untouched.
fn checkSkylakex(comptime T: type) !void {
    const K = skylakex.Skylakex(T);
    const L = K.lanes;
    const al = std.testing.allocator;
    const ms = [_]usize{ 1, 3, 4, 7, L / 2 + 1, 3 * L / 4, L + 3, 3 * L + L / 2, 7 * L + 3 };
    const ns = [_]usize{ 1, 2, 3, 7, 13 };
    const ms_wide = [_]usize{ 1, 3, 4, 9, 3 * L / 4 - 1 };
    const ns_wide = [_]usize{ 6 * L + L / 2 + 1, 2 * L + 1 };
    const ks = [_]usize{ 1, 2 * L - 1, 2 * L, 2 * L + L / 2 + 1 };
    const ld = 7 * L + 3 + 2;
    const kmax = 2 * L + L / 2 + 1;
    const size = ld * @max(ld, kmax);
    const buf = try al.alloc(T, 5 * size);
    defer al.free(buf);
    const a = buf[0..size];
    const b = buf[size..][0..size];
    const c = buf[2 * size ..][0..size];
    const want = buf[3 * size ..][0..size];
    const got = buf[4 * size ..][0..size];
    var prng = std.Random.DefaultPrng.init(0x5ca1ab1e);
    const rnd = prng.random();
    inline for (.{ false, true }) |exact| {
        for ([_][]T{ a, b, c }) |x| for (x) |*v| {
            v.* = if (exact) @floatFromInt(rnd.intRangeAtMost(i32, -4, 4)) else std.math.ldexp(rnd.float(T) * 2 - 1, rnd.intRangeAtMost(i32, -6, 6));
        };
        const alpha: T = if (exact) 2 else 0.7;
        const beta: T = if (exact) -3 else 1.3;
        inline for (.{ "nn", "nt", "tn", "tt" }) |op| inline for (.{ false, true }) |b0| {
            const f = @field(K, "gemm_small_kernel_" ++ (if (b0) "b0_" else "") ++ op);
            for ([_][]const usize{ &ms, &ms_wide }, [_][]const usize{ &ns, &ns_wide }) |mlist, nlist| for (mlist) |m| for (nlist) |n| for (ks) |k| {
                @memcpy(want, c);
                @memcpy(got, c);
                for (0..n) |j| for (0..m) |i| {
                    const p = i + j * ld;
                    if (exact) {
                        var s: T = 0;
                        for (0..k) |l| s += (if (op[0] == 't') a[l + i * ld] else a[i + l * ld]) * (if (op[1] == 't') b[j + l * ld] else b[l + j * ld]);
                        want[p] = alpha * s + (if (b0) 0 else beta * c[p]);
                    } else want[p] = skylakexModel(T, op, b0, m, n, k, a, ld, b, alpha, beta, c[p], i, j);
                };
                const lm: BlasLong = @intCast(ld);
                if (b0) _ = f(@intCast(m), @intCast(n), @intCast(k), a.ptr, lm, alpha, b.ptr, lm, got.ptr, lm) else _ = f(@intCast(m), @intCast(n), @intCast(k), a.ptr, lm, alpha, b.ptr, lm, beta, got.ptr, lm);
                if (!std.mem.eql(u8, std.mem.sliceAsBytes(want), std.mem.sliceAsBytes(got))) {
                    std.debug.print("{s} {s}{s} exact={} m={d} n={d} k={d}\n", .{ @typeName(T), if (b0) "b0_" else "", op, exact, m, n, k });
                    return error.TestUnexpectedResult;
                }
            };
        };
    }
}

test "Skylake-X small kernels: upstream's blocks, lanes and reductions, on any x86_64" {
    if (comptime sysinfo.target.cpu.arch != .x86_64) return error.SkipZigTest;
    try checkSkylakex(f32);
    try checkSkylakex(f64);
}

test "Skylake-X small kernels match the wide reference through gemm" {
    if (comptime sysinfo.target.cpu.arch != .x86_64) return error.SkipZigTest;
    var vt: Vtable = .{};
    portable.install(&vt);
    inline for (.{ .{ portable.Gemm.s, f32 }, .{ portable.Gemm.d, f64 } }) |p| {
        const K = skylakex.Kernels(p[1], gk.geometry(p[1]));
        inline for (.{ false, true }) |b0| inline for (0..2) |xa| inline for (0..2) |xb| {
            const s = comptime portable.GemmSlot.small(b0, xa, xb);
            portable.gemmSlotPtr(&vt, p[0], s).* = &@field(K, "gemm_" ++ @tagName(s));
        };
    }
    try cases(&vt);
}
