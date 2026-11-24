//! `SBGEMM`, `BGEMM`, `SHGEMM` (`interface/gemm.c` over `level3.c`'s loops) and
//! `SBGEMMT` (`interface/sbgemmt.c`, a half `gemv` per column), on
//! `ext/bf16/kernel.zig`'s kernels; how they round in `gemm/notes.md`.

const std = @import("std");
const blas = @import("../__root__.zig");
const driver = @import("driver.zig");
const loops = @import("loops.zig");
const interface = @import("interface.zig");
const threaded = @import("threaded.zig");
const threading = @import("../threading.zig");
const types = @import("../types.zig");
const half_precision = @import("../half_precision.zig");
const portable = @import("../portable.zig");
const util = @import("../util.zig");
const Vtable = @import("../Vtable.zig");
const bf = @import("../ext/bf16/__root__.zig");
const bf16_kernel = @import("../ext/bf16/kernel.zig");
const small_gemm = @import("gemm/small.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const BlasLong = types.BlasLong;
const Blocking = driver.Blocking;
const Trans = types.Trans;
const Uplo = types.Uplo;

/// What `gemmHalf` returns: `interface.zig`'s level-3 errors.
pub const Error = interface.Error;
/// What `gemmtHalf`, which allocates nothing, returns.
pub const HalfError = interface.ArgError || Io.Cancelable;

/// Which half-precision `gemm` (`half_precision.zig`).
pub const HalfPrecision = half_precision.HalfPrecision;

/// Points the vtable's half-precision gemm slots at the portable kernels.
pub const installGeneric = bf16_kernel.installGenericGemm;

/// A scalar as `f32`: bfloat16 through the vtable's conversion, as the
/// interfaces' `SBF16TOS_K` does.
fn toF32(comptime T: type, v: T) f32 {
    return if (T == f32) v else bf.scalarToF32(v);
}

// ---- the interface ---------------------------------------------------------------

/// `SBGEMM`, `BGEMM`, `SHGEMM` (`interface/gemm.c`): `C := alpha*op(A)*op(B)
/// + beta*C`, `op(A)` being `m x k`, `op(B)` `k x n`; `trans` `r` and `c`
/// are `n` and `t`. Arguments are checked as `gemm`'s; `C` is scaled by
/// `beta` first (a zero `beta` never reads it), and `A` and `B` are not read
/// for a zero `alpha` or `k`. Above `gemm`'s threshold it threads on `io`.
pub fn gemmHalf(allocator: Allocator, io: Io, comptime kind: HalfPrecision, transa: Trans, transb: Trans, m: isize, n: isize, k: isize, alpha: kind.Out(), a: [*]const kind.In(), lda: isize, b: [*]const kind.In(), ldb: isize, beta: kind.Out(), c: [*]kind.Out(), ldc: isize) Error!void {
    try interface.gemmArgs(transa, transb, m, n, k, lda, ldb, ldc);
    const ta = transa.transposed();
    const tb = transb.transposed();
    if (m == 0 or n == 0) return;
    try interface.gemmFits(transa, transb, m, n, k, lda, ldb, ldc, @sizeOf(kind.In()), @sizeOf(kind.Out()));
    // upstream excludes the half types from the per-core tables
    // (`gemm.c:232-242`): every core takes the default rule, as `f16` does
    const w = threaded.workers(.gemm, types.HF16, threading.product(&.{ @intCast(m), @intCast(n), @intCast(k) }));
    try gemm(allocator, io, kind, blas.vtable(), ta, tb, @intCast(m), @intCast(n), @intCast(k), alpha, a, @intCast(lda), b, @intCast(ldb), beta, c, @intCast(ldc), w);
}

/// `SBGEMMT` (`interface/sbgemmt.c`), and the same for `bgemm`'s and
/// `shgemm`'s types, which upstream does not export: `C := alpha*op(A)*op(B)
/// + beta*C` on the `uplo` triangle of the `m x m` `C`, `op(A)` being `m x
/// k`, arguments as `gemmt`'s. As upstream, each column of the triangle is
/// one half `gemv` (`gemvHalf`'s kernel and rounding); the other triangle is
/// never read or written. A zero `alpha` or `k` only scales the triangle by
/// `beta` (a zero `beta` never reads it), where upstream's `gemv` reads `A`
/// and `B` regardless. Serial, as the half `gemv`s are; nothing is allocated.
pub fn gemmtHalf(io: Io, comptime kind: HalfPrecision, uplo: Uplo, transa: Trans, transb: Trans, m: isize, k: isize, alpha: kind.Out(), a: [*]const kind.In(), lda: isize, b: [*]const kind.In(), ldb: isize, beta: kind.Out(), c: [*]kind.Out(), ldc: isize) HalfError!void {
    _ = io;
    try interface.gemmtArgs(uplo, transa, transb, m, k, lda, ldb, ldc);
    const ta = transa.transposed();
    const tb = transb.transposed();
    if (m == 0) return;
    try interface.gemmtFits(transa, transb, m, k, lda, ldb, ldc, @sizeOf(kind.In()), @sizeOf(kind.Out()));
    gemmt(kind, blas.vtable(), uplo == .u, ta, tb, @intCast(m), @intCast(k), alpha, a, @intCast(lda), b, @intCast(ldb), beta, c, @intCast(ldc));
}

// ---- the driver ------------------------------------------------------------------

/// `P`, `Q`, `R` and the unrolls from `kind`'s slots, with `Blocking.of`'s
/// fallbacks for sizes no rule set (mips64 sets none, arm64 no `sh`).
pub fn blocking(comptime kind: HalfPrecision, vt: *const Vtable) Blocking {
    return .fromSlots(vt, comptime kind.gemm());
}

/// One problem: `kind`'s kernels and the operands.
fn Problem(comptime kind: HalfPrecision) type {
    const I = kind.In();
    const O = kind.Out();
    return struct {
        const Self = @This();
        const W = driver.Work(I);

        vt: *const Vtable,
        blk: Blocking,
        ta: bool,
        tb: bool,
        k: usize,
        alpha: O,
        a: [*]I,
        lda: usize,
        b: [*]I,
        ldb: usize,
        /// Null for 1.
        beta: ?O,
        c: [*]O,
        ldc: usize,

        fn slot(p: *const Self, comptime s: portable.GemmSlot) util.SlotFn(portable.gemmSlotTag(kind.gemm(), s)) {
            return portable.gemmSlotPtr(p.vt, comptime kind.gemm(), s).*.?;
        }

        /// The part of `C` at rows `[r0, r0 + m)`, columns `[c0, c0 + n)`.
        fn part(p: *const Self, r0: usize, c0: usize) Self {
            var q = p.*;
            q.a = if (p.ta) p.a + r0 * p.lda else p.a + r0;
            q.b = if (p.tb) p.b + c0 else p.b + c0 * p.ldb;
            q.c = p.c + r0 + c0 * p.ldc;
            return q;
        }

        /// Packs rows `is..` of `op(A)` over depth `ls..` (`ICOPY_OPERATION`).
        fn packRows(p: *const Self, min_l: usize, min_i: usize, ls: usize, is: usize, dst: [*]I) void {
            if (p.ta) {
                _ = p.slot(.incopy)(@intCast(min_l), @intCast(min_i), p.a + ls + is * p.lda, @intCast(p.lda), dst);
            } else {
                _ = p.slot(.itcopy)(@intCast(min_l), @intCast(min_i), p.a + is + ls * p.lda, @intCast(p.lda), dst);
            }
        }

        /// Packs columns `jjs..` of `op(B)` over depth `ls..` (`OCOPY_OPERATION`).
        fn packCols(p: *const Self, min_l: usize, min_jj: usize, ls: usize, jjs: usize, dst: [*]I) void {
            if (p.tb) {
                _ = p.slot(.otcopy)(@intCast(min_l), @intCast(min_jj), p.b + jjs + ls * p.ldb, @intCast(p.ldb), dst);
            } else {
                _ = p.slot(.oncopy)(@intCast(min_l), @intCast(min_jj), p.b + ls + jjs * p.ldb, @intCast(p.ldb), dst);
            }
        }

        /// The kernel on packed rows `[ib, ib + mi)` and columns `[jb, jb +
        /// nj)`, the block scaled by `s` first unless null.
        fn block(p: *const Self, ib: usize, mi: usize, jb: usize, nj: usize, kl: usize, s: ?O, pa: [*]I, pb: [*]I) void {
            const at = p.c + ib + jb * p.ldc;
            if (s) |v| _ = p.slot(.beta)(@intCast(mi), @intCast(nj), 0, v, pa, 0, pb, 0, at, @intCast(p.ldc));
            _ = p.slot(.kernel)(@intCast(mi), @intCast(nj), @intCast(kl), p.alpha, pa, pb, at, @intCast(p.ldc));
        }

        /// `level3.c`'s js / ls / is loops over `C`'s `m x n` (`loops.zig`).
        fn drive(p: *const Self, m: usize, n: usize, work: W) void {
            const Ctx = struct {
                p: *const Self,
                m: usize,
                work: W,

                pub fn depth(s: *const @This(), rest: usize) usize {
                    return Blocking.split(rest, s.p.blk.q, s.p.blk.um);
                }
                pub fn chunk(s: *const @This(), rest: usize) usize {
                    return s.p.blk.chunk(rest);
                }
                pub fn slice(s: *const @This(), js: usize, min_j: usize, ls: usize, min_l: usize) void {
                    loops.slice(s, s.p.blk, 0, s.m, js, min_j, ls, min_l, if (ls == 0) s.p.beta else null, s.work.sa, s.work.sb);
                }
                pub fn packRows(s: *const @This(), min_l: usize, min_i: usize, ls: usize, is: usize, dst: [*]I) void {
                    s.p.packRows(min_l, min_i, ls, is, dst);
                }
                pub fn packCols(s: *const @This(), min_l: usize, min_jj: usize, ls: usize, jjs: usize, dst: [*]I) void {
                    s.p.packCols(min_l, min_jj, ls, jjs, dst);
                }
                pub fn block(s: *const @This(), ib: usize, mi: usize, jb: usize, nj: usize, kl: usize, b: ?O, pa: [*]I, pb: [*]I) void {
                    s.p.block(ib, mi, jb, nj, kl, b, pa, pb);
                }
            };
            const ctx: Ctx = .{ .p = p, .m = m, .work = work };
            loops.columns(p.blk, 0, n, p.k, &ctx);
        }

        /// Up to `want` workers, each on a range of columns (rows when `C`
        /// is taller than wide), starting on a panel.
        fn run(p: *const Self, allocator: Allocator, io: Io, m: usize, n: usize, want: u32) driver.Error!void {
            if (want <= 1) {
                const lens = W.lens(p.blk, m, n, p.k, false);
                const len = W.total(lens);
                var stack: [driver.stack_bytes]u8 align(W.alignment) = undefined;
                const sc = try util.Scratch(I, W.alignment, W.alignment).init(allocator, &stack, len);
                defer sc.deinit(allocator);
                return p.drive(m, n, W.at(sc.ptr, lens));
            }
            const Parts = struct {
                p: *const Self,
                m: usize,
                n: usize,
                cols: bool,

                fn dims(s: *const @This(), ch: threading.Chunk) [2]usize {
                    return if (s.cols) .{ s.m, ch.end - ch.start } else .{ ch.end - ch.start, s.n };
                }

                fn lens(s: *const @This(), ch: threading.Chunk) [3]usize {
                    const d = s.dims(ch);
                    return W.lens(s.p.blk, d[0], d[1], s.p.k, false);
                }

                fn len(s: *const @This(), ch: threading.Chunk) usize {
                    return W.total(s.lens(ch));
                }

                fn one(s: *const @This(), ch: threading.Chunk, buf: []I) void {
                    const d = s.dims(ch);
                    const q = if (s.cols) s.p.part(0, ch.start) else s.p.part(ch.start, 0);
                    q.drive(d[0], d[1], W.at(buf.ptr, s.lens(ch)));
                }
            };
            const cols = n >= m;
            const plan: threading.Plan = .{
                .chunks = .init(if (cols) n else m, want, .even, if (cols) p.blk.un else p.blk.um),
                .workers = want,
            };
            const parts: Parts = .{ .p = p, .m = m, .n = n, .cols = cols };
            try threading.forEachPart(allocator, io, &plan, I, W.alignment, &parts, Parts.len, Parts.one);
        }
    };
}

/// `C := alpha*op(A)*op(B) + beta*C` with `kind`'s kernels in `vt`; `m`,
/// `n` positive, arguments checked by the caller. Up to `want` workers.
pub fn gemm(allocator: Allocator, io: Io, comptime kind: HalfPrecision, vt: *const Vtable, ta: bool, tb: bool, m: usize, n: usize, k: usize, alpha: kind.Out(), a: [*]const kind.In(), lda: usize, b: [*]const kind.In(), ldb: usize, beta: kind.Out(), c: [*]kind.Out(), ldc: usize, want: u32) Error!void {
    const O = kind.Out();
    const beta_f = toF32(O, beta);
    if (k == 0 or toF32(O, alpha) == 0) {
        if (beta_f != 1) _ = portable.gemmSlotPtr(vt, comptime kind.gemm(), .beta).*.?(@intCast(m), @intCast(n), 0, beta, @constCast(a), 0, @constCast(b), 0, c, @intCast(ldc));
        return;
    }
    // `SMALL_MATRIX_OPT`, which upstream builds for `sbgemm` alone of these
    if (comptime kind == .bf16_to_f32) {
        if (small_gemm.gemmHalf(vt, ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)) return;
    }
    const p: Problem(kind) = .{
        .vt = vt,
        .blk = blocking(kind, vt),
        .ta = ta,
        .tb = tb,
        .k = k,
        .alpha = alpha,
        .a = @constCast(a),
        .lda = lda,
        .b = @constCast(b),
        .ldb = ldb,
        .beta = if (beta_f == 1) null else beta,
        .c = c,
        .ldc = ldc,
    };
    try p.run(allocator, io, m, n, want);
}

/// `sbgemmt.c`'s loop: column `i` of the triangle is `gemv_n` (or `gemv_t`
/// for `A^T`) of its rows of `op(A)` and column `i` of `op(B)`. `m` positive.
pub fn gemmt(comptime kind: HalfPrecision, vt: *const Vtable, upper: bool, ta: bool, tb: bool, m: usize, k: usize, alpha: kind.Out(), a: [*]const kind.In(), lda: usize, b: [*]const kind.In(), ldb: usize, beta: kind.Out(), c: [*]kind.Out(), ldc: usize) void {
    const O = kind.Out();
    const pa: [*]kind.In() = @constCast(a);
    const pb: [*]kind.In() = @constCast(b);
    const incb: BlasLong = @intCast(if (tb) ldb else 1);
    for (0..m) |i| {
        // the rows of column `i` in the triangle
        const r0 = if (upper) 0 else i;
        const len = if (upper) i + 1 else m - i;
        const cc = c + r0 + i * ldc;
        if (k == 0 or toF32(O, alpha) == 0) {
            if (toF32(O, beta) != 1) {
                if (comptime O == f32) {
                    _ = vt.sscal_k.?(@intCast(len), 0, 0, beta, cc, 1, cc, 0, cc, 0);
                } else {
                    _ = vt.bscal_k.?(@intCast(len), 0, 0, beta, cc, 1, cc, 0, cc, 0);
                }
            }
            continue;
        }
        const bb = if (tb) pb + i else pb + i * ldb;
        if (ta) {
            _ = portable.halfGemvPtr(vt, comptime kind.gemm(), true).*.?(@intCast(k), @intCast(len), alpha, pa + r0 * lda, @intCast(lda), bb, incb, beta, cc, 1);
        } else {
            _ = portable.halfGemvPtr(vt, comptime kind.gemm(), false).*.?(@intCast(len), @intCast(k), alpha, pa + r0, @intCast(lda), bb, incb, beta, cc, 1);
        }
    }
}
