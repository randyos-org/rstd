//! `SBGEMM`, `BGEMM`, `SHGEMM` (`interface/gemm.c` over `level3.c`'s loops) and
//! `SBGEMMT` (`interface/sbgemmt.c`, a half `gemv` per column), on
//! `ext/bf16/kernel.zig`'s kernels; how they round in `gemm/notes.md`.

const std = @import("std");
const sysinfo = @import("builtin");
const blas = @import("../__root__.zig");
const driver = @import("driver.zig");
const interface = @import("interface.zig");
const threaded = @import("threaded.zig");
const threading = @import("../threading.zig");
const test_io = @import("../_test/io.zig");
const types = @import("../types.zig");
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

/// The half-precision gemms: `sbgemm` is bfloat16 in and `f32` out, `bgemm`
/// bfloat16 throughout, `shgemm` `f16` in and `f32` out.
pub const HalfGemm = enum {
    sbgemm,
    bgemm,
    shgemm,

    /// Its blocking parameters' gemm.
    fn gemm(comptime k: HalfGemm) portable.Gemm {
        return switch (k) {
            .sbgemm => .sb,
            .bgemm => .b,
            .shgemm => .sh,
        };
    }

    /// The type of `A` and `B`.
    pub fn In(comptime k: HalfGemm) type {
        return if (k == .shgemm) types.HF16 else types.BF16;
    }

    /// The type of `C`, `alpha` and `beta`.
    pub fn Out(comptime k: HalfGemm) type {
        return if (k == .bgemm) types.BF16 else f32;
    }

    /// The `gemv` of the same types (`sbgemmt`'s kernel).
    fn gemv(comptime k: HalfGemm) blas.HalfGemv {
        return switch (k) {
            .sbgemm => .sbgemv,
            .bgemm => .bgemv,
            .shgemm => .shgemv,
        };
    }
};

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
pub fn gemmHalf(allocator: Allocator, io: Io, comptime kind: HalfGemm, transa: Trans, transb: Trans, m: isize, n: isize, k: isize, alpha: kind.Out(), a: [*]const kind.In(), lda: isize, b: [*]const kind.In(), ldb: isize, beta: kind.Out(), c: [*]kind.Out(), ldc: isize) Error!void {
    try interface.gemmArgs(transa, transb, m, n, k, lda, ldb, ldc);
    const ta = transa.transposed();
    const tb = transb.transposed();
    if (m == 0 or n == 0) return;
    try interface.gemmFits(transa, transb, m, n, k, lda, ldb, ldc, @sizeOf(kind.In()), @sizeOf(kind.Out()));
    // upstream excludes the half types from the per-core tables
    // (`gemm.c:229-239`): every core takes the default rule, as `f16` does
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
pub fn gemmtHalf(io: Io, comptime kind: HalfGemm, uplo: Uplo, transa: Trans, transb: Trans, m: isize, k: isize, alpha: kind.Out(), a: [*]const kind.In(), lda: isize, b: [*]const kind.In(), ldb: isize, beta: kind.Out(), c: [*]kind.Out(), ldc: isize) HalfError!void {
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
pub fn blocking(comptime kind: HalfGemm, vt: *const Vtable) Blocking {
    return .fromSlots(vt, comptime kind.gemm());
}

/// One problem: `kind`'s kernels and the operands.
fn Problem(comptime kind: HalfGemm) type {
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

        fn slot(p: *const Self, comptime name: []const u8) @typeInfo(@FieldType(Vtable, @tagName(kind) ++ name)).optional.child {
            return @field(p.vt, @tagName(kind) ++ name).?;
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
                _ = p.slot("_incopy")(@intCast(min_l), @intCast(min_i), p.a + ls + is * p.lda, @intCast(p.lda), dst);
            } else {
                _ = p.slot("_itcopy")(@intCast(min_l), @intCast(min_i), p.a + is + ls * p.lda, @intCast(p.lda), dst);
            }
        }

        /// Packs columns `jjs..` of `op(B)` over depth `ls..` (`OCOPY_OPERATION`).
        fn packCols(p: *const Self, min_l: usize, min_jj: usize, ls: usize, jjs: usize, dst: [*]I) void {
            if (p.tb) {
                _ = p.slot("_otcopy")(@intCast(min_l), @intCast(min_jj), p.b + jjs + ls * p.ldb, @intCast(p.ldb), dst);
            } else {
                _ = p.slot("_oncopy")(@intCast(min_l), @intCast(min_jj), p.b + ls + jjs * p.ldb, @intCast(p.ldb), dst);
            }
        }

        /// The kernel on packed rows `[ib, ib + mi)` and columns `[jb, jb +
        /// nj)`, the block scaled by `s` first unless null.
        fn block(p: *const Self, ib: usize, mi: usize, jb: usize, nj: usize, kl: usize, s: ?O, pa: [*]I, pb: [*]I) void {
            const at = p.c + ib + jb * p.ldc;
            if (s) |v| _ = p.slot("_beta")(@intCast(mi), @intCast(nj), 0, v, pa, 0, pb, 0, at, @intCast(p.ldc));
            _ = p.slot("_kernel")(@intCast(mi), @intCast(nj), @intCast(kl), p.alpha, pa, pb, at, @intCast(p.ldc));
        }

        /// `level3.c`'s js / ls / is loops over `C`'s `m x n`.
        fn drive(p: *const Self, m: usize, n: usize, work: W) void {
            const blk = p.blk;
            var js: usize = 0;
            while (js < n) : (js += blk.r) {
                const min_j = @min(blk.r, n - js);
                var ls: usize = 0;
                var min_l: usize = 0;
                while (ls < p.k) : (ls += min_l) {
                    min_l = Blocking.split(p.k - ls, blk.q, blk.um);
                    const s = if (ls == 0) p.beta else null;
                    // the first row block packs `B` as it goes
                    var min_i = Blocking.split(m, blk.p, blk.um);
                    const once = min_i == m;
                    p.packRows(min_l, min_i, ls, 0, work.sa);
                    var jjs = js;
                    while (jjs < js + min_j) {
                        const min_jj = blk.chunk(js + min_j - jjs);
                        const pb = if (once) work.sb else work.sb + min_l * (jjs - js);
                        p.packCols(min_l, min_jj, ls, jjs, pb);
                        p.block(0, min_i, jjs, min_jj, min_l, s, work.sa, pb);
                        jjs += min_jj;
                    }
                    var is = min_i;
                    while (is < m) : (is += min_i) {
                        min_i = Blocking.split(m - is, blk.p, blk.um);
                        p.packRows(min_l, min_i, ls, is, work.sa);
                        jjs = js;
                        while (jjs < js + min_j) {
                            const min_jj = blk.chunk(js + min_j - jjs);
                            p.block(is, min_i, jjs, min_jj, min_l, s, work.sa, work.sb + min_l * (jjs - js));
                            jjs += min_jj;
                        }
                    }
                }
            }
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
pub fn gemm(allocator: Allocator, io: Io, comptime kind: HalfGemm, vt: *const Vtable, ta: bool, tb: bool, m: usize, n: usize, k: usize, alpha: kind.Out(), a: [*]const kind.In(), lda: usize, b: [*]const kind.In(), ldb: usize, beta: kind.Out(), c: [*]kind.Out(), ldc: usize, want: u32) Error!void {
    const O = kind.Out();
    const beta_f = toF32(O, beta);
    if (k == 0 or toF32(O, alpha) == 0) {
        if (beta_f != 1) _ = @field(vt, @tagName(kind) ++ "_beta").?(@intCast(m), @intCast(n), 0, beta, @constCast(a), 0, @constCast(b), 0, c, @intCast(ldc));
        return;
    }
    // `SMALL_MATRIX_OPT`, which upstream builds for `sbgemm` alone of these
    if (comptime kind == .sbgemm) {
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
pub fn gemmt(comptime kind: HalfGemm, vt: *const Vtable, upper: bool, ta: bool, tb: bool, m: usize, k: usize, alpha: kind.Out(), a: [*]const kind.In(), lda: usize, b: [*]const kind.In(), ldb: usize, beta: kind.Out(), c: [*]kind.Out(), ldc: usize) void {
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

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../_test/check.zig");
const params = @import("gemm/params.zig");
const zigstd = @import("../../../zigstd.zig");

const kinds = std.enums.values(HalfGemm);

fn halfOf(comptime E: type, v: f32) E {
    return bf16_kernel.narrow(E, v);
}

fn valOf(comptime E: type, v: E) f32 {
    return bf16_kernel.wide(E, v);
}

/// Operands of one test problem, `A` and `B` stored with leading dimensions
/// past their rows, `C` with NaN past its rows (never read nor written).
fn Case(comptime kind: HalfGemm) type {
    const I = kind.In();
    const O = kind.Out();
    return struct {
        const Self = @This();
        a: []I,
        b: []I,
        c: []O,
        lda: usize,
        ldb: usize,
        ldc: usize,

        fn init(al: Allocator, ta: bool, tb: bool, m: usize, n: usize, k: usize, seed: u64) !Self {
            const lda = (if (ta) k else m) + 2;
            const ldb = (if (tb) n else k) + 1;
            const ldc = m + 3;
            const a = try al.alloc(I, lda * @max(1, if (ta) m else k));
            errdefer al.free(a);
            const b = try al.alloc(I, ldb * @max(1, if (tb) k else n));
            errdefer al.free(b);
            const c = try al.alloc(O, ldc * n);
            var prng = std.Random.DefaultPrng.init(seed);
            const r = prng.random();
            for (a) |*v| v.* = halfOf(I, r.float(f32) * 2 - 1);
            for (b) |*v| v.* = halfOf(I, r.float(f32) * 2 - 1);
            for (c, 0..) |*v, i| v.* = halfOf(O, if (i % ldc < m) r.float(f32) * 2 - 1 else std.math.nan(f32));
            return .{ .a = a, .b = b, .c = c, .lda = lda, .ldb = ldb, .ldc = ldc };
        }

        fn deinit(s: Self, al: Allocator) void {
            al.free(s.a);
            al.free(s.b);
            al.free(s.c);
        }

        fn opA(s: Self, ta: bool, i: usize, l: usize) f32 {
            return valOf(I, if (ta) s.a[l + i * s.lda] else s.a[i + l * s.lda]);
        }

        fn opB(s: Self, tb: bool, l: usize, j: usize) f32 {
            return valOf(I, if (tb) s.b[j + l * s.ldb] else s.b[l + j * s.ldb]);
        }
    };
}

/// The kernels' arithmetic for one element, slice by slice: `beta` first
/// (as the beta kernel), then per depth slice the `f32` sum in order and
/// `C + alpha*sum`, rounded into `O`. The products of half values are exact
/// in `f32`, so only the sums and the update round.
fn model(comptime kind: HalfGemm, blk: Blocking, slice: ?usize, cs: Case(kind), ta: bool, tb: bool, k: usize, alpha: f32, beta: f32, i: usize, j: usize, c0: kind.Out()) kind.Out() {
    const O = kind.Out();
    const one = struct {
        fn madd(x: f32, y: f32, z: f32) f32 {
            return util.rounded(f32, util.rounded(f32, x * y) + z);
        }
    };
    if (k == 0 or alpha == 0) return if (beta == 1) c0 else if (beta == 0) halfOf(O, 0) else halfOf(O, beta * valOf(O, c0));
    var cv: O = if (beta == 1) c0 else if (beta == 0) halfOf(O, 0) else halfOf(O, beta * valOf(O, c0));
    var ls: usize = 0;
    while (ls < k) {
        const min_l = if (slice) |q| @min(q, k - ls) else Blocking.split(k - ls, blk.q, blk.um);
        var acc: f32 = 0;
        for (ls..ls + min_l) |l| acc = one.madd(cs.opA(ta, i, l), cs.opB(tb, l, j), acc);
        cv = halfOf(O, one.madd(acc, alpha, valOf(O, cv)));
        ls += min_l;
    }
    return cv;
}

fn sameBits(comptime T: type, x: T, y: T) bool {
    const U = @Int(.unsigned, @bitSizeOf(T));
    const bx: U = if (T == types.BF16) @backingInt(x) else @bitCast(x);
    const by: U = if (T == types.BF16) @backingInt(y) else @bitCast(y);
    return bx == by;
}

/// `kind` on `vt` at `want` workers against `model`, bit for bit, and the
/// rows past `m` untouched; `slice`: the small path's depth slices, where
/// it takes the problem.
fn exact(comptime kind: HalfGemm, vt: *const Vtable, io: Io, slice: ?usize, ta: bool, tb: bool, m: usize, n: usize, k: usize, alpha: f32, beta: f32, want: u32) !void {
    const O = kind.Out();
    const al = testing.allocator;
    const cs = try Case(kind).init(al, ta, tb, m, n, k, m * 131 + n * 17 + k);
    defer cs.deinit(al);
    const orig = try al.dupe(O, cs.c);
    defer al.free(orig);
    const af = valOf(O, halfOf(O, alpha));
    const bt = valOf(O, halfOf(O, beta));
    try gemm(al, io, kind, vt, ta, tb, m, n, k, halfOf(O, alpha), cs.a.ptr, cs.lda, cs.b.ptr, cs.ldb, halfOf(O, beta), cs.c.ptr, cs.ldc, want);
    const blk = blocking(kind, vt);
    for (0..n) |j| for (0..cs.ldc) |i| {
        const at = i + j * cs.ldc;
        const want_v = if (i < m) model(kind, blk, slice, cs, ta, tb, k, af, bt, i, j, orig[at]) else orig[at];
        if (!sameBits(O, want_v, cs.c[at])) {
            check.say("{s} ta={} tb={} m={d} n={d} k={d} alpha={d} beta={d} want={d} C[{d},{d}]: got {d}, want {d}\n", .{ @tagName(kind), ta, tb, m, n, k, alpha, beta, want, i, j, valOf(O, cs.c[at]), valOf(O, want_v) });
            return error.TestExpectedEqual;
        }
    };
}

/// A portable vtable with small blocks, so a few dozen rows and columns
/// cross every block edge.
fn smallBlocks() Vtable {
    var vt: Vtable = .{};
    portable.install(&vt);
    inline for (kinds) |kind| {
        const g = comptime kind.gemm();
        portable.paramPtr(&vt, g, .p).* = 2 * portable.paramPtr(&vt, g, .unroll_m).*;
        portable.paramPtr(&vt, g, .q).* = 5;
        portable.paramPtr(&vt, g, .r).* = 3 * portable.paramPtr(&vt, g, .unroll_n).* + 1;
    }
    return vt;
}

const shapes = [_][3]usize{ .{ 1, 1, 1 }, .{ 3, 2, 7 }, .{ 17, 9, 7 }, .{ 33, 20, 13 }, .{ 8, 40, 26 }, .{ 50, 3, 4 } };
const scalars = [_][2]f32{ .{ 1, 0 }, .{ 1.5, 2 }, .{ -0.5, 1 }, .{ 0, 0.5 }, .{ 2, -1 } };

test "half gemm: the portable kernels give the model's bits, every trans, shape and scalar case" {
    inline for (.{ false, true }) |small| {
        var vt: Vtable = if (small) smallBlocks() else .{};
        if (!small) portable.install(&vt);
        inline for (kinds) |kind| for (shapes) |s| for (scalars) |ab| for ([_]bool{ false, true }) |ta| for ([_]bool{ false, true }) |tb| {
            try exact(kind, &vt, check.io, null, ta, tb, s[0], s[1], s[2], ab[0], ab[1], 1);
        };
    }
}

fn alwaysHalf(_: i32, _: i32, _: BlasLong, _: BlasLong, _: BlasLong, _: f32, _: f32) callconv(.c) i32 {
    return 1;
}

test "half gemm: Cooper Lake's small-matrix sbgemm gives the portable bits in its depth slices" {
    const Cooperlake = @import("gemm/small_x86_64.zig").Cooperlake;
    var vt: Vtable = .{};
    portable.install(&vt);
    vt.sbgemm_small_matrix_permit = &alwaysHalf;
    inline for (.{ false, true }) |b0| inline for (0..2) |xa| inline for (0..2) |xb| {
        const s = comptime portable.GemmSlot.small(b0, xa, xb);
        portable.gemmSlotPtr(&vt, .sb, s).* = Cooperlake.entry(comptime portable.gemmSlotTag(.sb, s));
    };
    // every row block and column tail, and depths past one and two slices
    const deep = [_][3]usize{ .{ 70, 17, 5 }, .{ 15, 8, Cooperlake.slice + 3 }, .{ 33, 2, 2 * Cooperlake.slice + 1 } };
    for (shapes ++ deep) |s| for (scalars) |ab| for ([_]bool{ false, true }) |ta| for ([_]bool{ false, true }) |tb| {
        try exact(.sbgemm, &vt, check.io, Cooperlake.slice, ta, tb, s[0], s[1], s[2], ab[0], ab[1], 1);
    };
    // Cooper Lake's own permit: 8, 16 or 32 columns up to `256^3`
    params.installCore(&vt, .x86_64, .cooperlake, .{});
    try exact(.sbgemm, &vt, check.io, Cooperlake.slice, false, true, 21, 16, 9, 1.5, 2, 1);
    try exact(.sbgemm, &vt, check.io, null, false, true, 21, 15, 9, 1.5, 2, 1);
}

test "half gemm: the same bits at 1 to 9 workers, on every Io" {
    var vt = smallBlocks();
    if (comptime test_io.hosted) {
        var t: Io.Threaded = .init(testing.allocator, .{ .async_limit = .limited(3), .concurrent_limit = .limited(3) });
        defer t.deinit();
        inline for (kinds) |kind| for ([_][3]usize{ .{ 41, 37, 12 }, .{ 70, 9, 11 }, .{ 5, 64, 6 } }) |s| for (1..10) |w| {
            try exact(kind, &vt, t.io(), null, true, false, s[0], s[1], s[2], 1.25, 0.5, @intCast(w));
            try exact(kind, &vt, t.io(), null, false, true, s[0], s[1], s[2], 1.25, 0.5, @intCast(w));
        };
    }
    inline for (kinds) |kind| for (1..5) |w| try exact(kind, &vt, Io.failing, null, false, false, 30, 29, 8, 1, 0, @intCast(w));
}

test "half gemm: a vector tile's lanes give the scalar bits, special values included" {
    // `widen` against `fromBf16` on every class: zeros, denormals, NaNs
    const bits = [_]u16{ 0x0000, 0x8000, 0x0001, 0x8075, 0x3f80, 0xc049, 0x7f80, 0xff80, 0x7f81, 0xffc1, 0x0080, 0x7f7f, 0x1234, 0xabcd, 0x0100, 0x8100 };
    var h: [bits.len]types.BF16 = undefined;
    for (&h, bits) |*d, s| d.* = @fromBackingInt(@intCast(s));
    inline for (.{ 2, 4, 8, 16 }) |n| {
        var at: usize = 0;
        while (at + n <= bits.len) : (at += n) {
            const v: [n]f32 = bf16_kernel.widen(types.BF16, n, h[at..].ptr);
            for (0..n) |l| try testing.expectEqual(@as(u32, @bitCast(bf.fromBf16(h[at + l]))), @as(u32, @bitCast(v[l])));
        }
    }
}

test "half gemm: alpha 0 or k 0 only scales C; beta 0 overwrites NaN; beta 1 leaves C" {
    const al = testing.allocator;
    var vt: Vtable = .{};
    portable.install(&vt);
    inline for (kinds) |kind| {
        const I = kind.In();
        const O = kind.Out();
        const nan_i = halfOf(I, std.math.nan(f32));
        var a = @as([4]I, @splat(nan_i));
        var c = [_]O{ halfOf(O, 1), halfOf(O, 2), halfOf(O, 3), halfOf(O, 4) };
        try gemm(al, Io.failing, kind, &vt, false, false, 2, 2, 2, halfOf(O, 0), &a, 2, &a, 2, halfOf(O, 2), &c, 2, 1);
        try testing.expectEqual(@as(f32, 8), valOf(O, c[3]));
        try gemm(al, Io.failing, kind, &vt, false, false, 2, 2, 0, halfOf(O, 1), &a, 2, &a, 2, halfOf(O, 1), &c, 2, 1);
        try testing.expectEqual(@as(f32, 8), valOf(O, c[3]));
        c = @splat(halfOf(O, std.math.nan(f32)));
        var one = @as([4]I, @splat(halfOf(I, 1)));
        try gemm(al, Io.failing, kind, &vt, true, true, 2, 2, 2, halfOf(O, 1), &one, 2, &one, 2, halfOf(O, 0), &c, 2, 1);
        for (c) |v| try testing.expectEqual(@as(f32, 2), valOf(O, v));
    }
}

test "interface: half gemm reports the first bad argument and does nothing for empty problems" {
    const fa = testing.failing_allocator;
    inline for (kinds) |kind| {
        const I = kind.In();
        const O = kind.Out();
        var a = @as([16]I, @splat(halfOf(I, 1)));
        var c = @as([16]O, @splat(halfOf(O, std.math.nan(f32))));
        const one = halfOf(O, 1);
        const zero = halfOf(O, 0);
        try testing.expectError(error.InvalidTrans, gemmHalf(fa, check.io, kind, .n, .o, -1, -1, -1, one, &a, 0, &a, 0, zero, &c, 0));
        try testing.expectError(error.InvalidM, gemmHalf(fa, check.io, kind, .n, .n, -1, -1, -1, one, &a, 0, &a, 0, zero, &c, 0));
        try testing.expectError(error.InvalidN, gemmHalf(fa, check.io, kind, .n, .n, 0, -1, -1, one, &a, 0, &a, 0, zero, &c, 0));
        try testing.expectError(error.InvalidK, gemmHalf(fa, check.io, kind, .n, .n, 0, 0, -1, one, &a, 0, &a, 0, zero, &c, 0));
        try testing.expectError(error.InvalidLda, gemmHalf(fa, check.io, kind, .n, .n, 3, 2, 4, one, &a, 2, &a, 1, zero, &c, 1));
        try testing.expectError(error.InvalidLda, gemmHalf(fa, check.io, kind, .t, .n, 3, 2, 4, one, &a, 3, &a, 1, zero, &c, 1));
        try testing.expectError(error.InvalidLdb, gemmHalf(fa, check.io, kind, .n, .n, 3, 2, 4, one, &a, 3, &a, 3, zero, &c, 1));
        try testing.expectError(error.InvalidLdb, gemmHalf(fa, check.io, kind, .n, .c, 3, 2, 4, one, &a, 3, &a, 1, zero, &c, 1));
        try testing.expectError(error.InvalidLdc, gemmHalf(fa, check.io, kind, .n, .n, 3, 2, 4, one, &a, 3, &a, 4, zero, &c, 2));
        const big = std.math.maxInt(isize) / 2;
        try testing.expectError(error.InvalidLda, gemmHalf(fa, check.io, kind, .n, .n, 2, 2, 3, one, &a, big, &a, 3, zero, &c, 2));
        // empty: nothing read or written, nothing allocated
        try gemmHalf(fa, check.io, kind, .n, .n, 0, 3, 2, one, &a, 1, &a, 2, zero, &c, 1);
        try gemmHalf(fa, check.io, kind, .n, .n, 3, 0, 2, one, &a, 3, &a, 2, zero, &c, 3);
        try testing.expect(std.math.isNan(valOf(O, c[0])));
        // `r` and `c` are `n` and `t`
        var x = [_]I{ halfOf(I, 1), halfOf(I, 2), halfOf(I, 3), halfOf(I, 4) };
        var y = @as([4]O, @splat(halfOf(O, 0)));
        try gemmHalf(fa, check.io, kind, .c, .r, 2, 2, 2, one, &x, 2, &x, 2, zero, &y, 2);
        // A^T B for A = B = [1 3; 2 4]
        try testing.expectEqual(@as(f32, 5), valOf(O, y[0]));
        try testing.expectEqual(@as(f32, 11), valOf(O, y[1]));
        try testing.expectEqual(@as(f32, 25), valOf(O, y[3]));
    }
}

test "interface: half gemm on a problem past the stack, at the global vtable" {
    const al = testing.allocator;
    inline for (kinds) |kind| {
        const O = kind.Out();
        const m = 70;
        const n = 45;
        const k = 300;
        const cs = try Case(kind).init(al, false, true, m, n, k, 7);
        defer cs.deinit(al);
        const orig = try al.dupe(O, cs.c);
        defer al.free(orig);
        try gemmHalf(al, check.io, kind, .n, .t, m, n, k, halfOf(O, 0.75), cs.a.ptr, @intCast(cs.lda), cs.b.ptr, @intCast(cs.ldb), halfOf(O, -1), cs.c.ptr, @intCast(cs.ldc));
        // `bgemm` rounds to bfloat16 per depth slice: 8 bits each, at
        // most `k` slices; `f32` out by the sums alone
        const rel: f32 = if (O == types.BF16) 0.004 * @as(f32, @floatFromInt(k)) else 2e-7 * @as(f32, @floatFromInt(k));
        for (0..n) |j| for (0..m) |i| {
            var want: f64 = -@as(f64, valOf(O, orig[i + j * cs.ldc]));
            var mag: f64 = @abs(want);
            for (0..k) |l| {
                const p = @as(f64, cs.opA(false, i, l)) * cs.opB(true, l, j) * 0.75;
                want += p;
                mag += @abs(p);
            }
            const got = valOf(O, cs.c[i + j * cs.ldc]);
            if (!(@abs(got - want) <= rel * mag)) {
                check.say("{s} C[{d},{d}]: got {d}, want {d}\n", .{ @tagName(kind), i, j, got, want });
                return error.TestExpectedApproxEqAbs;
            }
        };
    }
}

test "interface: gemmtHalf gives gemvHalf's bits per column, and leaves the other triangle" {
    const al = testing.allocator;
    inline for (kinds) |kind| {
        const I = kind.In();
        const O = kind.Out();
        for ([_]Uplo{ .u, .l }) |uplo| for ([_]Trans{ .n, .t }) |transa| for ([_]Trans{ .n, .c }) |transb| for ([_]usize{ 0, 1, 6 }) |k| for (scalars) |ab| {
            const m = 9;
            const ta = transa.transposed();
            const tb = transb.transposed();
            const cs = try Case(kind).init(al, ta, tb, m, m, k, k * 3 + m);
            defer cs.deinit(al);
            const want = try al.dupe(O, cs.c);
            defer al.free(want);
            const alpha = halfOf(O, ab[0]);
            const beta = halfOf(O, ab[1]);
            try gemmtHalf(check.io, kind, uplo, transa, transb, m, @intCast(k), alpha, cs.a.ptr, @intCast(cs.lda), cs.b.ptr, @intCast(cs.ldb), beta, cs.c.ptr, @intCast(cs.ldc));
            // the reference: one `gemvHalf` per column
            for (0..m) |i| {
                const r0 = if (uplo == .u) 0 else i;
                const len = if (uplo == .u) i + 1 else m - i;
                const bb: [*]I = if (tb) cs.b.ptr + i else cs.b.ptr + i * cs.ldb;
                const incb: isize = @intCast(if (tb) cs.ldb else 1);
                const col = want.ptr + r0 + i * cs.ldc;
                if (k == 0) {
                    if (ab[1] != 1) for (col[0..len]) |*v| {
                        v.* = if (ab[1] == 0) halfOf(O, 0) else halfOf(O, ab[1] * valOf(O, v.*));
                    };
                } else if (ta) {
                    try blas.gemvHalf(check.io, comptime kind.gemv(), .t, @intCast(k), @intCast(len), alpha, cs.a.ptr + r0 * cs.lda, @intCast(cs.lda), bb, incb, beta, col, 1);
                } else {
                    try blas.gemvHalf(check.io, comptime kind.gemv(), .n, @intCast(len), @intCast(k), alpha, cs.a.ptr + r0, @intCast(cs.lda), bb, incb, beta, col, 1);
                }
            }
            for (want, cs.c, 0..) |w, g, at| {
                if (!sameBits(O, w, g) and !(std.math.isNan(valOf(O, w)) and std.math.isNan(valOf(O, g)))) {
                    check.say("{s} {s} {s} {s} k={d} alpha={d} beta={d} C[{d}]: got {d}, want {d}\n", .{ @tagName(kind), @tagName(uplo), @tagName(transa), @tagName(transb), k, ab[0], ab[1], at, valOf(O, g), valOf(O, w) });
                    return error.TestExpectedEqual;
                }
            }
        };
    }
}

test "interface: gemmtHalf reports the first bad argument, and returns for empty problems" {
    inline for (kinds) |kind| {
        const I = kind.In();
        const O = kind.Out();
        var a = @as([16]I, @splat(halfOf(I, std.math.nan(f32))));
        var c = @as([16]O, @splat(halfOf(O, std.math.nan(f32))));
        const one = halfOf(O, 1);
        try testing.expectError(error.InvalidUplo, gemmtHalf(check.io, kind, .v, .o, .n, -1, -1, one, &a, 0, &a, 0, one, &c, 0));
        try testing.expectError(error.InvalidTrans, gemmtHalf(check.io, kind, .u, .o, .n, -1, -1, one, &a, 0, &a, 0, one, &c, 0));
        try testing.expectError(error.InvalidM, gemmtHalf(check.io, kind, .u, .n, .n, -1, -1, one, &a, 0, &a, 0, one, &c, 0));
        try testing.expectError(error.InvalidK, gemmtHalf(check.io, kind, .l, .n, .n, 0, -1, one, &a, 0, &a, 0, one, &c, 0));
        try testing.expectError(error.InvalidLda, gemmtHalf(check.io, kind, .l, .t, .n, 3, 4, one, &a, 3, &a, 4, one, &c, 3));
        try testing.expectError(error.InvalidLdb, gemmtHalf(check.io, kind, .u, .n, .c, 3, 4, one, &a, 3, &a, 2, one, &c, 3));
        try testing.expectError(error.InvalidLdc, gemmtHalf(check.io, kind, .u, .n, .n, 3, 4, one, &a, 3, &a, 4, one, &c, 2));
        try gemmtHalf(check.io, kind, .u, .n, .n, 0, 3, one, &a, 1, &a, 3, halfOf(O, 0), &c, 1);
        try testing.expect(std.math.isNan(valOf(O, c[0])));
        // a zero `alpha` scales the triangle without reading `A` or `B`
        try gemmtHalf(check.io, kind, .l, .n, .n, 2, 2, halfOf(O, 0), &a, 2, &a, 2, halfOf(O, 0), &c, 2);
        try testing.expectEqual(@as(f32, 0), valOf(O, c[0]));
        try testing.expectEqual(@as(f32, 0), valOf(O, c[1]));
        try testing.expect(std.math.isNan(valOf(O, c[2])));
        try testing.expectEqual(@as(f32, 0), valOf(O, c[3]));
    }
}

test "every core keeps the portable half gemm kernels and their unrolls" {
    const Each = struct {
        fn run(vt: *const Vtable) anyerror!void {
            inline for (bf16_kernel.families) |f| {
                const K = bf16_kernel.GemmKernels(f[1], f[2]);
                inline for (bf16_kernel.gemm_slots) |s| try zigstd.testExpectTrue(portable.gemmSlotPtr(vt, f[0], s).* == &@field(K, @tagName(s)));
                try zigstd.testExpectEqual(K.unroll_m, portable.paramPtr(vt, f[0], .unroll_m).*);
                try zigstd.testExpectEqual(K.unroll_n, portable.paramPtr(vt, f[0], .unroll_n).*);
            }
        }
    };
    try check.eachCore(Each.run);
}
