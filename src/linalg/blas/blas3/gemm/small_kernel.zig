//! The small-matrix kernels (`small.zig`): `kernel/generic/gemm_small_matrix_kernel_*.c`
//! and `zgemm_*`, each element summed over the depth in order, on the gemm
//! tile's `@Vector`s; in every core's object. Imports stay in `linalg/blas/`.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../../types.zig");
const coretypes = @import("../../cpu/coretypes.zig");
const gk = @import("kernel.zig");

const BlasLong = types.BlasLong;
const isComplex = types.isComplex;

/// The module whose `Kernels` `core`'s object exports: x86_64's for the
/// cores whose `KERNEL.<core>` names Skylake-X's `[SD]GEMM_SMALL_K_*`
/// (`KERNEL.SKYLAKEX`, included by `COOPERLAKE` and `SAPPHIRERAPIDS`), this
/// one's elsewhere.
pub fn forCore(comptime core: anytype) type {
    if (sysinfo.target.cpu.arch == .x86_64 and @TypeOf(core) == coretypes.X86CoreType) switch (core) {
        .skylakex, .cooperlake, .sapphirerapids => return @import("small_x86_64.zig"),
        else => {},
    };
    return @This();
}

/// The small kernels of `T` under the gemm tile `geo`; named as their slots
/// less the precision letter.
pub fn Kernels(comptime T: type, comptime geo: gk.Geometry) type {
    return if (comptime isComplex(T)) ComplexSmall(T) else RealSmall(T, geo);
}

fn RealSmall(comptime T: type, comptime geo: gk.Geometry) type {
    const cc = types.kernelCc(T);
    // `f16` sums in `f32`, in scalar lanes
    const A = types.Acc(T);
    const vl: usize = if (A == T) geo.vl else 1;
    const mv: usize = geo.mv;
    const nr: usize = geo.nr;
    return struct {
        pub fn gemm_small_kernel_nn(m: BlasLong, n: BlasLong, k: BlasLong, a: [*]T, lda: BlasLong, alpha: T, b: [*]T, ldb: BlasLong, beta: T, c: [*]T, ldc: BlasLong) callconv(cc) i32 {
            return run(false, false, false, m, n, k, a, lda, alpha, b, ldb, beta, c, ldc);
        }
        pub fn gemm_small_kernel_nt(m: BlasLong, n: BlasLong, k: BlasLong, a: [*]T, lda: BlasLong, alpha: T, b: [*]T, ldb: BlasLong, beta: T, c: [*]T, ldc: BlasLong) callconv(cc) i32 {
            return run(false, true, false, m, n, k, a, lda, alpha, b, ldb, beta, c, ldc);
        }
        pub fn gemm_small_kernel_tn(m: BlasLong, n: BlasLong, k: BlasLong, a: [*]T, lda: BlasLong, alpha: T, b: [*]T, ldb: BlasLong, beta: T, c: [*]T, ldc: BlasLong) callconv(cc) i32 {
            return run(true, false, false, m, n, k, a, lda, alpha, b, ldb, beta, c, ldc);
        }
        pub fn gemm_small_kernel_tt(m: BlasLong, n: BlasLong, k: BlasLong, a: [*]T, lda: BlasLong, alpha: T, b: [*]T, ldb: BlasLong, beta: T, c: [*]T, ldc: BlasLong) callconv(cc) i32 {
            return run(true, true, false, m, n, k, a, lda, alpha, b, ldb, beta, c, ldc);
        }
        pub fn gemm_small_kernel_b0_nn(m: BlasLong, n: BlasLong, k: BlasLong, a: [*]T, lda: BlasLong, alpha: T, b: [*]T, ldb: BlasLong, c: [*]T, ldc: BlasLong) callconv(cc) i32 {
            return run(false, false, true, m, n, k, a, lda, alpha, b, ldb, 0, c, ldc);
        }
        pub fn gemm_small_kernel_b0_nt(m: BlasLong, n: BlasLong, k: BlasLong, a: [*]T, lda: BlasLong, alpha: T, b: [*]T, ldb: BlasLong, c: [*]T, ldc: BlasLong) callconv(cc) i32 {
            return run(false, true, true, m, n, k, a, lda, alpha, b, ldb, 0, c, ldc);
        }
        pub fn gemm_small_kernel_b0_tn(m: BlasLong, n: BlasLong, k: BlasLong, a: [*]T, lda: BlasLong, alpha: T, b: [*]T, ldb: BlasLong, c: [*]T, ldc: BlasLong) callconv(cc) i32 {
            return run(true, false, true, m, n, k, a, lda, alpha, b, ldb, 0, c, ldc);
        }
        pub fn gemm_small_kernel_b0_tt(m: BlasLong, n: BlasLong, k: BlasLong, a: [*]T, lda: BlasLong, alpha: T, b: [*]T, ldb: BlasLong, c: [*]T, ldc: BlasLong) callconv(cc) i32 {
            return run(true, true, true, m, n, k, a, lda, alpha, b, ldb, 0, c, ldc);
        }

        fn Vec(comptime n: usize) type {
            return if (n == 1) A else @Vector(n, A);
        }

        /// The operands, unpacked: `op(A)` is `ta`-transposed, `op(B)` `tb`.
        const Ops = struct {
            a: [*]const T,
            lda: usize,
            b: [*]const T,
            ldb: usize,
            k: usize,
            alpha: T,
            beta: T,
            /// `beta` zero: `C` is written, not read.
            b0: bool,
            c: [*]T,
            ldc: usize,
        };

        /// Columns in groups of `nr`, then one at a time; rows in groups of
        /// `mv` vectors, then of one, then single.
        fn run(comptime ta: bool, comptime tb: bool, b0: bool, m: BlasLong, n: BlasLong, k: BlasLong, a: [*]T, lda: BlasLong, alpha: T, b: [*]T, ldb: BlasLong, beta: T, c: [*]T, ldc: BlasLong) i32 {
            const o: Ops = .{ .a = a, .lda = @intCast(lda), .b = b, .ldb = @intCast(ldb), .k = @intCast(k), .alpha = alpha, .beta = beta, .b0 = b0, .c = c, .ldc = @intCast(ldc) };
            const rows: usize = @intCast(m);
            const cols: usize = @intCast(n);
            var j: usize = 0;
            while (j + nr <= cols) : (j += nr) panel(ta, tb, nr, o, rows, j);
            while (j < cols) : (j += 1) panel(ta, tb, 1, o, rows, j);
            return 0;
        }

        fn panel(comptime ta: bool, comptime tb: bool, comptime wn: usize, o: Ops, rows: usize, j0: usize) void {
            var i: usize = 0;
            if (vl > 1) {
                while (i + vl * mv <= rows) : (i += vl * mv) block(ta, tb, vl, mv, wn, o, i, j0);
                while (i + vl <= rows) : (i += vl) block(ta, tb, vl, 1, wn, o, i, j0);
            }
            while (i < rows) : (i += 1) block(ta, tb, 1, 1, wn, o, i, j0);
        }

        /// `lanes` rows of `op(A)` from row `i`, at depth `l`.
        inline fn rowsAt(comptime ta: bool, comptime lanes: usize, o: Ops, i: usize, l: usize) Vec(lanes) {
            if (lanes == 1) return if (ta) o.a[l + i * o.lda] else o.a[i + l * o.lda];
            if (!ta) return o.a[i + l * o.lda ..][0..lanes].*;
            var v: [lanes]A = undefined;
            inline for (0..lanes) |r| v[r] = o.a[l + (i + r) * o.lda];
            return v;
        }

        /// `nv` vectors of `lanes` rows from `r0` by `wn` columns from
        /// `j0`, each element's sum in order over the depth, then `beta*C +
        /// alpha*sum` as the reference's C contracts it.
        fn block(comptime ta: bool, comptime tb: bool, comptime lanes: usize, comptime nv: usize, comptime wn: usize, o: Ops, r0: usize, j0: usize) void {
            @setEvalBranchQuota(100_000);
            const V = Vec(lanes);
            var acc: [wn][nv]V = undefined;
            inline for (0..wn) |j| inline for (0..nv) |v| {
                acc[j][v] = if (lanes == 1) 0 else @splat(0);
            };
            for (0..o.k) |l| {
                var av: [nv]V = undefined;
                inline for (0..nv) |v| av[v] = rowsAt(ta, lanes, o, r0 + v * lanes, l);
                inline for (0..wn) |j| {
                    const y: A = if (tb) o.b[j0 + j + l * o.ldb] else o.b[l + (j0 + j) * o.ldb];
                    const bv: V = if (lanes == 1) y else @splat(y);
                    inline for (0..nv) |v| acc[j][v] = gk.madd(A, lanes, av[v], bv, acc[j][v]);
                }
            }
            const al: V = if (lanes == 1) o.alpha else @splat(o.alpha);
            const be: V = if (lanes == 1) o.beta else @splat(o.beta);
            inline for (0..wn) |j| inline for (0..nv) |v| {
                const p = o.c + r0 + v * lanes + (j0 + j) * o.ldc;
                const prod = al * acc[j][v];
                if (lanes == 1) {
                    const r = if (o.b0) prod else gk.madd(A, 1, p[0], be, prod);
                    p[0] = if (A == T) r else @floatCast(r);
                } else {
                    p[0..lanes].* = if (o.b0) prod else gk.madd(A, lanes, p[0..lanes].*, be, prod);
                }
            };
        }
    };
}

fn ComplexSmall(comptime T: type) type {
    const F = types.Float(T);
    const cc = types.kernelCc(F);
    return struct {
        /// Every op pair, by value: unreachable but under a permit no core
        /// grants, so plain loops.
        fn run(comptime ta: bool, comptime tb: bool, ca: bool, cb: bool, b0: bool, m: BlasLong, n: BlasLong, k: BlasLong, a: [*]F, lda: BlasLong, ar: F, ai: F, b: [*]F, ldb: BlasLong, br: F, bi: F, c: [*]F, ldc: BlasLong) i32 {
            const la: usize = @intCast(lda);
            const lb: usize = @intCast(ldb);
            const lc: usize = @intCast(ldc);
            for (0..@intCast(m)) |i| for (0..@intCast(n)) |j| {
                var re: F = 0;
                var im: F = 0;
                for (0..@intCast(k)) |l| {
                    const pa = 2 * (if (ta) l + i * la else i + l * la);
                    const pb = 2 * (if (tb) j + l * lb else l + j * lb);
                    const x = a[pa];
                    const y = if (ca) -a[pa + 1] else a[pa + 1];
                    const u = b[pb];
                    const v = if (cb) -b[pb + 1] else b[pb + 1];
                    re += gk.madd(F, 1, x, u, -(y * v));
                    im += gk.madd(F, 1, x, v, y * u);
                }
                const p = c + 2 * (i + j * lc);
                if (b0) {
                    p[0] = gk.madd(F, 1, ar, re, -(ai * im));
                    p[1] = gk.madd(F, 1, ar, im, re * ai);
                } else {
                    const t0 = gk.madd(F, 1, br, p[0], -(bi * p[1]));
                    const t1 = gk.madd(F, 1, br, p[1], bi * p[0]);
                    p[0] = gk.madd(F, 1, -ai, im, gk.madd(F, 1, ar, re, t0));
                    p[1] = gk.madd(F, 1, re, ai, gk.madd(F, 1, ar, im, t1));
                }
            };
            return 0;
        }

        fn Full(comptime ops: []const u8) type {
            const ta = ops[0] == 't' or ops[0] == 'c';
            const tb = ops[1] == 't' or ops[1] == 'c';
            const ca = ops[0] == 'r' or ops[0] == 'c';
            const cb = ops[1] == 'r' or ops[1] == 'c';
            return struct {
                fn f(m: BlasLong, n: BlasLong, k: BlasLong, a: [*]F, lda: BlasLong, ar: F, ai: F, b: [*]F, ldb: BlasLong, br: F, bi: F, c: [*]F, ldc: BlasLong) callconv(cc) i32 {
                    return run(ta, tb, ca, cb, false, m, n, k, a, lda, ar, ai, b, ldb, br, bi, c, ldc);
                }
                fn b0(m: BlasLong, n: BlasLong, k: BlasLong, a: [*]F, lda: BlasLong, ar: F, ai: F, b: [*]F, ldb: BlasLong, c: [*]F, ldc: BlasLong) callconv(cc) i32 {
                    return run(ta, tb, ca, cb, true, m, n, k, a, lda, ar, ai, b, ldb, 0, 0, c, ldc);
                }
            };
        }

        pub const gemm_small_kernel_nn = Full("nn").f;
        pub const gemm_small_kernel_nt = Full("nt").f;
        pub const gemm_small_kernel_nr = Full("nr").f;
        pub const gemm_small_kernel_nc = Full("nc").f;
        pub const gemm_small_kernel_tn = Full("tn").f;
        pub const gemm_small_kernel_tt = Full("tt").f;
        pub const gemm_small_kernel_tr = Full("tr").f;
        pub const gemm_small_kernel_tc = Full("tc").f;
        pub const gemm_small_kernel_rn = Full("rn").f;
        pub const gemm_small_kernel_rt = Full("rt").f;
        pub const gemm_small_kernel_rr = Full("rr").f;
        pub const gemm_small_kernel_rc = Full("rc").f;
        pub const gemm_small_kernel_cn = Full("cn").f;
        pub const gemm_small_kernel_ct = Full("ct").f;
        pub const gemm_small_kernel_cr = Full("cr").f;
        pub const gemm_small_kernel_cc = Full("cc").f;
        pub const gemm_small_kernel_b0_nn = Full("nn").b0;
        pub const gemm_small_kernel_b0_nt = Full("nt").b0;
        pub const gemm_small_kernel_b0_nr = Full("nr").b0;
        pub const gemm_small_kernel_b0_nc = Full("nc").b0;
        pub const gemm_small_kernel_b0_tn = Full("tn").b0;
        pub const gemm_small_kernel_b0_tt = Full("tt").b0;
        pub const gemm_small_kernel_b0_tr = Full("tr").b0;
        pub const gemm_small_kernel_b0_tc = Full("tc").b0;
        pub const gemm_small_kernel_b0_rn = Full("rn").b0;
        pub const gemm_small_kernel_b0_rt = Full("rt").b0;
        pub const gemm_small_kernel_b0_rr = Full("rr").b0;
        pub const gemm_small_kernel_b0_rc = Full("rc").b0;
        pub const gemm_small_kernel_b0_cn = Full("cn").b0;
        pub const gemm_small_kernel_b0_ct = Full("ct").b0;
        pub const gemm_small_kernel_b0_cr = Full("cr").b0;
        pub const gemm_small_kernel_b0_cc = Full("cc").b0;
    };
}
