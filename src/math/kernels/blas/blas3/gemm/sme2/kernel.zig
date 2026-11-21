//! The `vortexm4` object's SME2 direct real `gemm` (`sme_sgemm_kernel`,
//! `sme_dgemm_kernel`): `s2_drive` of OpenBLAS `kernel/arm64/sme2_gemm_impl.h`
//! and its micro-kernels (`s2_kernel`, `sme2_gemm_tile.h`), `v0.3.34-794-gdbf2c24c9`,
//! BSD-3-Clause, the OpenBLAS Project, on one job of `job.zig`. The loops and
//! the packing of `A` and of `B` with contiguous columns are Zig; each run of
//! micro-kernels is one `asm` block that commits a caller's lazy ZA save,
//! enters streaming mode and leaves it (`smstart`, `smstop`), so no compiled
//! code runs streaming. Upstream moves `A`'s rows through ZA to transpose
//! them; the values are the same. `B` with contiguous rows is packed online
//! by the first row panel's kernels, as upstream. Tests in `_test/sme2.zig`.

const std = @import("std");
const util = @import("../../../util.zig");
const job = @import("job.zig");

const Slot = util.Slot;

/// The kernel of `slot` (`registry.Port`).
pub fn entry(comptime slot: Slot) util.SlotFn(slot) {
    return switch (slot) {
        .sme_sgemm_kernel => &Kernel(f32).run,
        .sme_dgemm_kernel => &Kernel(f64).run,
        else => @compileError("no SME2 kernel for " ++ @tagName(slot)),
    };
}

/// How a block of `C` comes into ZA (`S2_ZERO`, `S2_LOAD`, `S2_SCALE`).
const Mode = enum(u64) { zero = 0, load = 1, scale = 2 };

/// One run of micro-kernels, read by the `asm` from memory: `count` panels of
/// `B` side by side, the first at `br` (packed) and `bs` (its source when
/// online), `cols` columns left from it; `rows` rows of `C` from `cp`; byte
/// strides `astr` (between `A` panels), `ldb`, `ldc`.
const Args = extern struct {
    kb: u64,
    ar: u64,
    astr: u64,
    br: u64,
    bs: u64,
    ldb: u64,
    cols: u64,
    cp: u64,
    ldc: u64,
    rows: u64,
    mode: u64,
    beta: u64,
    count: u64,
};

fn Kernel(comptime F: type) type {
    const S = job.Shape(F);
    const vl = S.vl;
    const nt = S.nt;
    const nr = S.nr;
    const w2 = S.w2;
    const Bits = if (F == f32) u32 else u64;
    return struct {
        fn run(jb: *const job.Job(F)) callconv(.c) void {
            drive(jb);
        }

        /// `s2_drive`: the Goto nest over `mc` rows, `kc` deep slices and `nc`
        /// columns, then the main kernel over the full rows of tiles, the
        /// half-width one over tails past a vector, the edge one over the rest.
        fn drive(jb: *const job.Job(F)) void {
            const blk = jb.blk;
            const online = !jb.b_cols;
            var i: usize = 0;
            while (i < jb.m) : (i += blk.mc) {
                const mb = @min(jb.m - i, blk.mc);
                var k: usize = 0;
                while (k < jb.k) : (k += blk.kc) {
                    const kb = @min(jb.k - k, blk.kc);
                    if (jb.a_cols) {
                        packACols(mb, kb, jb.a + i + k * jb.lda, jb.lda, jb.alpha, jb.ac);
                    } else {
                        var p0: usize = 0;
                        while (p0 < mb) : (p0 += vl) packARows(@min(mb - p0, vl), kb, jb.a + (i + p0) * jb.lda + k, jb.lda, jb.alpha, jb.ac + p0 * kb);
                    }
                    const mode: Mode = if (k > 0) .load else if (jb.beta == 0) .zero else if (jb.beta == 1) .load else .scale;
                    var j: usize = 0;
                    while (j < jb.n) : (j += blk.nc) {
                        const nb = @min(jb.n - j, blk.nc);
                        const bs = if (jb.b_cols) jb.b + j * jb.ldb + k else jb.b + k * jb.ldb + j;
                        const cb = jb.c + i * jb.ldc + j;
                        columns(jb, mb, kb, nb, bs, cb, mode, online);
                    }
                }
            }
        }

        /// One `nc` block of columns over the packed `mb x kb` block of `A`.
        fn columns(jb: *const job.Job(F), mb: usize, kb: usize, nb: usize, bs: [*]const F, cb: [*]F, mode: Mode, online: bool) void {
            const nmain = nb / nr * nr;
            // tails wider than a vector go to the half-width kernel, the
            // last vector or less to the edge kernel
            const tail = nb - nmain;
            const ecol = if (tail > vl) nmain + (tail - 1) / w2 * w2 + (if ((tail - 1) % w2 >= vl) w2 else 0) else nmain;
            if (!online) {
                var jj: usize = 0;
                while (jj < nmain) : (jj += nr) packBCols(kb, nr, nr, bs + jj * jb.ldb, jb.ldb, jb.bc + jj * kb);
                while (jj < ecol) : (jj += w2) packBCols(kb, @min(nb - jj, w2), w2, bs + jj * jb.ldb, jb.ldb, jb.bc + jj * kb);
                while (jj < nb) : (jj += vl) packBCols(kb, @min(nb - jj, vl), vl, bs + jj * jb.ldb, jb.ldb, jb.bc + jj * kb);
            }
            var args: Args = .{
                .kb = kb,
                .ar = 0,
                .astr = vl * kb * @sizeOf(F),
                .br = 0,
                .bs = 0,
                .ldb = jb.ldb * @sizeOf(F),
                .cols = 0,
                .cp = 0,
                .ldc = jb.ldc * @sizeOf(F),
                .rows = 0,
                .mode = @backingInt(mode),
                .beta = @as(Bits, @bitCast(jb.beta)),
                .count = 0,
            };
            if (nmain > 0) {
                var ii: usize = 0;
                while (ii < mb) : (ii += vl) {
                    set(&args, jb, ii, 0, kb, bs, cb, nmain, @min(mb - ii, vl), nmain / nr);
                    if (online and ii == 0) block(1, nt, true, &args) else block(1, nt, false, &args);
                }
            }
            var jj = nmain;
            while (jj < ecol) : (jj += w2) {
                var ii: usize = 0;
                while (ii < mb) : (ii += 2 * vl) {
                    const rows = @min(mb - ii, 2 * vl);
                    set(&args, jb, ii, jj, kb, bs, cb, @min(nb - jj, w2), rows, 1);
                    switch (rows > vl) {
                        inline else => |two| if (online and ii == 0) block(if (two) 2 else 1, nt / 2, true, &args) else block(if (two) 2 else 1, nt / 2, false, &args),
                    }
                }
            }
            if (ecol < nb) {
                var ii: usize = 0;
                while (ii < mb) : (ii += nt * vl) {
                    const rows = @min(mb - ii, nt * vl);
                    set(&args, jb, ii, ecol, kb, bs, cb, nb - ecol, rows, (nb - ecol + vl - 1) / vl);
                    switch ((rows + vl - 1) / vl) {
                        inline 1...nt => |tm| if (online and ii == 0) block(tm, 1, true, &args) else block(tm, 1, false, &args),
                        else => unreachable,
                    }
                }
            }
        }

        /// The run at row `ii`, column `jj` of the block.
        fn set(args: *Args, jb: *const job.Job(F), ii: usize, jj: usize, kb: usize, bs: [*]const F, cb: [*]F, cols: usize, rows: usize, count: usize) void {
            args.ar = @intFromPtr(jb.ac + ii * kb);
            args.br = @intFromPtr(jb.bc + jj * kb);
            args.bs = @intFromPtr(bs + jj);
            args.cp = @intFromPtr(cb + ii * jb.ldc + jj);
            args.cols = cols;
            args.rows = rows;
            args.count = count;
        }

        /// `A` with contiguous columns into panels of `vl` rows, `alpha`
        /// folded in unless 1 (`s2_pack_a_cols`); rows past `mb` are zero.
        fn packACols(mb: usize, kb: usize, a: [*]const F, ld: usize, alpha: F, ac: [*]F) void {
            const scale = alpha != 1;
            var p0: usize = 0;
            while (p0 < mb) : (p0 += vl) {
                const rows = @min(mb - p0, vl);
                for (0..kb) |l| {
                    const s = a + p0 + l * ld;
                    const d = ac + p0 * kb + l * vl;
                    for (0..rows) |r| d[r] = if (scale) s[r] * alpha else s[r];
                    for (rows..vl) |r| d[r] = 0;
                }
            }
        }

        /// `rows` rows of `A` with contiguous rows into one panel, `alpha`
        /// folded in unless 1 (`s2_pack_a_rows`, `_scaled`); rows past `rows`
        /// are zero.
        fn packARows(rows: usize, kb: usize, a: [*]const F, ld: usize, alpha: F, d: [*]F) void {
            const scale = alpha != 1;
            for (0..rows) |r| {
                const s = a + r * ld;
                for (0..kb) |l| d[l * vl + r] = if (scale) s[l] * alpha else s[l];
            }
            if (rows < vl) for (0..kb) |l| @memset(d[l * vl + rows .. (l + 1) * vl], 0);
        }

        /// A panel of `ncols` (at most `w`) columns of `B` with contiguous
        /// columns, `kb` rows of `w` (`s2_pack_b_cols`); columns past `ncols`
        /// are zero.
        fn packBCols(kb: usize, ncols: usize, w: usize, b: [*]const F, ld: usize, d: [*]F) void {
            for (0..w) |cl| {
                if (cl < ncols) {
                    const s = b + cl * ld;
                    for (0..kb) |l| d[l * w + cl] = s[l];
                } else for (0..kb) |l| {
                    d[l * w + cl] = 0;
                }
            }
        }

        /// The micro-kernels of `tm x tn` tiles over `args`' run (`s2_kernel`,
        /// `S2_DEF_KERNEL`), `online` packing `B` as it goes.
        fn block(comptime tm: usize, comptime tn: usize, comptime online: bool, args: *const Args) void {
            asm volatile (text(tm, tn, online)
                :
                : [args] "{x0}" (args),
                : .{ .x1 = true, .x2 = true, .x3 = true, .x4 = true, .x5 = true, .x6 = true, .x7 = true, .x8 = true, .x9 = true, .x10 = true, .x11 = true, .x12 = true, .x13 = true, .x14 = true, .x15 = true, .x16 = true, .x17 = true, .x19 = true, .x20 = true, .x21 = true, .x22 = true, .x23 = true, .x24 = true, .x25 = true, .z0 = true, .z1 = true, .z2 = true, .z3 = true, .z4 = true, .z5 = true, .z6 = true, .z7 = true, .z8 = true, .z9 = true, .z10 = true, .z11 = true, .z12 = true, .z13 = true, .z14 = true, .z15 = true, .z16 = true, .z17 = true, .z18 = true, .z19 = true, .z20 = true, .z21 = true, .z22 = true, .z23 = true, .z24 = true, .z25 = true, .z26 = true, .z27 = true, .z28 = true, .z29 = true, .z30 = true, .z31 = true, .d8 = true, .d9 = true, .d10 = true, .d11 = true, .d12 = true, .d13 = true, .d14 = true, .d15 = true, .p0 = true, .p1 = true, .p2 = true, .p3 = true, .p4 = true, .p5 = true, .p6 = true, .p7 = true, .p8 = true, .p9 = true, .p10 = true, .p11 = true, .p12 = true, .p13 = true, .p14 = true, .p15 = true, .ffr = true, .nzcv = true, .fpsr = true, .memory = true });
        }

        const es = if (F == f32) "s" else "d";
        const ld1 = if (F == f32) "ld1w" else "ld1d";
        const st1 = if (F == f32) "st1w" else "st1d";
        const width_bytes = vl * @sizeOf(F);

        /// The predicate of column vector `q`: `p0`, all lanes, for the main
        /// kernel's whole panels.
        fn colPred(comptime tn: usize, comptime q: usize) usize {
            return if (tn == nt) 0 else 1 + q;
        }

        fn line(comptime fmt: []const u8, args: anytype) []const u8 {
            return std.fmt.comptimePrint(fmt ++ "\n", args);
        }

        /// The `asm` of `block`. Registers: `x25` the depth, `x1` `A`,
        /// `x2` its panel stride, `x5` `ldb`, `x8` `ldc`, `x9` the rows,
        /// `x10` the mode, `x19` the panels left, `x20` to `x23` the panel's
        /// packed `B`, source `B`, `C` and columns left; `z16` and up `A`,
        /// `z0` to `z15` `B` and `C`, `z31` `beta`.
        fn text(comptime tm: usize, comptime tn: usize, comptime online: bool) []const u8 {
            comptime {
                @setEvalBranchQuota(100_000);
                var t: []const u8 = ".arch_extension sme2\n.arch_extension sme-f64f64\n";
                t = t ++ line("ldp x25, x1, [x0, #0]", .{}) ++ line("ldp x2, x3, [x0, #16]", .{}) ++ line("ldp x4, x5, [x0, #32]", .{}) ++
                    line("ldp x6, x7, [x0, #48]", .{}) ++ line("ldp x8, x9, [x0, #64]", .{}) ++ line("ldp x10, x11, [x0, #80]", .{}) ++ line("ldr x19, [x0, #96]", .{});
                // a caller's ZA under a lazy save goes to its buffer first
                // (AAPCS64's SME support routine `__arm_tpidr2_save`)
                t = t ++ "mrs x13, tpidr2_el0\ncbz x13, 1f\nldr x14, [x13]\nldrh w15, [x13, #8]\ncbz x14, 2f\nmov w12, #0\n3:\ncmp w12, w15\nb.hs 2f\nstr za[w12, 0], [x14]\naddsvl x14, x14, #1\nadd w12, w12, #1\nb 3b\n2:\nmsr tpidr2_el0, xzr\n1:\nsmstart\n";
                t = t ++ line("ptrue p0.{s}", .{es}) ++ line("ptrue pn12.{s}", .{es}) ++ line("dup z31.{s}, {s}11", .{ es, if (F == f32) "w" else "x" });
                t = t ++ "mov x20, x3\nmov x21, x4\nmov x22, x7\nmov x23, x6\n4:\n";
                // the main kernel takes whole panels only, and governing predicates are p0 to p7
                if (tn < nt) for (0..tn) |q| {
                    t = t ++ line("mov x13, #{d}", .{q * vl}) ++ line("whilelt p{d}.{s}, x13, x23", .{ 1 + q, es });
                };
                // `C` in (`s2_c_load`)
                t = t ++ "cbnz x10, 5f\nzero {za}\nb 6f\n5:\n";
                for (0..tm) |p| {
                    t = t ++ rowsHead(p, 100, 200);
                    for (0..tn) |q| t = t ++ line("{s} {{z{d}.{s}}}, p{d}/z, [x16, #{d}, mul vl]", .{ ld1, q, es, colPred(tn, q), q });
                    t = t ++ line("cmp x10, #2", .{}) ++ line("b.ne {d}f", .{120 + p});
                    for (0..tn) |q| t = t ++ line("fmul z{d}.{s}, z{d}.{s}, z31.{s}", .{ q, es, q, es, es });
                    t = t ++ line("{d}:", .{120 + p});
                    for (0..tn) |q| t = t ++ line("mova za{d}h.{s}[w12, 0], p0/m, z{d}.{s}", .{ p * tn + q, es, q, es });
                    t = t ++ rowsEnd(p, 100, 200);
                }
                t = t ++ "6:\nmov x13, x25\nmov x14, x1\nmov x15, x20\nmov x16, x21\n";
                // four depth steps at a time from one panel of `A` (`s2_kernel`)
                if (tm == 1 and !online) {
                    t = t ++ "7:\ncmp x13, #4\nb.lt 8f\n";
                    t = t ++ line("{s} {{z16.{s}-z19.{s}}}, pn12/z, [x14]", .{ ld1, es, es });
                    for (0..4) |u| {
                        const rb = (u % 2) * 8;
                        t = t ++ loadB(tn, rb, false);
                        for (0..tn) |q| t = t ++ line("fmopa za{d}.{s}, p0/m, p0/m, z{d}.{s}, z{d}.{s}", .{ q, es, 16 + u, es, rb + q, es });
                    }
                    t = t ++ "addvl x14, x14, #4\nsub x13, x13, #4\nb 7b\n8:\n";
                }
                t = t ++ "9:\ncbz x13, 10f\n" ++ loadB(tn, 0, online) ++ "mov x17, x14\n";
                for (0..tm) |p| {
                    t = t ++ line("{s} {{z{d}.{s}}}, p0/z, [x17]", .{ ld1, 16 + p, es });
                    if (p + 1 < tm) t = t ++ "add x17, x17, x2\n";
                }
                for (0..tm) |p| for (0..tn) |q| {
                    t = t ++ line("fmopa za{d}.{s}, p0/m, p0/m, z{d}.{s}, z{d}.{s}", .{ p * tn + q, es, 16 + p, es, q, es });
                };
                t = t ++ "addvl x14, x14, #1\nsub x13, x13, #1\nb 9b\n10:\n";
                // `C` out (`s2_c_store`)
                for (0..tm) |p| {
                    t = t ++ rowsHead(p, 140, 220);
                    for (0..tn) |q| t = t ++ line("mova z{d}.{s}, p0/m, za{d}h.{s}[w12, 0]", .{ q, es, p * tn + q, es });
                    for (0..tn) |q| t = t ++ line("{s} {{z{d}.{s}}}, p{d}, [x16, #{d}, mul vl]", .{ st1, q, es, colPred(tn, q), q });
                    t = t ++ rowsEnd(p, 140, 220);
                }
                // the next panel
                t = t ++ "subs x19, x19, #1\nb.eq 11f\n" ++ line("mov x17, #{d}", .{tn * width_bytes}) ++ "madd x20, x25, x17, x20\n" ++
                    line("add x21, x21, #{d}", .{tn * width_bytes}) ++ line("add x22, x22, #{d}", .{tn * width_bytes}) ++ line("sub x23, x23, #{d}", .{tn * vl}) ++ "b 4b\n11:\nsmstop\n";
                return t;
            }
        }

        /// The rows of tile row `p` (`x14`, at most a vector's) from `C`'s
        /// row `p * vl` (`x16`), `w12` counting them; skips to `skip + p`
        /// when there are none. `rowsEnd` closes the loop.
        fn rowsHead(comptime p: usize, comptime top: usize, comptime skip: usize) []const u8 {
            return line("subs x14, x9, #{d}", .{p * vl}) ++ line("b.le {d}f", .{skip + p}) ++ line("cmp x14, #{d}", .{vl}) ++ line("mov x15, #{d}", .{vl}) ++
                "csel x14, x14, x15, lt\n" ++ line("mov x15, #{d}", .{p * vl}) ++ "madd x16, x15, x8, x22\nmov w12, #0\n" ++ line("{d}:", .{top + p});
        }

        fn rowsEnd(comptime p: usize, comptime top: usize, comptime skip: usize) []const u8 {
            return "add x16, x16, x8\nadd w12, w12, #1\ncmp x12, x14\n" ++ line("b.lt {d}b", .{top + p}) ++ line("{d}:", .{skip + p});
        }

        /// One depth step of `B` into `z<rb>` and up: from the source
        /// (`x16`, a row of `ldb`) written on to the packed panel (`x15`)
        /// when `online`, else from the packed panel, in groups of four (or
        /// two) vectors where the panel has them.
        fn loadB(comptime tn: usize, comptime rb: usize, comptime online: bool) []const u8 {
            comptime {
                var t: []const u8 = "";
                if (online) {
                    for (0..tn) |q| t = t ++ line("{s} {{z{d}.{s}}}, p{d}/z, [x16, #{d}, mul vl]", .{ ld1, rb + q, es, colPred(tn, q), q });
                    for (0..tn) |q| t = t ++ line("{s} {{z{d}.{s}}}, p0, [x15, #{d}, mul vl]", .{ st1, rb + q, es, q });
                    t = t ++ "add x16, x16, x5\n";
                } else if (tn >= 4) {
                    for (0..tn / 4) |g| t = t ++ line("{s} {{z{d}.{s}-z{d}.{s}}}, pn12/z, [x15, #{d}, mul vl]", .{ ld1, rb + 4 * g, es, rb + 4 * g + 3, es, 4 * g });
                } else if (tn == 2) {
                    t = t ++ line("{s} {{z{d}.{s}-z{d}.{s}}}, pn12/z, [x15]", .{ ld1, rb, es, rb + 1, es });
                } else {
                    t = t ++ line("{s} {{z{d}.{s}}}, p0/z, [x15]", .{ ld1, rb, es });
                }
                return t ++ line("addvl x15, x15, #{d}", .{tn});
            }
        }
    };
}
