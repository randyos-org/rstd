//! The Zig gemm kernels for real `T`: the packing copies, `beta`, and the
//! register-tile kernel `C += alpha * A * B` on packed panels (the slots
//! `<p>gemm_{incopy,itcopy,oncopy,otcopy,beta,kernel}`). Compiled into the
//! main library as the portable kernels (`installGeneric`), and into every
//! per-core object (`Registry`), where the core's features set the vector
//! width and register count, so one source serves every core.
//!
//! Packed layout, shared by the copies and the kernel: `n / u` panels of
//! width `u`, then one of each power of two below `u` set in `n % u`. Panel
//! `p` at column `c0` with width `w` sits at `m * c0`, as `m` groups of `w`.
//! For a power-of-two `u` this is OpenBLAS's `kernel/generic` layout.
//!
//! The tile is `mr x nr`: `mv` vectors of `vl` lanes by `nr` columns, the
//! most accumulators the registers hold beside the `mv` loads of A and a
//! broadcast of B (`geometry`). Imports stay inside `linalg/blas/`: the
//! per-core objects compile this file on their own.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../../types.zig");
const fma = @import("../../fma.zig");
const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");

const BlasLong = types.BlasLong;
const XDouble = types.XDouble;

/// A kernel set's tile: `mv` vectors of `vl` lanes (`mr` rows) by `nr`
/// columns; `mr` and `nr` are also the copies' unrolls.
pub const Geometry = struct {
    vl: u16,
    mv: u16,
    nr: u16,

    pub fn mr(g: Geometry) u16 {
        return g.vl * g.mv;
    }
};

/// Bits of the widest vector of `T` the target computes on natively; 0 for
/// none. Not `std.simd.suggestVectorLength`: that counts SVE as 256-bit
/// fixed vectors (LLVM splits them into NEON pairs) and leaves AVX at 128.
fn vectorBits(comptime T: type, comptime cpu: std.Target.Cpu) u16 {
    const double = T == f64;
    if (T != f32 and T != f64) return 0;
    return switch (cpu.arch) {
        .x86_64, .x86 => if (cpu.has(.x86, .avx512f) and sysinfo.zig_backend != .stage2_x86_64)
            512
        else if (cpu.has(.x86, .avx))
            256
        else if (cpu.has(.x86, if (double) .sse2 else .sse))
            128
        else
            0,
        .aarch64, .aarch64_be => if (cpu.has(.aarch64, .neon)) 128 else 0,
        // ARMv7 NEON has no double lanes
        .arm, .armeb, .thumb, .thumbeb => if (!double and cpu.has(.arm, .neon)) 128 else 0,
        // AltiVec has no double lanes; VSX has
        .powerpc, .powerpcle, .powerpc64, .powerpc64le => if (cpu.has(.powerpc, .vsx) or (!double and cpu.has(.powerpc, .altivec))) 128 else 0,
        .riscv32, .riscv32be, .riscv64, .riscv64be => if (!cpu.has(.riscv, .v))
            0
        else if (cpu.has(.riscv, .zvl512b))
            512
        else if (cpu.has(.riscv, .zvl256b))
            256
        else
            128,
        .loongarch32, .loongarch64 => if (cpu.has(.loongarch, .lasx)) 256 else if (cpu.has(.loongarch, .lsx)) 128 else 0,
        .mips, .mipsel, .mips64, .mips64el => if (cpu.has(.mips, .msa)) 128 else 0,
        // z13 vectors hold doubles; singles need the z14 enhancements
        .s390x => if (cpu.has(.s390x, .vector) and (double or cpu.has(.s390x, .vector_enhancements_1))) 128 else 0,
        .wasm32, .wasm64 => if (cpu.has(.wasm, .simd128)) 128 else 0,
        else => 0,
    };
}

/// Registers that hold one such vector (or one scalar where `bits` is 0).
fn registers(comptime bits: u16, comptime cpu: std.Target.Cpu) u16 {
    return switch (cpu.arch) {
        .x86_64 => if (bits == 512) 32 else 16,
        .x86 => 8,
        .arm, .armeb, .thumb, .thumbeb => if (bits == 0 and cpu.has(.arm, .d32)) 32 else 16,
        .powerpc, .powerpcle, .powerpc64, .powerpc64le => if (cpu.has(.powerpc, .vsx)) 64 else 32,
        .s390x => if (bits == 0) 16 else 32,
        .wasm32, .wasm64 => 16,
        .mips, .mipsel, .sparc => 16,
        else => 32,
    };
}

/// Whether LLVM multiplies by a lane of a register (AArch64 `fmla .s[i]`,
/// NEON `vmla` by scalar): B then comes as whole vectors, `nr / vl` registers,
/// rather than one broadcast at a time.
fn byLane(comptime cpu: std.Target.Cpu, comptime bits: u16) bool {
    return bits > 0 and switch (cpu.arch) {
        .aarch64, .aarch64_be, .arm, .armeb, .thumb, .thumbeb => true,
        else => false,
    };
}

/// The tile for `T` on `cpu`: the most accumulators (`mv * nr`) that fit
/// beside `mv` A vectors, B (one broadcast, or `nr / vl` vectors where
/// `byLane`) and one spare register, with no more loads per step than three
/// quarters of the multiply-adds (relaxed where nothing meets that); ties
/// go to fewer loads, then to wider `nr`. Without the spare, LLVM spills
/// accumulators in the loop (`notes.md`).
pub fn geometryFor(comptime T: type, comptime cpu: std.Target.Cpu) Geometry {
    const bits = vectorBits(T, cpu);
    const vl: u16 = if (bits == 0) 1 else bits / @bitSizeOf(T);
    const regs = registers(bits, cpu);
    @setEvalBranchQuota(10_000);
    inline for (.{ true, false }) |strict| {
        var best: ?Geometry = null;
        for (1..5) |mv| for (1..17) |nr| {
            const acc = mv * nr;
            const b_regs = if (byLane(cpu, bits)) (nr + vl - 1) / vl else 1;
            if (acc + mv + b_regs + 1 > regs) continue;
            if (strict and 4 * (mv + nr) > 3 * acc) continue;
            const g: Geometry = .{ .vl = vl, .mv = mv, .nr = nr };
            if (best) |b| {
                const bacc = @as(usize, b.mv) * b.nr;
                if (acc < bacc) continue;
                if (acc == bacc and mv + nr > b.mv + b.nr) continue;
                if (acc == bacc and mv + nr == b.mv + b.nr and nr < b.nr) continue;
            }
            best = g;
        };
        if (best) |b| return b;
    }
    unreachable;
}

/// The compile target's tile.
pub fn geometry(comptime T: type) Geometry {
    return comptime geometryFor(T, sysinfo.cpu);
}

/// The powers of two below `u`, largest first: the remainder panels' widths.
fn tails(comptime u: usize) []const usize {
    comptime {
        var out: []const usize = &.{};
        if (u < 2) return out;
        var w = std.math.floorPowerOfTwo(usize, u - 1);
        while (w >= 1) : (w /= 2) out = out ++ [_]usize{w};
        return out;
    }
}

fn Lanes(comptime T: type, comptime n: usize) type {
    return if (n == 1) T else @Vector(n, T);
}

inline fn load(comptime T: type, comptime n: usize, p: [*]const T) Lanes(T, n) {
    return if (n == 1) p[0] else util.load(@Vector(n, T), p);
}

inline fn store(comptime T: type, comptime n: usize, p: [*]T, v: Lanes(T, n)) void {
    if (n == 1) p[0] = v else util.store(@Vector(n, T), p, v);
}

inline fn splat(comptime T: type, comptime n: usize, x: T) Lanes(T, n) {
    return if (n == 1) x else @splat(x);
}

/// `a*b + c` under the build's contraction policy (`fma.zig`).
inline fn madd(comptime T: type, comptime n: usize, a: Lanes(T, n), b: Lanes(T, n), c: Lanes(T, n)) Lanes(T, n) {
    const fused = (if (n == 1) fma.fuse else fma.fuse_vector) and fma.fusible(T);
    return if (fused) @mulAdd(Lanes(T, n), a, b, c) else a * b + c;
}

/// The kernel set of `T` with tile `geo`; its functions have the vtable's
/// signatures.
pub fn Kernels(comptime T: type, comptime geo: Geometry) type {
    const vl = geo.vl;
    // C for the vtable's types, Zig's own for the rest
    const cc = types.kernelCc(T);
    // `f16` sums in `f32`; vectors are only ever `f32` or `f64`, summed as is
    const A = types.Acc(T);
    if (A != T and vl > 1) @compileError("a widened sum needs scalar lanes");
    const mr = geo.mr();
    const nr = geo.nr;
    return struct {
        pub const unroll_m: i32 = mr;
        pub const unroll_n: i32 = nr;

        /// `gemm_beta`: `C := beta*C` on `m x n`; a zero `beta` writes zeros
        /// without reading `C`.
        pub fn beta(m: BlasLong, n: BlasLong, _: BlasLong, b: T, _: [*]T, _: BlasLong, _: [*]T, _: BlasLong, c: [*]T, ldc: BlasLong) callconv(cc) i32 {
            const height: usize = @intCast(m);
            const ld: usize = @intCast(ldc);
            for (0..@intCast(n)) |j| {
                const col = c[j * ld ..][0..height];
                if (b == 0) {
                    @memset(col, 0);
                } else {
                    for (col) |*v| v.* *= b;
                }
            }
            return 0;
        }

        /// Packs `n` lines of depth `m` in panels of `u`. `across`: line `j`
        /// is element `j` of each of `m` rows `lda` apart (`tcopy`); else it is
        /// column `j`, `lda` apart (`ncopy`).
        fn pack(comptime u: usize, comptime across: bool, m: BlasLong, n: BlasLong, a: [*]T, lda: BlasLong, b: [*]T) void {
            const depth: usize = @intCast(m);
            const lines: usize = @intCast(n);
            const ld: usize = @intCast(lda);
            var c0: usize = 0;
            while (c0 + u <= lines) : (c0 += u) panel(u, across, depth, a, ld, c0, b + depth * c0);
            inline for (comptime tails(u)) |w| {
                if ((lines - c0) & w != 0) {
                    panel(w, across, depth, a, ld, c0, b + depth * c0);
                    c0 += w;
                }
            }
        }

        inline fn panel(comptime w: usize, comptime across: bool, depth: usize, a: [*]const T, ld: usize, c0: usize, dst: [*]T) void {
            if (across) {
                for (0..depth) |l| dst[l * w ..][0..w].* = a[l * ld + c0 ..][0..w].*;
            } else {
                var cols: [w][*]const T = undefined;
                inline for (0..w) |c| cols[c] = a + (c0 + c) * ld;
                for (0..depth) |l| {
                    inline for (0..w) |c| dst[l * w + c] = cols[c][l];
                }
            }
        }

        /// `gemm_incopy`: rows of `op(A) = A^T` (`m` deep, `n` rows), `mr` a panel.
        pub fn incopy(m: BlasLong, n: BlasLong, a: [*]T, lda: BlasLong, b: [*]T) callconv(cc) i32 {
            pack(mr, false, m, n, a, lda, b);
            return 0;
        }
        /// `gemm_itcopy`: rows of `op(A) = A`.
        pub fn itcopy(m: BlasLong, n: BlasLong, a: [*]T, lda: BlasLong, b: [*]T) callconv(cc) i32 {
            pack(mr, true, m, n, a, lda, b);
            return 0;
        }
        /// `gemm_oncopy`: columns of `op(B) = B`, `nr` a panel.
        pub fn oncopy(m: BlasLong, n: BlasLong, a: [*]T, lda: BlasLong, b: [*]T) callconv(cc) i32 {
            pack(nr, false, m, n, a, lda, b);
            return 0;
        }
        /// `gemm_otcopy`: columns of `op(B) = B^T`.
        pub fn otcopy(m: BlasLong, n: BlasLong, a: [*]T, lda: BlasLong, b: [*]T) callconv(cc) i32 {
            pack(nr, true, m, n, a, lda, b);
            return 0;
        }

        /// `gemm_kernel`: `C += alpha * A * B` for `A` packed `bm x bk` (rows
        /// in panels of `mr`) and `B` packed `bk x bn` (columns in panels of
        /// `nr`). Columns outermost, so a B panel stays in L1 while the A
        /// panels stream from L2.
        pub fn kernel(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: T, ba: [*]T, bb: [*]T, c: [*]T, ldc: BlasLong) callconv(cc) i32 {
            const m: usize = @intCast(bm);
            const n: usize = @intCast(bn);
            const k: usize = @intCast(bk);
            const ld: usize = @intCast(ldc);
            var q0: usize = 0;
            while (q0 + nr <= n) : (q0 += nr) rows(nr, m, k, alpha, ba, bb + q0 * k, c + q0 * ld, ld);
            inline for (comptime tails(nr)) |w| {
                if ((n - q0) & w != 0) {
                    rows(w, m, k, alpha, ba, bb + q0 * k, c + q0 * ld, ld);
                    q0 += w;
                }
            }
            return 0;
        }

        /// Every row panel against one B panel of width `wn`: the full tile
        /// inline, the row tails of whole vectors as calls of their own, and
        /// the narrower ones one row at a time (`thin`), so the unrolled tiles
        /// stay few.
        fn rows(comptime wn: usize, m: usize, k: usize, alpha: T, ba: [*]const T, pb: [*]const T, c: [*]T, ld: usize) void {
            var r0: usize = 0;
            while (r0 + mr <= m) : (r0 += mr) @call(.always_inline, tile, .{ mr, wn, k, alpha, ba + r0 * k, pb, c + r0, ld });
            inline for (comptime tails(mr)) |w| {
                if ((m - r0) & w != 0) {
                    if (w >= vl) {
                        @call(.never_inline, tile, .{ w, wn, k, alpha, ba + r0 * k, pb, c + r0, ld });
                    } else {
                        @call(.never_inline, thin, .{ wn, w, k, alpha, ba + r0 * k, pb, c + r0, ld });
                    }
                    r0 += w;
                }
            }
        }

        /// A row panel narrower than a vector, `w` rows: each row against
        /// the `wn` columns as one vector, in two chains over `k`; one
        /// function for every row tail of a column width.
        fn thin(comptime wn: usize, w: usize, k: usize, alpha: T, pa: [*]const T, pb: [*]const T, c: [*]T, ld: usize) void {
            const V = @Vector(wn, T);
            const fused = fma.fuse_vector and fma.fusible(T);
            for (0..w) |r| {
                var acc: [2]V = .{ @splat(0), @splat(0) };
                var l: usize = 0;
                while (l + 2 <= k) : (l += 2) {
                    inline for (0..2) |h| {
                        const x: V = @splat(pa[(l + h) * w + r]);
                        const y = util.load(V, pb + (l + h) * wn);
                        acc[h] = if (fused) @mulAdd(V, x, y, acc[h]) else x * y + acc[h];
                    }
                }
                if (l < k) {
                    const x: V = @splat(pa[l * w + r]);
                    const y = util.load(V, pb + l * wn);
                    acc[0] = if (fused) @mulAdd(V, x, y, acc[0]) else x * y + acc[0];
                }
                const sum: [wn]T = acc[0] + acc[1];
                inline for (0..wn) |j| c[r + j * ld] = madd(T, 1, sum[j], alpha, c[r + j * ld]);
            }
        }

        /// One `wm x wn` tile: `wm / lanes` vectors per column in registers
        /// over all of `k`, then `C += alpha * acc`.
        fn tile(comptime wm: usize, comptime wn: usize, k: usize, alpha: T, pa: [*]const T, pb: [*]const T, c: [*]T, ld: usize) void {
            @setEvalBranchQuota(100_000);
            const lanes = @min(vl, wm);
            const nv = wm / lanes;
            const V = Lanes(A, lanes);
            var acc: [wn][nv]V = undefined;
            inline for (0..wn) |j| inline for (0..nv) |v| {
                acc[j][v] = splat(A, lanes, 0);
            };
            var a = pa;
            var b = pb;
            for (0..k) |_| {
                var av: [nv]V = undefined;
                inline for (0..nv) |v| av[v] = wide(lanes, a + v * lanes);
                inline for (0..wn) |j| {
                    const bv = splat(A, lanes, b[j]);
                    inline for (0..nv) |v| acc[j][v] = madd(A, lanes, av[v], bv, acc[j][v]);
                }
                a += wm;
                b += wn;
            }
            const al = splat(A, lanes, alpha);
            inline for (0..wn) |j| inline for (0..nv) |v| {
                const p = c + j * ld + v * lanes;
                const r = madd(A, lanes, acc[j][v], al, wide(lanes, p));
                if (A == T) store(T, lanes, p, r) else p[0] = @floatCast(r);
            };
        }

        /// `lanes` elements at `p` as `A`.
        inline fn wide(comptime lanes: usize, p: [*]const T) Lanes(A, lanes) {
            return if (A == T) load(T, lanes, p) else p[0];
        }
    };
}

/// The precisions these kernels serve, by vtable prefix.
const precisions = .{ .{ "s", f32 }, .{ "d", f64 }, .{ "q", XDouble } };
/// Those a core's object may replace: `q` slots keep the portable kernels.
const core_precisions = .{ precisions[0], precisions[1] };
const roles = .{ "kernel", "beta", "incopy", "itcopy", "oncopy", "otcopy" };

/// Points `vt`'s real gemm slots at `K(T)` for each of `list`, and its unrolls
/// at the tiles; `unroll(prefix, which)` gives the tile size to store.
fn fill(vt: *Vtable, comptime list: anytype, comptime source: anytype) void {
    inline for (list) |p| {
        inline for (roles) |r| @field(vt, p[0] ++ "gemm_" ++ r) = source.function(p[0], p[1], r);
        @field(vt, p[0] ++ "gemm_unroll_m") = source.unroll(p[0], p[1], "m");
        @field(vt, p[0] ++ "gemm_unroll_n") = source.unroll(p[0], p[1], "n");
        @field(vt, p[0] ++ "gemm_unroll_mn") = @max(source.unroll(p[0], p[1], "m"), source.unroll(p[0], p[1], "n"));
    }
}

/// Puts back the portable kernels' unrolls where `vt`'s slots still hold
/// them: a core's blocking (`params.installCore`), set after them, writes
/// `param.h`'s, and no per-core object may come to replace both.
pub fn settle(vt: *Vtable) void {
    inline for (precisions) |p| {
        const K = Kernels(p[1], geometry(p[1]));
        if (@field(vt, p[0] ++ "gemm_kernel") == &K.kernel) {
            @field(vt, p[0] ++ "gemm_unroll_m") = K.unroll_m;
            @field(vt, p[0] ++ "gemm_unroll_n") = K.unroll_n;
            @field(vt, p[0] ++ "gemm_unroll_mn") = @max(K.unroll_m, K.unroll_n);
        }
    }
}

/// The portable kernels: the main library's own build of these.
pub fn installGeneric(vt: *Vtable) void {
    fill(vt, precisions, struct {
        fn function(comptime _: []const u8, comptime T: type, comptime r: []const u8) @TypeOf(&@field(Kernels(T, geometry(T)), r)) {
            return &@field(Kernels(T, geometry(T)), r);
        }
        fn unroll(comptime _: []const u8, comptime T: type, comptime which: []const u8) i32 {
            return @field(Kernels(T, geometry(T)), "unroll_" ++ which);
        }
    });
}

/// The per-core kernels, for an architecture's registry (`kernels/<arch>.zig`):
/// every per-core object exports them, built for its own features, with
/// their unrolls as data (the tile depends on the object's features, which
/// the main library does not have). They never make a core need an object.
pub fn Registry(comptime Core: type) type {
    return struct {
        pub fn has(comptime core: Core) bool {
            _ = core;
            return false;
        }

        fn Unrolls(comptime T: type) type {
            return struct {
                const m: i32 = Kernels(T, geometry(T)).unroll_m;
                const n: i32 = Kernels(T, geometry(T)).unroll_n;
            };
        }

        pub fn exportCore(comptime core: Core) void {
            inline for (core_precisions) |p| {
                const K = Kernels(p[1], geometry(p[1]));
                inline for (roles) |r| {
                    @export(&@field(K, r), .{ .name = util.kernelSymbol(p[0] ++ "gemm_" ++ r, @tagName(core)), .visibility = .hidden });
                }
                inline for (.{ "m", "n" }) |w| {
                    @export(&@field(Unrolls(p[1]), w), .{ .name = util.kernelSymbol(p[0] ++ "gemm_unroll_" ++ w, @tagName(core)), .visibility = .hidden });
                }
            }
        }

        pub fn install(vt: *Vtable, comptime core: Core) void {
            fill(vt, core_precisions, struct {
                fn function(comptime pre: []const u8, comptime _: type, comptime r: []const u8) @typeInfo(@FieldType(Vtable, pre ++ "gemm_" ++ r)).optional.child {
                    const Fn = @typeInfo(@FieldType(Vtable, pre ++ "gemm_" ++ r)).optional.child;
                    return @extern(Fn, .{ .name = util.kernelSymbol(pre ++ "gemm_" ++ r, @tagName(core)), .visibility = .hidden });
                }
                fn unroll(comptime pre: []const u8, comptime _: type, comptime which: []const u8) i32 {
                    return @extern(*const i32, .{ .name = util.kernelSymbol(pre ++ "gemm_unroll_" ++ which, @tagName(core)), .visibility = .hidden }).*;
                }
            });
        }
    };
}

// ---- tests -----------------------------------------------------------------

fn fillSeq(comptime T: type, buf: []T, seed: usize) void {
    for (buf, 0..) |*v, i| v.* = @floatFromInt(@as(i64, @intCast((i * 7 + seed * 13) % 23)) - 11);
}

/// The packed kernel against plain arithmetic for every panel remainder, the
/// operands packed both ways (exact: small integers). Storage follows the
/// tile, so every size runs whatever the tile.
fn checkLayout(comptime T: type, comptime geo: Geometry) !void {
    const K = Kernels(T, geo);
    const al = std.testing.allocator;
    const sizes = [_]usize{ 1, 2, 3, 5, geo.mr() - 1, geo.mr(), geo.mr() + 1, 2 * geo.mr() + 3, geo.nr - 1, geo.nr + 1, 2 * geo.nr + 1 };
    const ld = std.mem.max(usize, &sizes);
    const kmax = 13;
    const a = try al.alloc(T, ld * kmax);
    defer al.free(a);
    const at = try al.alloc(T, kmax * ld);
    defer al.free(at);
    const b = try al.alloc(T, kmax * ld);
    defer al.free(b);
    const bt = try al.alloc(T, ld * kmax);
    defer al.free(bt);
    const c = try al.alloc(T, ld * ld);
    defer al.free(c);
    const want = try al.alloc(T, ld * ld);
    defer al.free(want);
    const got = try al.alloc(T, ld * ld);
    defer al.free(got);
    const pa = try al.alloc(T, ld * kmax);
    defer al.free(pa);
    const pb = try al.alloc(T, ld * kmax);
    defer al.free(pb);
    fillSeq(T, a, 1);
    fillSeq(T, b, 2);
    fillSeq(T, c, 3);
    // `at` is A transposed (`k x m`, ld `kmax`), `bt` B transposed (`n x k`)
    for (0..ld) |i| for (0..kmax) |l| {
        at[l + i * kmax] = a[i + l * ld];
        bt[i + l * ld] = b[l + i * kmax];
    };
    for (sizes) |m| for (sizes) |n| for ([_]usize{ 1, 4, kmax }) |k| {
        if (m == 0 or n == 0) continue;
        @memcpy(want, c);
        for (0..n) |j| for (0..m) |i| {
            var s: T = 0;
            for (0..k) |l| s += a[i + l * ld] * b[l + j * kmax];
            want[i + j * ld] += 2 * s;
        };
        inline for (.{ false, true }) |trans| {
            @memcpy(got, c);
            const lm: BlasLong = @intCast(ld);
            if (trans) {
                _ = K.incopy(@intCast(k), @intCast(m), at.ptr, kmax, pa.ptr);
                _ = K.otcopy(@intCast(k), @intCast(n), bt.ptr, lm, pb.ptr);
            } else {
                _ = K.itcopy(@intCast(k), @intCast(m), a.ptr, lm, pa.ptr);
                _ = K.oncopy(@intCast(k), @intCast(n), b.ptr, kmax, pb.ptr);
            }
            _ = K.kernel(@intCast(m), @intCast(n), @intCast(k), 2, pa.ptr, pb.ptr, got.ptr, lm);
            for (0..n) |j| for (0..m) |i| {
                if (want[i + j * ld] != got[i + j * ld]) return error.TestUnexpectedResult;
            };
        }
    };
}

test "packed kernels match plain arithmetic for every panel remainder" {
    try checkLayout(f64, geometry(f64));
    try checkLayout(f32, geometry(f32));
    try checkLayout(f64, .{ .vl = 1, .mv = 3, .nr = 6 });
    try checkLayout(f32, .{ .vl = 4, .mv = 2, .nr = 5 });
    // AVX-512 and RVV zvl512b tiles, wider than the old fixed storage
    try checkLayout(f32, .{ .vl = 16, .mv = 2, .nr = 14 });
    try checkLayout(f64, .{ .vl = 8, .mv = 2, .nr = 14 });
}

test "geometry: the most accumulators the registers hold" {
    const cpu = comptime blk: {
        var c = std.Target.Cpu.baseline(.x86_64, sysinfo.os);
        c.features.addFeature(@backingInt(std.Target.x86.Feature.avx2));
        c.features.populateDependencies(&std.Target.x86.all_features);
        break :blk c;
    };
    try std.testing.expectEqual(Geometry{ .vl = 8, .mv = 2, .nr = 6 }, comptime geometryFor(f32, cpu));
    try std.testing.expectEqual(Geometry{ .vl = 4, .mv = 2, .nr = 6 }, comptime geometryFor(f64, cpu));
    const arm64 = comptime std.Target.Cpu.baseline(.aarch64, sysinfo.os);
    try std.testing.expectEqual(Geometry{ .vl = 4, .mv = 4, .nr = 6 }, comptime geometryFor(f32, arm64));
    try std.testing.expectEqual(Geometry{ .vl = 2, .mv = 4, .nr = 6 }, comptime geometryFor(f64, arm64));
    try std.testing.expectEqualSlices(usize, &.{ 4, 2, 1 }, comptime tails(6));
    try std.testing.expectEqualSlices(usize, &.{ 8, 4, 2, 1 }, comptime tails(9));
}

test "beta zero does not read c; other values scale the block only" {
    const K = Kernels(f64, geometry(f64));
    var c = [_]f64{ std.math.nan(f64), 2, 3, 4, 5, 6, 7, 8 };
    _ = K.beta(2, 2, 0, 0, &c, 0, &c, 0, &c, 4);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 3, 4, 0, 0, 7, 8 }, &c);
    c[4] = 5;
    _ = K.beta(2, 2, 0, 2, &c, 0, &c, 0, &c, 4);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 3, 4, 10, 0, 7, 8 }, &c);
}
