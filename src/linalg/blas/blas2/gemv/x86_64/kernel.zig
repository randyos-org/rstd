//! x86_64 gemv kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`). A core's kernel is a `.S` kernel
//! (`asm/`) or the C driver `{s,d,c,z}gemv_{n,t}_4.c` (`../kernel.zig`) around
//! the blocks `uk/index_n.zig` and `uk/index_t.zig` name: the reference's
//! inline-assembly microkernels, or for the AVX-512 cores' `dgemv_n`, whose
//! blocks are intrinsics, the same arithmetic on `@Vector`s. Skylake-X's
//! `sgemv` takes a second, intrinsics block besides Haswell's
//! (`HAVE_SGEMV_N_SKYLAKE_KERNEL`, `_T_`), which `uk_block.Kernels` leaves
//! out, so it is a port here (`ports`). Tests live in `gemv/__root__.zig`;
//! the AVX-512 kernels run only on such hardware (T43).

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../../../types.zig");
const coretypes = @import("../../../cpu/_coretypes.zig");
const fma = @import("../../../fma.zig");
const util = @import("../../../util.zig");
const kernel = @import("../kernel.zig");
const uk = @import("../../../blas1/uk_block.zig");
const index_n = @import("uk/index_n.zig");
const index_t = @import("uk/index_t.zig");
const registry = @import("../../../asm_registry.zig");
const vector = @import("../../../vector.zig");
const Vtable = @import("../../../Vtable.zig");

const BlasLong = types.BlasLong;
const Slot = util.Slot;

pub const Core = coretypes.X86CoreType;

/// The reference compiles the blocks with a 512-bit pass first when
/// `__AVX512CD__` is defined.
const has_avx512cd = std.Target.x86.featureSetHas(sysinfo.target.cpu.features, .avx512cd);

/// `dgemv_n_microk_skylakex-4.c`: `y += alpha * (a0*x0 + a1*x1 + ...)` per row,
/// as vector expressions the compiler contracts (`fuse`): the sum's first
/// product fused with the second, each later one into it, and `y` fused with
/// `sum * alpha`. 512-bit vectors first where `wide`, then 256-bit ones.
pub fn VectorBlockWith(comptime slot: Slot, comptime wide: bool, comptime fuse: bool) type {
    if (slot != .dgemv_n) @compileError(@tagName(slot) ++ " has no vector block");
    const F = f64;
    const V4 = @Vector(4, F);
    const V8 = @Vector(8, F);
    return struct {
        fn rows(comptime V: type, comptime cols: comptime_int, ap: [*]const [*]F, x: [*]F, y: [*]F, alpha: F, i: usize) void {
            var sum = fma.madd(fuse, V, util.load(V, ap[0] + i), @splat(x[0]), util.load(V, ap[1] + i) * @as(V, @splat(x[1])));
            inline for (2..cols) |c| sum = fma.madd(fuse, V, util.load(V, ap[c] + i), @splat(x[c]), sum);
            util.store(V, y + i, fma.madd(fuse, V, sum, @splat(alpha), util.load(V, y + i)));
        }

        fn run(comptime cols: comptime_int, n: BlasLong, ap: [*]const [*]F, x: [*]F, y: [*]F, alpha: F) void {
            const len: usize = @intCast(n);
            var i: usize = 0;
            if (wide and cols == 4) {
                const n5 = len & ~@as(usize, 7);
                while (i < n5) : (i += 8) rows(V8, cols, ap, x, y, alpha, i);
            }
            while (i < len) : (i += 4) rows(V4, cols, ap, x, y, alpha, i);
        }

        pub fn dgemv_kernel_4x4(n: BlasLong, ap: [*]const [*]F, x: [*]F, y: [*]F, alpha: [*]const F) void {
            run(4, n, ap, x, y, alpha[0]);
        }

        pub fn dgemv_kernel_4x2(n: BlasLong, ap: [*]const [*]F, x: [*]F, y: [*]F, alpha: [*]const F) void {
            run(2, n, ap, x, y, alpha[0]);
        }
    };
}

fn VectorBlock(comptime slot: Slot) type {
    return VectorBlockWith(slot, has_avx512cd, fma.fuse);
}

/// Both indexes as one.
const ukidx = struct {
    pub const entries = index_n.entries ++ index_t.entries;
    pub fn find(comptime slot: Slot, comptime core: []const u8) ?uk.Entry {
        return index_n.find(slot, core) orelse index_t.find(slot, core);
    }
};

const K = uk.Kernels(Core, kernel, @import("asm/__root__.zig"), ukidx, VectorBlock);

/// `sgemv_n_microk_skylakex-8.c`'s `sgemv_kernel_n_128` and `_n_64`: per row
/// of `y`, `a*x` fused into a sum column by column, then `y` fused with
/// `sum*alpha`. Each row's arithmetic is its own, so vectors of 16 rows and
/// the masked last one as single rows give the intrinsics' bits.
const SkylakeBlock = struct {
    const V = @Vector(16, f32);

    fn rows(m: BlasLong, n: BlasLong, alpha: f32, a: [*]f32, lda: BlasLong, x: [*]f32, y: [*]f32) void {
        const mu: usize = @intCast(m);
        const nu: usize = @intCast(n);
        const ld: usize = @intCast(lda);
        var i: usize = 0;
        while (i + 16 <= mu) : (i += 16) {
            var acc: V = @splat(0);
            for (0..nu) |j| acc = @mulAdd(V, util.load(V, a + j * ld + i), @splat(x[j]), acc);
            util.store(V, y + i, @mulAdd(V, acc, @splat(alpha), util.load(V, y + i)));
        }
        while (i < mu) : (i += 1) {
            var acc: f32 = 0;
            for (0..nu) |j| acc = @mulAdd(f32, a[j * ld + i], x[j], acc);
            y[i] = @mulAdd(f32, acc, alpha, y[i]);
        }
    }

    pub const n128 = rows;
    pub const n64 = rows;
};

/// `sgemv_t_microk_skylakex_template.c`'s `sgemv_kernel_t_<k>`: `y += alpha*A'*x`
/// for a contiguous `A` of `k` rows, 1 to 8, each output a dot of `k`. Each
/// shape sums its outputs one way in its widest blocks (`f(0, ...)`) and
/// others in its narrower ones and tails, by block (`section`); the
/// permutes only gather the operands. The widest blocks take 16 outputs a
/// vector, the rest one at a time, with the same arithmetic.
const SkylakeTBlock = struct {
    const L = 16;
    const V = @Vector(L, f32);

    pub fn t(m: usize, n: usize, alpha: f32, a: [*]f32, x: [*]f32, y: [*]f32) void {
        switch (m) {
            inline 1...8 => |k| Rows(k).run(n, alpha, a, x, y),
            else => unreachable,
        }
    }

    fn Rows(comptime k: usize) type {
        return struct {
            /// Where output `j` of `n` falls: its shape's block sizes
            /// (`tag_m_16x`, `_8x`, `_4x`, `_2x`) and which sum it takes.
            fn section(j: usize, n: usize) u2 {
                const t16 = n & ~@as(usize, 15);
                const t8 = n & ~@as(usize, 7);
                const t4 = n & ~@as(usize, 3);
                const t2 = n & ~@as(usize, 1);
                return switch (k) {
                    1 => 0,
                    2 => if (j < t16) 0 else 1,
                    3 => if (j < t8) 0 else if (j < t4) 1 else if (j < t2) 2 else 3,
                    4, 8 => if (j < t4) 0 else if (j < t2) 1 else 2,
                    5 => if (j < t4) 0 else 1,
                    6, 7 => if (j < t8) 0 else if (j < t4) 1 else 2,
                    else => unreachable,
                };
            }

            /// The outputs below which every one is in section 0.
            fn bound0(n: usize) usize {
                const bits: usize = switch (k) {
                    1 => 0,
                    2 => 15,
                    3, 6, 7 => 7,
                    4, 5, 8 => 3,
                    else => unreachable,
                };
                return n & ~bits;
            }

            /// Lanes `c`, `k + c`, `2k + c`, ...: element `c` of 16 rows.
            fn gather(comptime c: usize) @Vector(L, i32) {
                var r: [L]i32 = undefined;
                for (&r, 0..) |*v, i| v.* = @intCast(i * k + c);
                return r;
            }

            /// The new `y` of section `sec` from the row `e` (a vector of
            /// rows, or one): `p` its products with `x`, `q` with `x*alpha`.
            inline fn f(comptime sec: u2, comptime T: type, e: [k]T, xs: [k]f32, alpha: f32, y: T) T {
                const S = struct {
                    inline fn v(s: f32) T {
                        return if (T == f32) s else @splat(s);
                    }
                };
                const al = S.v(alpha);
                const z = S.v(0);
                var p: [k]T = undefined;
                var q: [k]T = undefined;
                inline for (0..k) |c| {
                    p[c] = e[c] * S.v(xs[c]);
                    q[c] = e[c] * S.v(xs[c] * alpha);
                }
                switch (k) {
                    1 => return @mulAdd(T, e[0], S.v(alpha * xs[0]), y),
                    2 => return switch (sec) {
                        0 => @mulAdd(T, S.v(xs[1] * alpha), e[1], q[0]) + y,
                        else => (p[1] * al + p[0] * al) + y,
                    },
                    3 => return switch (sec) {
                        0 => @mulAdd(T, S.v(xs[2] * alpha), e[2], @mulAdd(T, S.v(xs[1] * alpha), e[1], q[0])) + y,
                        1 => @mulAdd(T, e[2], S.v(xs[2] * alpha), @mulAdd(T, e[0], S.v(xs[0] * alpha), q[1])) + y,
                        2 => fma.madd(fma.fuse, T, al, (p[0] + p[1]) + (p[2] + z), y),
                        else => fma.madd(fma.fuse, T, al, fma.madd(fma.fuse, T, e[2], S.v(xs[2]), fma.madd(fma.fuse, T, e[0], S.v(xs[0]), p[1])), y),
                    },
                    4 => return switch (sec) {
                        0 => ((q[0] + q[1]) + (q[2] + q[3])) + y,
                        1 => fma.madd(fma.fuse, T, al, (p[0] + p[1]) + (p[2] + p[3]), y),
                        else => fma.madd(fma.fuse, T, al, ((p[0] + p[1]) + p[2]) + p[3], y),
                    },
                    5, 6, 7 => {
                        if (sec == 0) {
                            var acc = @mulAdd(T, e[1], S.v(xs[1]), p[0]);
                            inline for (2..k) |c| acc = @mulAdd(T, e[c], S.v(xs[c]), acc);
                            return @mulAdd(T, acc, al, y);
                        }
                        if (k == 5) return y + ((((q[0] + q[1]) + q[2]) + q[3]) + q[4]);
                        if (sec == 1) {
                            const u = if (k == 6) ((p[0] + z) + (p[1] + p[2])) * al else ((p[0] + p[1]) + (p[4] + p[5])) * al;
                            const w = if (k == 6) ((p[3] + z) + (p[4] + p[5])) * al else ((p[2] + p[3]) + (p[6] + z)) * al;
                            return (u + w) + y;
                        }
                        var r = y;
                        inline for (0..k) |c| r += q[c];
                        return r;
                    },
                    8 => return switch (sec) {
                        0 => (((q[0] + q[1]) + (q[2] + q[3])) + ((q[4] + q[5]) + (q[6] + q[7]))) + y,
                        1 => y + (((q[0] + q[1]) + (q[4] + q[5])) + ((q[2] + q[3]) + (q[6] + q[7]))),
                        else => y + ((((q[0] + q[4]) + (q[1] + q[5])) + (q[2] + q[6])) + (q[3] + q[7])),
                    },
                    else => unreachable,
                }
            }

            fn run(n: usize, alpha: f32, a: [*]f32, x: [*]f32, y: [*]f32) void {
                const xs: [k]f32 = x[0..k].*;
                const b0 = bound0(n);
                var j: usize = 0;
                while (j + L <= b0) : (j += L) {
                    const blk: @Vector(L * k, f32) = (a + j * k)[0 .. L * k].*;
                    var e: [k]V = undefined;
                    inline for (0..k) |c| e[c] = @shuffle(f32, blk, undefined, gather(c));
                    util.store(V, y + j, f(0, V, e, xs, alpha, util.load(V, y + j)));
                }
                while (j < n) : (j += 1) {
                    const e: [k]f32 = (a + j * k)[0..k].*;
                    y[j] = switch (section(j, n)) {
                        inline else => |sec| f(sec, f32, e, xs, alpha, y[j]),
                    };
                }
            }
        };
    }
};

/// The entry of `slot` for Skylake-X in `idx` that is Haswell's block.
fn haswellBlock(comptime idx: type, comptime slot: Slot) uk.Entry {
    return comptime for (idx.entries) |x| {
        if (x.slot == slot and x.block == .assembly and uk.find(&[_]uk.Entry{x}, slot, "skylakex") != null) break x;
    } else unreachable;
}

/// Skylake-X's `sgemv_n` and `sgemv_t`: Haswell's blocks, which its index
/// entries name beside the intrinsics, and `SkylakeBlock`, `SkylakeTBlock`.
pub const Skylake = struct {
    pub fn entry(comptime slot: Slot) util.SlotFn(slot) {
        return switch (slot) {
            .sgemv_n => kernel.skylakeN(haswellBlock(index_n, .sgemv_n), SkylakeBlock),
            .sgemv_t => kernel.skylakeT(haswellBlock(index_t, .sgemv_t), SkylakeTBlock),
            else => @compileError(@tagName(slot) ++ " is not Skylake-X's sgemv"),
        };
    }
};

/// The hand ports: Skylake-X's `sgemv`, whose index names two blocks, and
/// the AVX-512 BF16 cores' `sbgemv_n.c` (`sbgemv_n_microk_cooperlake.c`) as
/// `vector.Half`, which keeps the portable kernels' bits.
pub const ports = [_]registry.Port(Core){
    .{ .module = Skylake, .cores = &.{ .skylakex, .cooperlake, .sapphirerapids }, .slots = &.{ .sgemv_n, .sgemv_t }, .written = .vector },
    .{ .module = vector.ported, .cores = &.{ .cooperlake, .sapphirerapids }, .slots = &.{.sbgemv_n}, .written = .vector },
};

const P = registry.Ports(Core, "x86_64", &ports);

pub fn hasSlot(comptime core: Core, comptime slot: Slot) bool {
    return K.hasSlot(core, slot) or P.hasSlot(core, slot);
}

pub fn has(comptime core: Core) bool {
    return K.has(core) or P.has(core);
}

pub const skips = K.skips;

pub fn exportCore(comptime core: Core) void {
    K.exportCore(core);
    P.exportCore(core);
}

pub fn install(vt: *Vtable, comptime core: Core) void {
    K.install(vt, core);
    P.install(vt, core);
}
