//! Each library's level-1, `gemv` and `gemm` in the ABI of rstd's kernels, for the kernel races' checks.

const rstd = @import("rstd");
const bb = @import("../blas_bench.zig");
const kernels = @import("../kernels.zig");
const lb = @import("../libs.zig");

const blas = rstd.math.kernels.blas;
const Fn = bb.Fn;
const Prec = bb.Prec;
const L1 = kernels.L1;
const l1Slot = kernels.l1Slot;
const Num = lb.Num;
const Src = lb.Src;
const cast = lb.cast;
const col_major = lb.col_major;
const fail = lb.fail;
const lookup = lb.lookup;
const matrix = lb.matrix;
const no_trans = lb.no_trans;
const rstd_alloc = lb.rstd_alloc;
const trans_t = lb.trans_t;
const view = lb.view;

/// Level-1 routine `r` over `p` from `src`, in the kernel ABI of its vtable
/// slot, so `checkL1` and `L1Call` take it as they take a kernel: `f` the
/// library's symbol, `k` the wrapper. `from` makes each library's struct
/// a type of its own; without it they would share one `f`.
pub fn L1Via(comptime src: Src, comptime r: L1, comptime p: Prec) type {
    const F = p.F;
    const T = Num(p);
    const Cx = blas.types.CComplex(F);
    const lib = src != .rstd;
    const s = p.name;
    // the real letter, for `scnrm2`, `csrot` and the like
    const re = if (F == f32) "s" else "d";
    const I = isize;
    const cx = p.complex;
    return switch (r) {
        .axpy => if (!cx) struct {
            const sym = "cblas_" ++ s ++ "axpy";
            const from = src;
            const Sig = fn (c_int, F, [*]F, c_int, [*]F, c_int) callconv(.c) void;
            var f: *const Sig = undefined;
            fn k(n: I, _: I, _: I, a: F, x: [*]F, ix: I, y: [*]F, iy: I, _: [*]F, _: I) callconv(.c) i32 {
                if (lib) f(@intCast(n), a, x, @intCast(ix), y, @intCast(iy)) else blas.axpy(lb.rstd_io, T, a, view(T, x, n, ix).toConst(), view(T, y, n, iy)) catch |e| fail(e);
                return 0;
            }
        } else struct {
            const sym = "cblas_" ++ s ++ "axpy";
            const from = src;
            const Sig = fn (c_int, *const [2]F, [*]F, c_int, [*]F, c_int) callconv(.c) void;
            var f: *const Sig = undefined;
            fn k(n: I, _: I, _: I, ar: F, ai: F, x: [*]F, ix: I, y: [*]F, iy: I, _: [*]F, _: I) callconv(.c) i32 {
                if (lib) f(@intCast(n), &.{ ar, ai }, x, @intCast(ix), y, @intCast(iy)) else blas.axpy(lb.rstd_io, T, .init(ar, ai), view(T, cast(p, x), n, ix).toConst(), view(T, cast(p, y), n, iy)) catch |e| fail(e);
                return 0;
            }
        },
        .dot => if (!cx) struct {
            const sym = "cblas_" ++ s ++ "dot";
            const from = src;
            const Sig = fn (c_int, [*]F, c_int, [*]F, c_int) callconv(.c) F;
            var f: *const Sig = undefined;
            fn k(n: I, x: [*]F, ix: I, y: [*]F, iy: I) callconv(.c) F {
                if (lib) return f(@intCast(n), x, @intCast(ix), y, @intCast(iy));
                return blas.dot(lb.rstd_io, T, view(T, x, n, ix).toConst(), view(T, y, n, iy).toConst()) catch |e| {
                    fail(e);
                    return 0;
                };
            }
        } else struct {
            const sym = "cblas_" ++ s ++ "dotu_sub";
            const from = src;
            const Sig = fn (c_int, [*]F, c_int, [*]F, c_int, *[2]F) callconv(.c) void;
            var f: *const Sig = undefined;
            fn k(n: I, x: [*]F, ix: I, y: [*]F, iy: I) callconv(.c) Cx {
                if (lib) {
                    var o: [2]F = undefined;
                    f(@intCast(n), x, @intCast(ix), y, @intCast(iy), &o);
                    return .{ .re = o[0], .im = o[1] };
                }
                const v = blas.dot(lb.rstd_io, T, view(T, cast(p, x), n, ix).toConst(), view(T, cast(p, y), n, iy).toConst()) catch |e| {
                    fail(e);
                    return .{ .re = 0, .im = 0 };
                };
                return .{ .re = v.re, .im = v.im };
            }
        },
        .scal => if (!cx) struct {
            const sym = "cblas_" ++ s ++ "scal";
            const from = src;
            const Sig = fn (c_int, F, [*]F, c_int) callconv(.c) void;
            var f: *const Sig = undefined;
            fn k(n: I, _: I, _: I, a: F, x: [*]F, ix: I, _: [*]F, _: I, _: [*]F, _: I) callconv(.c) i32 {
                if (lib) f(@intCast(n), a, x, @intCast(ix)) else blas.scal(lb.rstd_io, T, a, view(T, x, n, ix)) catch |e| fail(e);
                return 0;
            }
        } else struct {
            const sym = "cblas_" ++ s ++ "scal";
            const from = src;
            const Sig = fn (c_int, *const [2]F, [*]F, c_int) callconv(.c) void;
            var f: *const Sig = undefined;
            fn k(n: I, _: I, _: I, ar: F, ai: F, x: [*]F, ix: I, _: [*]F, _: I, _: [*]F, _: I) callconv(.c) i32 {
                if (lib) f(@intCast(n), &.{ ar, ai }, x, @intCast(ix)) else blas.scal(lb.rstd_io, T, .init(ar, ai), view(T, cast(p, x), n, ix)) catch |e| fail(e);
                return 0;
            }
        },
        .nrm2, .asum, .sum => struct {
            const sym = "cblas_" ++ (if (cx) re else "") ++ s ++ @tagName(r);
            const from = src;
            const Sig = fn (c_int, [*]F, c_int) callconv(.c) F;
            var f: *const Sig = undefined;
            fn k(n: I, x: [*]F, ix: I) callconv(.c) F {
                if (lib) return f(@intCast(n), x, @intCast(ix));
                const xv = view(T, cast(p, x), n, ix).toConst();
                if (comptime r == .nrm2) return blas.nrm2(lb.rstd_io, T, xv);
                return (if (comptime r == .sum) blas.sum(lb.rstd_io, T, xv) else blas.asum(lb.rstd_io, T, xv)) catch |e| {
                    fail(e);
                    return 0;
                };
            }
        },
        .iamax => struct {
            const sym = "cblas_i" ++ s ++ "amax";
            const from = src;
            // OpenBLAS returns `size_t`, Accelerate `int`: both from 0
            const Sig = fn (c_int, [*]F, c_int) callconv(.c) c_int;
            var f: *const Sig = undefined;
            fn k(n: I, x: [*]F, ix: I) callconv(.c) I {
                if (lib) return @as(I, f(@intCast(n), x, @intCast(ix))) + 1;
                const i = blas.index(lb.rstd_io, T, .largest_magnitude, view(T, cast(p, x), n, ix).toConst()) catch |e| {
                    fail(e);
                    return 0;
                };
                // 1-based, as the libraries' results are made above
                return if (i) |j| @as(I, @intCast(j)) + 1 else 0;
            }
        },
        .copy => struct {
            const sym = "cblas_" ++ s ++ "copy";
            const from = src;
            const Sig = fn (c_int, [*]F, c_int, [*]F, c_int) callconv(.c) void;
            var f: *const Sig = undefined;
            fn k(n: I, x: [*]F, ix: I, y: [*]F, iy: I) callconv(.c) i32 {
                if (lib) f(@intCast(n), x, @intCast(ix), y, @intCast(iy)) else blas.copy(lb.rstd_io, T, view(T, cast(p, x), n, ix).toConst(), view(T, cast(p, y), n, iy)) catch |e| fail(e);
                return 0;
            }
        },
        .swap => if (!cx) struct {
            const sym = "cblas_" ++ s ++ "swap";
            const from = src;
            const Sig = fn (c_int, [*]F, c_int, [*]F, c_int) callconv(.c) void;
            var f: *const Sig = undefined;
            fn k(n: I, _: I, _: I, _: F, x: [*]F, ix: I, y: [*]F, iy: I, _: [*]F, _: I) callconv(.c) i32 {
                if (lib) f(@intCast(n), x, @intCast(ix), y, @intCast(iy)) else blas.swap(lb.rstd_io, T, view(T, x, n, ix), view(T, y, n, iy)) catch |e| fail(e);
                return 0;
            }
        } else struct {
            const sym = "cblas_" ++ s ++ "swap";
            const from = src;
            const Sig = fn (c_int, [*]F, c_int, [*]F, c_int) callconv(.c) void;
            var f: *const Sig = undefined;
            fn k(n: I, _: I, _: I, _: F, _: F, x: [*]F, ix: I, y: [*]F, iy: I, _: [*]F, _: I) callconv(.c) i32 {
                if (lib) f(@intCast(n), x, @intCast(ix), y, @intCast(iy)) else blas.swap(lb.rstd_io, T, view(T, cast(p, x), n, ix), view(T, cast(p, y), n, iy)) catch |e| fail(e);
                return 0;
            }
        },
        .rot => struct {
            const sym = "cblas_" ++ s ++ (if (cx) re else "") ++ "rot";
            const from = src;
            const Sig = fn (c_int, [*]F, c_int, [*]F, c_int, F, F) callconv(.c) void;
            var f: *const Sig = undefined;
            fn k(n: I, x: [*]F, ix: I, y: [*]F, iy: I, c: F, sn: F) callconv(.c) i32 {
                if (lib) f(@intCast(n), x, @intCast(ix), y, @intCast(iy), c, sn) else blas.rot(lb.rstd_io, T, view(T, cast(p, x), n, ix), view(T, cast(p, y), n, iy), c, sn) catch |e| fail(e);
                return 0;
            }
        },
    };
}

/// `src`'s routine `r` as a kernel; null where the library lacks it.
pub fn libL1(comptime src: Src, comptime r: L1, comptime p: Prec) ?Fn(l1Slot(r, p)) {
    const W = L1Via(src, r, p);
    if (src != .rstd) W.f = lookup(src, *const W.Sig, W.sym) orelse return null;
    return &W.k;
}

/// `gemv` over `p` from `src` in the kernel ABI of `?gemv_n`/`?gemv_t`:
/// `y += alpha*op(A)*x` (`beta` 1).
pub fn GemvVia(comptime src: Src, comptime p: Prec, comptime t: []const u8) type {
    const F = p.F;
    const T = Num(p);
    const lib = src != .rstd;
    const I = isize;
    const tr = if (t[0] == 't') trans_t else no_trans;
    const rt: blas.Trans = if (t[0] == 't') .t else .n;
    return if (!p.complex) struct {
        const sym = "cblas_" ++ p.name ++ "gemv";
        const from = src;
        const Sig = fn (c_int, c_int, c_int, c_int, F, [*]F, c_int, [*]F, c_int, F, [*]F, c_int) callconv(.c) void;
        var f: *const Sig = undefined;
        pub fn k(m: I, n: I, _: I, a: F, ma: [*]F, lda: I, x: [*]F, ix: I, y: [*]F, iy: I, _: [*]F) callconv(.c) i32 {
            if (lib) f(col_major, tr, @intCast(m), @intCast(n), a, ma, @intCast(lda), x, @intCast(ix), 1, y, @intCast(iy)) else blas.gemv(rstd_alloc, lb.rstd_io, T, rt, a, matrix(T, ma, m, n, lda), @intCast(lda), view(T, x, if (rt.transposed()) m else n, ix).toConst(), 1, view(T, y, if (rt.transposed()) n else m, iy)) catch |e| fail(e);
            return 0;
        }
    } else struct {
        const sym = "cblas_" ++ p.name ++ "gemv";
        const from = src;
        const Sig = fn (c_int, c_int, c_int, c_int, *const [2]F, [*]F, c_int, [*]F, c_int, *const [2]F, [*]F, c_int) callconv(.c) void;
        var f: *const Sig = undefined;
        pub fn k(m: I, n: I, _: I, ar: F, ai: F, ma: [*]F, lda: I, x: [*]F, ix: I, y: [*]F, iy: I, _: [*]F) callconv(.c) i32 {
            if (lib) f(col_major, tr, @intCast(m), @intCast(n), &.{ ar, ai }, ma, @intCast(lda), x, @intCast(ix), &.{ 1, 0 }, y, @intCast(iy)) else blas.gemv(rstd_alloc, lb.rstd_io, T, rt, .init(ar, ai), matrix(T, cast(p, ma), m, n, lda), @intCast(lda), view(T, cast(p, x), if (rt.transposed()) m else n, ix).toConst(), .init(1, 0), view(T, cast(p, y), if (rt.transposed()) n else m, iy)) catch |e| fail(e);
            return 0;
        }
    };
}

pub fn libGemv(comptime src: Src, comptime p: Prec, comptime t: []const u8) ?Fn(p.name ++ "gemv_" ++ t) {
    const W = GemvVia(src, p, t);
    if (src != .rstd) W.f = lookup(src, *const W.Sig, W.sym) orelse return null;
    return &W.k;
}

/// Square `NN` gemm, `C += A*B`, column-major with leading dimension `n`.
pub fn GemmFn(comptime p: Prec) type {
    return *const fn (n: usize, a: [*]p.F, b: [*]p.F, c: [*]p.F) void;
}

pub fn GemmVia(comptime src: Src, comptime p: Prec) type {
    const F = p.F;
    const lib = src != .rstd;
    return if (!p.complex) struct {
        const sym = "cblas_" ++ p.name ++ "gemm";
        const from = src;
        const Sig = fn (c_int, c_int, c_int, c_int, c_int, c_int, F, [*]F, c_int, [*]F, c_int, F, [*]F, c_int) callconv(.c) void;
        var f: *const Sig = undefined;
        pub fn k(n: usize, a: [*]F, b: [*]F, c: [*]F) void {
            const ni: c_int = @intCast(n);
            if (lib) f(col_major, no_trans, no_trans, ni, ni, ni, 1, a, ni, b, ni, 1, c, ni) else blas.gemm(rstd_alloc, lb.rstd_io, F, .n, .n, n, n, n, 1, a[0 .. n * n], n, b[0 .. n * n], n, 1, c[0 .. n * n], n) catch |e| fail(e);
        }
    } else struct {
        const sym = "cblas_" ++ p.name ++ "gemm";
        const from = src;
        const Sig = fn (c_int, c_int, c_int, c_int, c_int, c_int, *const [2]F, [*]F, c_int, [*]F, c_int, *const [2]F, [*]F, c_int) callconv(.c) void;
        var f: *const Sig = undefined;
        pub fn k(n: usize, a: [*]F, b: [*]F, c: [*]F) void {
            const ni: c_int = @intCast(n);
            const one: Num(p) = .init(1, 0);
            if (lib) f(col_major, no_trans, no_trans, ni, ni, ni, &.{ 1, 0 }, a, ni, b, ni, &.{ 1, 0 }, c, ni) else blas.gemm(rstd_alloc, lb.rstd_io, Num(p), .n, .n, n, n, n, one, cast(p, a)[0 .. n * n], n, cast(p, b)[0 .. n * n], n, one, cast(p, c)[0 .. n * n], n) catch |e| fail(e);
        }
    };
}

/// `src`'s gemm; null where the library has none.
pub fn libGemm(comptime src: Src, comptime p: Prec) ?GemmFn(p) {
    const W = GemmVia(src, p);
    if (src != .rstd) W.f = lookup(src, *const W.Sig, W.sym) orelse return null;
    return &W.k;
}
