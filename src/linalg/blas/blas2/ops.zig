//! Level 2's element arithmetic and its kernel calls through a vtable: one set
//! for the triangular, banded, rank, gemv and symv drivers.

const blas = @import("../__root__.zig");
const types = @import("../types.zig");
const util = @import("../util.zig");
const Vtable = @import("../Vtable.zig");
const portable = @import("../portable.zig");

const BlasLong = types.BlasLong;
const Trans = types.Trans;
const Op = portable.Op;
const isComplex = types.isComplex;
const Float = types.Float;
const raw = util.raw;

/// Arithmetic on `T` that reads the same for real and complex.
pub fn Num(comptime T: type) type {
    const cplx = isComplex(T);
    return struct {
        pub fn add(a: T, b: T) T {
            return if (cplx) a.add(b) else a + b;
        }
        pub fn sub(a: T, b: T) T {
            return if (cplx) a.sub(b) else a - b;
        }
        pub fn neg(a: T) T {
            return if (cplx) .init(-a.re, -a.im) else -a;
        }
        /// `a`, conjugated for `c`.
        pub fn conj(a: T, comptime c: bool) T {
            return if (cplx and c) a.conjugate() else a;
        }
        /// `a * b`, with `a` conjugated for `c`.
        pub fn mul(a: T, b: T, comptime c: bool) T {
            return if (cplx) conj(a, c).mul(b) else a * b;
        }
        /// `a * b`, each step rounded on the x87 (`util.cmul`): the rank
        /// updates' scalars.
        pub fn mulRounded(a: T, b: T) T {
            if (!cplx) return a * b;
            const p = util.cmul(Float(T), a.re, a.im, b.re, b.im);
            return .init(p[0], p[1]);
        }
        /// `a * f`, `f` real.
        pub fn scale(a: T, f: Float(T)) T {
            return if (cplx) .init(a.re * f, a.im * f) else a * f;
        }
        /// `b / a`, with `a` conjugated for `c`; complex by `util.cdiv`, as
        /// level-3 `trsm` divides.
        pub fn div(b: T, a: T, comptime c: bool) T {
            if (!cplx) return b / a;
            const q = util.cdiv(Float(T), b.re, b.im, a.re, if (c) -a.im else a.im);
            return T.init(q[0], q[1]);
        }
    };
}

/// The kernels for `T` from `vt` (`portable.kernel`). The level-1 calls do
/// nothing for `n <= 0` (a dot gives zero).
pub fn Ops(comptime T: type) type {
    const F = Float(T);
    return struct {
        const Self = @This();

        vt: *const Vtable,
        /// `dtb_entries`: the triangular drivers' block size.
        dtb: BlasLong = 0,

        /// The running vtable, initialized first if need be (`blas.tableOf`).
        pub fn init() Self {
            return .{ .vt = blas.tableOf(T) };
        }

        pub fn copy(s: Self, n: BlasLong, x: [*]T, incx: BlasLong, y: [*]T, incy: BlasLong) void {
            if (n <= 0) return;
            _ = portable.kernel(T, s.vt, .copy_k)(n, raw(T, x), incx, raw(T, y), incy);
        }

        /// `y += alpha*x`, unit stride, `x` conjugated for `cj`.
        pub fn axpy(s: Self, comptime cj: bool, n: BlasLong, alpha: T, x: [*]T, y: [*]T) void {
            if (n <= 0) return;
            const op: Op = comptime if (cj and isComplex(T)) .axpyc_k else .axpy_k;
            _ = @call(.auto, portable.kernel(T, s.vt, op), .{ n, 0, 0 } ++ util.scalarArgs(T, alpha) ++ .{ raw(T, x), 1, raw(T, y), 1, raw(T, y), 0 });
        }

        /// `x . y`, unit stride, `x` conjugated for `cj`.
        pub fn dot(s: Self, comptime cj: bool, n: BlasLong, x: [*]T, y: [*]T) T {
            if (n <= 0) return util.zero(T);
            if (comptime isComplex(T)) {
                const k = portable.kernel(T, s.vt, if (cj) .dotc_k else .dotu_k);
                return types.complexOf(k(n, raw(T, x), 1, raw(T, y), 1));
            }
            return portable.kernel(T, s.vt, .dot_k)(n, x, 1, y, 1);
        }

        /// `y := beta*y` with the scal kernel's flag 0, as the reference's
        /// interfaces scale: a zero `beta` stores zeros, over NaN too.
        pub fn scal(s: Self, n: BlasLong, beta: T, y: [*]T, incy: BlasLong) void {
            if (n <= 0) return;
            const p = raw(T, y);
            _ = @call(.auto, portable.kernel(T, s.vt, .scal_k), .{ n, 0, 0 } ++ util.scalarArgs(T, beta) ++ .{ p, incy, p, 0, p, 0 });
        }

        /// The gemv kernel `op` (`portable.gemvOp`): `y += alpha*op(A)*x`,
        /// `A` of `m` by `n`.
        pub fn gemvK(s: Self, comptime op: Op, m: BlasLong, n: BlasLong, alpha: T, a: [*]T, lda: BlasLong, x: [*]T, incx: BlasLong, y: [*]T, incy: BlasLong, buf: [*]F) void {
            _ = @call(.auto, portable.kernel(T, s.vt, op), .{ m, n, 0 } ++ util.scalarArgs(T, alpha) ++ .{ raw(T, a), lda, raw(T, x), incx, raw(T, y), incy, buf });
        }

        /// `gemvK` of `op`, unit stride, a real `alpha`.
        pub fn gemv(s: Self, comptime op: Trans, m: BlasLong, n: BlasLong, alpha: F, a: [*]T, lda: BlasLong, x: [*]T, y: [*]T, buf: [*]F) void {
            const k = comptime @field(Op, "gemv_" ++ @tagName(op));
            s.gemvK(k, m, n, if (comptime isComplex(T)) .init(alpha, 0) else alpha, a, lda, x, 1, y, 1, buf);
        }

        /// The symv or hemv kernel `op` on the columns `offset` selects of an
        /// `m` by `m` triangle.
        pub fn mv(s: Self, comptime op: Op, m: BlasLong, offset: BlasLong, alpha: T, a: [*]T, lda: BlasLong, x: [*]T, incx: BlasLong, y: [*]T, incy: BlasLong, buf: [*]F) void {
            _ = @call(.auto, portable.kernel(T, s.vt, op), .{ m, offset } ++ util.scalarArgs(T, alpha) ++ .{ raw(T, a), lda, raw(T, x), incx, raw(T, y), incy, buf });
        }
    };
}
