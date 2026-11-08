//! `?TRMV`, `?TRSV`, `?TBMV`, `?TBSV`, `?TPMV`, `?TPSV`: the checks of
//! `interface/{trmv,trsv,tbmv,tbsv,tpmv,tpsv}.c` over
//! `driver/level2/{,z}{trmv,trsv,tbmv,tbsv,tpmv,tpsv}_{U,L}.c`, loops over the
//! copy, axpy, dot and gemv kernels. `_U` also does the transposed lower
//! matrix, as upstream; the products thread (`Split`), the solves are serial.

const std = @import("std");
const sysinfo = @import("builtin");
const blas = @import("../__root__.zig");
const types = @import("../types.zig");
const util = @import("../util.zig");
const Vtable = @import("../Vtable.zig");
const portable = @import("../portable.zig");
const accel = @import("../apple_accel.zig");
const threading = @import("../threading.zig");
const interface = @import("interface.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const BlasLong = types.BlasLong;
const Trans = types.Trans;
const Uplo = types.Uplo;
const Diag = types.Diag;

pub const Error = interface.Error;

/// Work space kept on the stack; bigger problems use the allocator.
pub const stack_bytes = interface.stack_bytes;

/// `dtb_entries` of a hand-built vtable that has none.
const default_dtb = 64;

// ---- element helpers -------------------------------------------------------

const isComplex = types.isComplex;
const Float = types.Float;

const raw = util.raw;

const at = util.at;

const off = util.off;

const first = util.first;

/// Arithmetic on `T` that reads the same for real and complex.
fn Num(comptime T: type) type {
    const cplx = isComplex(T);
    return struct {
        fn add(a: T, b: T) T {
            return if (cplx) a.add(b) else a + b;
        }
        fn sub(a: T, b: T) T {
            return if (cplx) a.sub(b) else a - b;
        }
        fn neg(a: T) T {
            return if (cplx) .init(-a.re, -a.im) else -a;
        }
        fn conj(a: T, comptime c: bool) T {
            return if (cplx and c) a.conjugate() else a;
        }
        /// `a * b`, with `a` conjugated for `c`.
        fn mul(a: T, b: T, comptime c: bool) T {
            return if (cplx) conj(a, c).mul(b) else a * b;
        }
        /// `b / a`, with `a` conjugated for `c`; complex by `util.cdiv`, as
        /// level-3 `trsm` divides.
        fn div(b: T, a: T, comptime c: bool) T {
            if (!cplx) return b / a;
            const q = util.cdiv(Float(T), b.re, b.im, a.re, if (c) -a.im else a.im);
            return T.init(q[0], q[1]);
        }
    };
}

// ---- kernels ---------------------------------------------------------------

/// The kernels for `T` (`portable.kernel`), on unit-stride vectors.
fn Ops(comptime T: type) type {
    const F = Float(T);
    return struct {
        const Self = @This();

        vt: *const Vtable,
        dtb: BlasLong,

        fn copy(s: Self, n: BlasLong, x: [*]T, incx: BlasLong, y: [*]T, incy: BlasLong) void {
            _ = portable.kernel(T, s.vt, .copy_k)(n, raw(T, x), incx, raw(T, y), incy);
        }

        /// `y += alpha*x`, `x` conjugated for `cj`.
        fn axpy(s: Self, comptime cj: bool, n: BlasLong, alpha: T, x: [*]T, y: [*]T) void {
            const op: portable.Op = comptime if (cj and isComplex(T)) .axpyc_k else .axpy_k;
            _ = @call(.auto, portable.kernel(T, s.vt, op), .{ n, 0, 0 } ++ util.scalarArgs(T, alpha) ++ .{ raw(T, x), 1, raw(T, y), 1, raw(T, y), 0 });
        }

        /// `x . y`, `x` conjugated for `cj`.
        fn dot(s: Self, comptime cj: bool, n: BlasLong, x: [*]T, y: [*]T) T {
            if (comptime isComplex(T)) {
                const k = portable.kernel(T, s.vt, if (cj) .dotc_k else .dotu_k);
                return types.complexOf(k(n, raw(T, x), 1, raw(T, y), 1));
            }
            return portable.kernel(T, s.vt, .dot_k)(n, x, 1, y, 1);
        }

        /// `y += alpha*op(A)*x`, `A` of `m` by `n`.
        fn gemv(s: Self, comptime op: Trans, m: BlasLong, n: BlasLong, alpha: F, a: [*]T, lda: BlasLong, x: [*]T, y: [*]T, buf: [*]F) void {
            const k = portable.kernel(T, s.vt, @field(portable.Op, "gemv_" ++ @tagName(op)));
            if (comptime isComplex(T)) {
                _ = k(m, n, 0, alpha, 0, raw(T, a), lda, raw(T, x), 1, raw(T, y), 1, buf);
            } else {
                _ = k(m, n, 0, alpha, a, lda, x, 1, y, 1, buf);
            }
        }
    };
}

/// Real types have no `r` and `c`: they are `n` and `t`.
fn canon(comptime T: type, op: Trans) Trans {
    if (isComplex(T)) return op;
    return if (op.transposed()) .t else .n;
}

// ---- drivers ---------------------------------------------------------------

/// Which reference driver file: `_U` or `_L`.
const File = enum { u, l };

/// Shared by every driver: `m` is `n` for the full and packed ones.
fn Drive(comptime T: type, comptime file: File, comptime op: Trans, comptime unit: bool) type {
    const F = Float(T);
    const N = Num(T);
    const tr = comptime op.transposed();
    const cj = comptime op.conjA();
    return struct {
        /// `trmv_{U,L}.c`: `b := op(A) * b`.
        fn trmv(o: Ops(T), m: BlasLong, a: [*]T, lda: BlasLong, b: [*]T, g: [*]F) void {
            const dtb = o.dtb;
            if (file == .u) {
                var is: BlasLong = 0;
                while (is < m) : (is += dtb) {
                    const min_i = @min(m - is, dtb);
                    if (!tr and is > 0) o.gemv(op, is, min_i, 1, off(T, a, is * lda), lda, off(T, b, is), b, g);
                    var i: BlasLong = 0;
                    while (i < min_i) : (i += 1) {
                        const aa = off(T, a, is + (i + is) * lda);
                        const bb = off(T, b, is);
                        const bi = at(T, bb, i);
                        if (!tr and i > 0) o.axpy(cj, i, bi.*, aa, bb);
                        if (!unit) bi.* = N.mul(at(T, aa, i).*, bi.*, cj);
                        if (tr and i < min_i - 1) bi.* = N.add(bi.*, o.dot(cj, min_i - i - 1, off(T, aa, i + 1), off(T, bb, i + 1)));
                    }
                    if (tr and m - is > min_i) {
                        o.gemv(op, m - is - min_i, min_i, 1, off(T, a, is + min_i + is * lda), lda, off(T, b, is + min_i), off(T, b, is), g);
                    }
                }
            } else {
                var is = m;
                while (is > 0) : (is -= dtb) {
                    const min_i = @min(is, dtb);
                    if (!tr and m - is > 0) {
                        o.gemv(op, m - is, min_i, 1, off(T, a, is + (is - min_i) * lda), lda, off(T, b, is - min_i), off(T, b, is), g);
                    }
                    var i: BlasLong = 0;
                    while (i < min_i) : (i += 1) {
                        const aa = off(T, a, (is - i - 1) + (is - i - 1) * lda);
                        const bb = off(T, b, is - i - 1);
                        if (!tr and i > 0) o.axpy(cj, i, bb[0], off(T, aa, 1), off(T, bb, 1));
                        if (!unit) bb[0] = N.mul(aa[0], bb[0], cj);
                        if (tr and i < min_i - 1) {
                            const len = min_i - i - 1;
                            bb[0] = N.add(bb[0], o.dot(cj, len, off(T, aa, -len), off(T, bb, -len)));
                        }
                    }
                    if (tr and is - min_i > 0) o.gemv(op, is - min_i, min_i, 1, off(T, a, (is - min_i) * lda), lda, b, off(T, b, is - min_i), g);
                }
            }
        }

        /// `trsv_{U,L}.c`: `b := inv(op(A)) * b`.
        fn trsv(o: Ops(T), m: BlasLong, a: [*]T, lda: BlasLong, b: [*]T, g: [*]F) void {
            const dtb = o.dtb;
            if (file == .u) {
                var is = m;
                while (is > 0) : (is -= dtb) {
                    const min_i = @min(is, dtb);
                    if (tr and m - is > 0) {
                        o.gemv(op, m - is, min_i, -1, off(T, a, is + (is - min_i) * lda), lda, off(T, b, is), off(T, b, is - min_i), g);
                    }
                    var i: BlasLong = 0;
                    while (i < min_i) : (i += 1) {
                        const aa = off(T, a, (is - i - 1) + (is - i - 1) * lda);
                        const bb = off(T, b, is - i - 1);
                        if (tr and i > 0) bb[0] = N.sub(bb[0], o.dot(cj, i, off(T, aa, 1), off(T, bb, 1)));
                        if (!unit) bb[0] = N.div(bb[0], aa[0], cj);
                        if (!tr and i < min_i - 1) {
                            const len = min_i - i - 1;
                            o.axpy(cj, len, N.neg(bb[0]), off(T, aa, -len), off(T, bb, -len));
                        }
                    }
                    if (!tr and is - min_i > 0) o.gemv(op, is - min_i, min_i, -1, off(T, a, (is - min_i) * lda), lda, off(T, b, is - min_i), b, g);
                }
            } else {
                var is: BlasLong = 0;
                while (is < m) : (is += dtb) {
                    const min_i = @min(m - is, dtb);
                    if (tr and is > 0) o.gemv(op, is, min_i, -1, off(T, a, is * lda), lda, b, off(T, b, is), g);
                    var i: BlasLong = 0;
                    while (i < min_i) : (i += 1) {
                        const aa = off(T, a, is + (i + is) * lda);
                        const bb = off(T, b, is);
                        const bi = at(T, bb, i);
                        if (tr and i > 0) bi.* = N.sub(bi.*, o.dot(cj, i, aa, bb));
                        if (!unit) bi.* = N.div(bi.*, at(T, aa, i).*, cj);
                        if (!tr and i < min_i - 1) o.axpy(cj, min_i - i - 1, N.neg(bi.*), off(T, aa, i + 1), off(T, bb, i + 1));
                    }
                    if (!tr and m - is > min_i) {
                        o.gemv(op, m - is - min_i, min_i, -1, off(T, a, is + min_i + is * lda), lda, off(T, b, is), off(T, b, is + min_i), g);
                    }
                }
            }
        }

        /// `tbmv_{U,L}.c`: `b := op(A) * b`, `A` of `n` with `k` off-diagonals.
        fn tbmv(o: Ops(T), n: BlasLong, k: BlasLong, a0: [*]T, lda: BlasLong, b: [*]T) void {
            var a = a0;
            if (file == .u) {
                var i: BlasLong = 0;
                while (i < n) : (i += 1) {
                    const bi = at(T, b, i);
                    if (!tr) {
                        const len = @min(i, k);
                        if (len > 0) o.axpy(cj, len, bi.*, off(T, a, k - len), off(T, b, i - len));
                        if (!unit) bi.* = N.mul(at(T, a, k).*, bi.*, cj);
                    } else {
                        if (!unit) bi.* = N.mul(a[0], bi.*, cj);
                        const len = @min(n - i - 1, k);
                        if (len > 0) bi.* = N.add(bi.*, o.dot(cj, len, off(T, a, 1), off(T, b, i + 1)));
                    }
                    a = off(T, a, lda);
                }
            } else {
                a = off(T, a, (n - 1) * lda);
                var i = n - 1;
                while (i >= 0) : (i -= 1) {
                    const bi = at(T, b, i);
                    if (!tr) {
                        const len = @min(n - i - 1, k);
                        if (len > 0) o.axpy(cj, len, bi.*, off(T, a, 1), off(T, b, i + 1));
                        if (!unit) bi.* = N.mul(a[0], bi.*, cj);
                    } else {
                        if (!unit) bi.* = N.mul(at(T, a, k).*, bi.*, cj);
                        const len = @min(i, k);
                        if (len > 0) bi.* = N.add(bi.*, o.dot(cj, len, off(T, a, k - len), off(T, b, i - len)));
                    }
                    a = off(T, a, -lda);
                }
            }
        }

        /// `tbsv_{U,L}.c`: `b := inv(op(A)) * b`.
        fn tbsv(o: Ops(T), n: BlasLong, k: BlasLong, a0: [*]T, lda: BlasLong, b: [*]T) void {
            var a = a0;
            if (file == .u) {
                a = off(T, a, (n - 1) * lda);
                var i = n - 1;
                while (i >= 0) : (i -= 1) {
                    const bi = at(T, b, i);
                    if (tr) {
                        const len = @min(n - i - 1, k);
                        if (len > 0) bi.* = N.sub(bi.*, o.dot(cj, len, off(T, a, 1), off(T, b, i + 1)));
                        if (!unit) bi.* = N.div(bi.*, a[0], cj);
                    } else {
                        if (!unit) bi.* = N.div(bi.*, at(T, a, k).*, cj);
                        const len = @min(i, k);
                        if (len > 0) o.axpy(cj, len, N.neg(bi.*), off(T, a, k - len), off(T, b, i - len));
                    }
                    a = off(T, a, -lda);
                }
            } else {
                var i: BlasLong = 0;
                while (i < n) : (i += 1) {
                    const bi = at(T, b, i);
                    if (tr) {
                        const len = @min(i, k);
                        if (len > 0) bi.* = N.sub(bi.*, o.dot(cj, len, off(T, a, k - len), off(T, b, i - len)));
                        if (!unit) bi.* = N.div(bi.*, at(T, a, k).*, cj);
                    } else {
                        if (!unit) bi.* = N.div(bi.*, a[0], cj);
                        const len = @min(n - i - 1, k);
                        if (len > 0) o.axpy(cj, len, N.neg(bi.*), off(T, a, 1), off(T, b, i + 1));
                    }
                    a = off(T, a, lda);
                }
            }
        }

        /// `tpmv_{U,L}.c`: `b := op(A) * b`, `A` packed by columns.
        fn tpmv(o: Ops(T), m: BlasLong, a0: [*]T, b: [*]T) void {
            var a = a0;
            if (file == .u) {
                var i: BlasLong = 0;
                while (i < m) : (i += 1) {
                    const bi = at(T, b, i);
                    if (!tr and i > 0) o.axpy(cj, i, bi.*, a, b);
                    if (!unit) bi.* = N.mul(at(T, a, if (tr) 0 else i).*, bi.*, cj);
                    if (tr and i < m - 1) bi.* = N.add(bi.*, o.dot(cj, m - i - 1, off(T, a, 1), off(T, b, i + 1)));
                    a = off(T, a, if (tr) m - i else i + 1);
                }
            } else {
                a = off(T, a, @divExact(m * (m + 1), 2) - 1);
                var i: BlasLong = 0;
                while (i < m) : (i += 1) {
                    const bi = at(T, b, m - i - 1);
                    if (!tr and i > 0) o.axpy(cj, i, bi.*, off(T, a, 1), off(T, b, m - i));
                    if (!unit) bi.* = N.mul(a[0], bi.*, cj);
                    if (tr and i < m - 1) bi.* = N.add(bi.*, o.dot(cj, m - i - 1, off(T, a, -(m - i - 1)), b));
                    a = off(T, a, -(if (tr) m - i else i + 2));
                }
            }
        }

        /// `tpsv_{U,L}.c`: `b := inv(op(A)) * b`.
        fn tpsv(o: Ops(T), m: BlasLong, a0: [*]T, b: [*]T) void {
            var a = a0;
            if (file == .u) {
                a = off(T, a, @divExact(m * (m + 1), 2) - 1);
                var i: BlasLong = 0;
                while (i < m) : (i += 1) {
                    const bi = at(T, b, m - i - 1);
                    if (tr and i > 0) bi.* = N.sub(bi.*, o.dot(cj, i, off(T, a, 1), off(T, b, m - i)));
                    if (!unit) bi.* = N.div(bi.*, a[0], cj);
                    if (!tr and i < m - 1) o.axpy(cj, m - i - 1, N.neg(bi.*), off(T, a, -(m - i - 1)), b);
                    a = off(T, a, -(if (tr) i + 2 else m - i));
                }
            } else {
                var i: BlasLong = 0;
                while (i < m) : (i += 1) {
                    const bi = at(T, b, i);
                    if (tr and i > 0) bi.* = N.sub(bi.*, o.dot(cj, i, a, b));
                    if (!unit) bi.* = N.div(bi.*, at(T, a, if (tr) i else 0).*, cj);
                    if (!tr and i < m - 1) o.axpy(cj, m - i - 1, N.neg(bi.*), off(T, a, 1), off(T, b, i + 1));
                    a = off(T, a, if (tr) i + 1 else m - i);
                }
            }
        }
    };
}

// ---- interface -------------------------------------------------------------

const Kind = enum { trmv, trsv, tbmv, tbsv, tpmv, tpsv };

/// Elements of `F` a copy of `x` takes, rounded to a 64-byte multiple so
/// the `gemv` buffer after it stays aligned (`f80` on x86 is 12 bytes).
fn copyElems(comptime T: type, n: usize) usize {
    const per = 64 / std.math.gcd(64, @sizeOf(Float(T)));
    const comp: usize = if (isComplex(T)) 2 else 1;
    return std.mem.alignForward(usize, n * comp, per);
}

fn run(
    comptime T: type,
    comptime kind: Kind,
    comptime file: File,
    comptime op: Trans,
    comptime unit: bool,
    o: Ops(T),
    n: BlasLong,
    k: BlasLong,
    a: [*]T,
    lda: BlasLong,
    x: [*]T,
    incx: BlasLong,
    buf: [*]Float(T),
) void {
    const D = Drive(T, file, op, unit);
    var b = x;
    var g = buf;
    if (incx != 1) {
        b = @ptrCast(@alignCast(buf));
        g = buf + copyElems(T, @intCast(n));
        o.copy(n, x, incx, b, 1);
    }
    switch (kind) {
        .trmv => D.trmv(o, n, a, lda, b, g),
        .trsv => D.trsv(o, n, a, lda, b, g),
        .tbmv => D.tbmv(o, n, k, a, lda, b),
        .tbsv => D.tbsv(o, n, k, a, lda, b),
        .tpmv => D.tpmv(o, n, a, b),
        .tpsv => D.tpsv(o, n, a, b),
    }
    if (incx != 1) o.copy(n, b, 1, x, incx);
}

/// The public routines: Accelerate where it runs them (`apple_accel.zig`),
/// after `call`'s checks, else `call` on the selected kernels.
fn public(allocator: Allocator, io: Io, comptime T: type, comptime kind: Kind, uplo: Uplo, trans: Trans, diag: Diag, n: isize, k: isize, a: [*]T, lda: isize, x: [*]T, incx: isize) Error!void {
    if (comptime accel.covers(T)) {
        if (!try checks(T, kind, uplo, trans, n, k, lda, incx)) return;
        if (accel.tri(T, @field(accel.TriKind, @tagName(kind)), uplo, trans, diag, n, k, a, lda, x, incx)) return;
    }
    return call(allocator, io, T, kind, blas.vtable(), uplo, trans, diag, n, k, a, lda, x, incx);
}

/// The reference's checks, and nothing past the address space; false for
/// nothing to do.
fn checks(comptime T: type, comptime kind: Kind, uplo: Uplo, trans: Trans, n: isize, k: isize, lda: isize, incx: isize) interface.ArgError!bool {
    const banded = kind == .tbmv or kind == .tbsv;
    const packed_ = kind == .tpmv or kind == .tpsv;
    if (uplo.isRev()) return error.InvalidUplo;
    if (trans.conjX()) return error.InvalidTrans;
    if (n < 0) return error.InvalidN;
    if (banded) {
        if (k < 0) return error.InvalidK;
        if (lda < k + 1) return error.InvalidLda;
    } else if (!packed_ and lda < @max(1, n)) return error.InvalidLda;
    if (incx == 0) return error.InvalidIncx;
    if (n == 0) return false;
    // nothing past the address space (`util.fits`)
    if (packed_) {
        if (!util.fitsMatrix(n, n, n, @sizeOf(T))) return error.InvalidN;
    } else if (!util.fitsMatrix(if (banded) lda else n, n, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fits(n, incx, @sizeOf(T))) return error.InvalidIncx;
    return true;
}

/// The checks, work space and driver choice common to the six routines;
/// `vt` is the kernels' vtable.
fn call(
    allocator: Allocator,
    io: Io,
    comptime T: type,
    comptime kind: Kind,
    vt: *const Vtable,
    uplo: Uplo,
    trans: Trans,
    diag: Diag,
    n: isize,
    k: isize,
    a: [*]T,
    lda: isize,
    x: [*]T,
    incx: isize,
) Error!void {
    const F = Float(T);
    const banded = kind == .tbmv or kind == .tbsv;
    const packed_ = kind == .tpmv or kind == .tpsv;
    if (!try checks(T, kind, uplo, trans, n, k, lda, incx)) return;

    const kn: BlasLong = n;
    const o: Ops(T) = .{ .vt = vt, .dtb = if (vt.dtb_entries > 0) vt.dtb_entries else default_dtb };
    const px = first(T, x, kn, incx);
    if (comptime kind == .trmv or kind == .tbmv or kind == .tpmv) {
        const nu: usize = @intCast(n);
        // `{trmv,tbmv,tpmv}_thread.c` split the columns on multiples of 8,
        // balanced by length but for the band
        const shape: threading.Shape = if (kind == .tbmv) .even else if (uplo == .l) .falling else .rising;
        const p = threading.plan(if (kind == .trmv) .trmv else .tpmv, T, nu, nu, shape, 8);
        if (p.chunks.count > 1) return split(allocator, io, T, kind, o, &p, uplo, trans, diag, kn, k, a, lda, px, incx);
    }

    // x's copy, then `gemv`'s buffer (the full routines, past one block); the
    // 40-element pad is ours (`interface/ztrmv.c:245` pads its small case by 40)
    const un: usize = @intCast(n);
    var elems: usize = if (incx != 1) copyElems(T, un) else 0;
    if (!banded and !packed_ and kn > o.dtb) {
        const comp: usize = if (isComplex(T)) 2 else 1;
        elems += comp * (2 * un + @as(usize, @intCast(o.dtb))) + 40;
    }
    // the copy of `x` is rounded to 64 bytes, so the gemv buffer after it
    // keeps the stack's alignment
    var stack: [stack_bytes]u8 align(64) = undefined;
    const sc = try util.Scratch(F, 64, 64).init(allocator, &stack, elems);
    defer sc.deinit(allocator);
    const buf = sc.slice();
    // `zgemv_t.S` and its siblings return NaN if the buffer holds NaN or
    // infinity; cleared as `interface.gemv` does
    if (comptime isComplex(T) and sysinfo.target.cpu.arch == .x86_64) @memset(buf, 0);

    // `_U` for upper untransposed or lower transposed
    const fl: File = if ((uplo == .u) != trans.transposed()) .u else .l;
    switch (fl) {
        inline else => |f| switch (trans) {
            inline .n, .t, .r, .c => |t| switch (diag) {
                inline else => |d| run(T, kind, f, canon(T, t), d == .u, o, kn, k, a, lda, px, incx, buf.ptr),
            },
            .o, .u, .s, .d => unreachable,
        },
    }
}

/// The products over more than one chunk of columns: each chunk's part of
/// `op(A) * x` into a partial vector of its own, the partials then folded into
/// `x` in chunk order (`threading.Partials`). With a transpose a chunk's part
/// is its own rows, so nothing overlaps.
fn split(allocator: Allocator, io: Io, comptime T: type, comptime kind: Kind, o: Ops(T), p: *const threading.Plan, uplo: Uplo, trans: Trans, diag: Diag, n: BlasLong, k: BlasLong, a: [*]T, lda: BlasLong, x: [*]T, incx: BlasLong) (Allocator.Error || Io.Cancelable)!void {
    const F = Float(T);
    const un: usize = @intCast(n);
    const comp: usize = if (isComplex(T)) 2 else 1;
    // a strided `x`'s copy, then each chunk's `gemv` buffer, as `run` lays
    // them out
    const xe: usize = if (incx != 1) copyElems(T, un) else 0;
    const ge: usize = if (kind == .trmv) copyElems(F, comp * (2 * un + @as(usize, @intCast(o.dtb))) + 40) else 0;
    const buf = try allocator.alignedAlloc(F, .@"64", xe + p.chunks.count * ge);
    defer allocator.free(buf);
    // as `call` clears it
    if (comptime isComplex(T) and sysinfo.target.cpu.arch == .x86_64) @memset(buf, 0);
    const xs: [*]T = if (incx != 1) @ptrCast(buf.ptr) else x;
    if (incx != 1) o.copy(n, x, incx, xs, 1);
    switch (uplo == .l) {
        inline else => |lower| switch (trans) {
            inline .n, .t, .r, .c => |t| switch (diag) {
                inline else => |d| {
                    const S = Split(T, kind, lower, canon(T, t), d == .u);
                    var s: S = .{ .o = o, .n = n, .k = k, .a = a, .lda = lda, .x = xs, .g = buf.ptr + xe, .ge = ge, .parts = undefined };
                    const parts: threading.Partials(T) = try .init(allocator, p, @as(*const S, &s), S.rows);
                    defer parts.deinit(allocator);
                    s.parts = &parts;
                    try threading.forEach(io, p, @as(*const S, &s), S.one);
                    // `x` is the output: zeroed, then the partials added
                    for (0..un) |i| at(T, x, @as(BlasLong, @intCast(i)) * incx).* = std.mem.zeroes(T);
                    try parts.fold(io, x, incx);
                },
            },
            .o, .u, .s, .d => unreachable,
        },
    }
}

/// A chunk of columns of `op(A) * x`, out of place: the `trmv_kernel`s of
/// `driver/level2/trmv_thread.c`, `tbmv_thread.c` and `tpmv_thread.c`.
fn Split(comptime T: type, comptime kind: Kind, comptime lower: bool, comptime op: Trans, comptime unit: bool) type {
    const F = Float(T);
    const N = Num(T);
    const tr = comptime op.transposed();
    const cj = comptime op.conjA();
    return struct {
        const Self = @This();

        o: Ops(T),
        n: BlasLong,
        k: BlasLong,
        a: [*]T,
        lda: BlasLong,
        /// Unit stride.
        x: [*]T,
        /// The `gemv` buffers, `ge` elements each.
        g: [*]F,
        ge: usize,
        parts: *const threading.Partials(T),

        /// The rows of `op(A) * x` columns `c.start..c.end` reach.
        fn rows(s: *const Self, c: threading.Chunk) threading.Rows {
            if (tr) return .{ .lo = c.start, .hi = c.end };
            const nu: usize = @intCast(s.n);
            const k: usize = if (kind == .tbmv) @intCast(s.k) else nu;
            return if (lower) .{ .lo = c.start, .hi = @min(nu, c.end + k) } else .{ .lo = c.start -| k, .hi = c.end };
        }

        fn one(s: *const Self, c: threading.Chunk) void {
            const y = s.parts.take(c.index).ptr;
            const yb: BlasLong = @intCast(s.parts.rows[c.index].lo);
            const c0: BlasLong = @intCast(c.start);
            const c1: BlasLong = @intCast(c.end);
            if (kind == .trmv) s.full(y, yb, c0, c1, s.g + c.index * s.ge) else s.columns(y, yb, c0, c1);
        }

        /// `trmv_thread.c`: blocks of `dtb` columns from `c0`, each a `gemv`
        /// off the block and `axpy` or `dot` calls inside it.
        fn full(s: *const Self, y: [*]T, yb: BlasLong, c0: BlasLong, c1: BlasLong, g: [*]F) void {
            const o = s.o;
            const a = s.a;
            const lda = s.lda;
            const x = s.x;
            const m = s.n;
            var is = c0;
            while (is < c1) : (is += o.dtb) {
                const min_i = @min(c1 - is, o.dtb);
                if (!lower and is > 0) {
                    const pa = off(T, a, is * lda);
                    if (tr) o.gemv(op, is, min_i, 1, pa, lda, x, off(T, y, is - yb), g) else o.gemv(op, is, min_i, 1, pa, lda, off(T, x, is), off(T, y, -yb), g);
                }
                var i = is;
                while (i < is + min_i) : (i += 1) {
                    const col = off(T, a, i * lda);
                    const yi = at(T, y, i - yb);
                    const xi = at(T, x, i).*;
                    if (!lower and i > is) {
                        if (tr) yi.* = N.add(yi.*, o.dot(cj, i - is, off(T, col, is), off(T, x, is))) else o.axpy(cj, i - is, xi, off(T, col, is), off(T, y, is - yb));
                    }
                    yi.* = N.add(yi.*, if (unit) xi else N.mul(at(T, col, i).*, xi, cj));
                    if (lower and is + min_i > i + 1) {
                        const len = is + min_i - i - 1;
                        if (tr) yi.* = N.add(yi.*, o.dot(cj, len, off(T, col, i + 1), off(T, x, i + 1))) else o.axpy(cj, len, xi, off(T, col, i + 1), off(T, y, i + 1 - yb));
                    }
                }
                if (lower and m > is + min_i) {
                    const pa = off(T, a, is + min_i + is * lda);
                    if (tr) o.gemv(op, m - is - min_i, min_i, 1, pa, lda, off(T, x, is + min_i), off(T, y, is - yb), g) else o.gemv(op, m - is - min_i, min_i, 1, pa, lda, off(T, x, is), off(T, y, is + min_i - yb), g);
                }
            }
        }

        /// `tbmv_thread.c`, `tpmv_thread.c`: one column at a time, the
        /// off-diagonal part before the diagonal above it, after it below.
        fn columns(s: *const Self, y: [*]T, yb: BlasLong, c0: BlasLong, c1: BlasLong) void {
            const n = s.n;
            const x = s.x;
            var i = c0;
            while (i < c1) : (i += 1) {
                const room = if (lower) n - 1 - i else i;
                const len = if (kind == .tbmv) @min(room, s.k) else room;
                const col: BlasLong = if (kind == .tbmv) i * s.lda else if (lower) i * n - @divExact(i * (i - 1), 2) else @divExact(i * (i + 1), 2);
                const dp = off(T, s.a, col + (if (lower) 0 else if (kind == .tbmv) s.k else i));
                const yi = at(T, y, i - yb);
                const xi = at(T, x, i).*;
                if (!lower and len > 0) s.part(y, yb, yi, xi, len, off(T, dp, -len), i - len);
                yi.* = N.add(yi.*, if (unit) xi else N.mul(dp[0], xi, cj));
                if (lower and len > 0) s.part(y, yb, yi, xi, len, off(T, dp, 1), i + 1);
            }
        }

        /// Column `i`'s `len` off-diagonal elements at `ap`, rows `r0..`.
        fn part(s: *const Self, y: [*]T, yb: BlasLong, yi: *T, xi: T, len: BlasLong, ap: [*]T, r0: BlasLong) void {
            if (tr) yi.* = N.add(yi.*, s.o.dot(cj, len, ap, off(T, s.x, r0))) else s.o.axpy(cj, len, xi, ap, off(T, y, r0 - yb));
        }
    };
}

/// `?TRMV` (`interface/trmv.c`, `ztrmv.c`): `x := op(A) * x`, `A` upper or
/// lower triangular of `n`, column-major with leading dimension `lda`. The
/// other triangle (and the diagonal for `unit`) is not read. A negative
/// `incx` walks `x` from its end.
pub fn trmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, n: isize, a: [*]T, lda: isize, x: [*]T, incx: isize) Error!void {
    return public(allocator, io, T, .trmv, uplo, trans, diag, n, 0, a, lda, x, incx);
}

/// `?TRSV` (`interface/trsv.c`, `ztrsv.c`): solves `op(A) * y = x`, `y`
/// replacing `x`. No singularity check, as in the reference.
pub fn trsv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, n: isize, a: [*]T, lda: isize, x: [*]T, incx: isize) Error!void {
    return public(allocator, io, T, .trsv, uplo, trans, diag, n, 0, a, lda, x, incx);
}

/// `trsv` over `vt`'s kernels at unit stride, never routed to Accelerate:
/// LAPACK's `getrs` with one right-hand side (`lapack/getrs/getrs_single.c`),
/// which calls the level-3 drivers alike. The caller checks the arguments.
pub fn trsvTable(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, uplo: Uplo, trans: Trans, diag: Diag, n: usize, a: [*]const T, lda: usize, x: [*]T) (Allocator.Error || Io.Cancelable)!void {
    call(allocator, io, T, .trsv, vt, uplo, trans, diag, @intCast(n), 0, @constCast(a), @intCast(lda), x, 1) catch |e| switch (e) {
        error.OutOfMemory, error.Canceled => |f| return f,
        else => unreachable,
    };
}

/// `?TBMV` (`interface/tbmv.c`, `ztbmv.c`): `x := op(A) * x`, `A` triangular
/// band with `k` off-diagonals in LAPACK band storage: `A(i,j)` is
/// `a[k + i - j + j*lda]` (upper) or `a[i - j + j*lda]` (lower).
pub fn tbmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, n: isize, k: isize, a: [*]T, lda: isize, x: [*]T, incx: isize) Error!void {
    return public(allocator, io, T, .tbmv, uplo, trans, diag, n, k, a, lda, x, incx);
}

/// `?TBSV` (`interface/tbsv.c`, `ztbsv.c`): solves `op(A) * y = x` for band `A`.
pub fn tbsv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, n: isize, k: isize, a: [*]T, lda: isize, x: [*]T, incx: isize) Error!void {
    return public(allocator, io, T, .tbsv, uplo, trans, diag, n, k, a, lda, x, incx);
}

/// `?TPMV` (`interface/tpmv.c`, `ztpmv.c`): `x := op(A) * x`, `A` packed by
/// columns: `A(i,j)` is `ap[j*(j+1)/2 + i]` (upper, `i <= j`) or
/// `ap[j*(2n-j+1)/2 + i - j]` (lower, `i >= j`).
pub fn tpmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, n: isize, ap: [*]T, x: [*]T, incx: isize) Error!void {
    return public(allocator, io, T, .tpmv, uplo, trans, diag, n, 0, ap, 1, x, incx);
}

/// `?TPSV` (`interface/tpsv.c`, `ztpsv.c`): solves `op(A) * y = x` for packed `A`.
pub fn tpsv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, n: isize, ap: [*]T, x: [*]T, incx: isize) Error!void {
    return public(allocator, io, T, .tpsv, uplo, trans, diag, n, 0, ap, 1, x, incx);
}

// ---- tests -----------------------------------------------------------------

const check = @import("../_test/check.zig");
const guard = @import("../_test/guard.zig");
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;

const real = .{ f32, f64 } ++ check.x_real;
const complex = .{ Complex32, Complex64 } ++ check.x_complex;
const ops = [_]Trans{ .n, .t, .r, .c };
const sizes = [_]usize{ 1, 2, 3, 5, 8, 17, 70 };
const incs = [_]isize{ 1, 2, -1, -3 };
const max_n = 70;

const Layout = enum { full, band, packed_ };

fn layoutOf(comptime kind: Kind) Layout {
    return switch (kind) {
        .trmv, .trsv => .full,
        .tbmv, .tbsv => .band,
        .tpmv, .tpsv => .packed_,
    };
}

fn isSolve(comptime kind: Kind) bool {
    return kind == .trsv or kind == .tbsv or kind == .tpsv;
}

/// Where `A(i,j)` lives in the storage, for `i` and `j` inside the triangle.
fn index(layout: Layout, uplo: Uplo, n: usize, k: usize, lda: usize, i: usize, j: usize) usize {
    return switch (layout) {
        .full => i + j * lda,
        .band => (if (uplo == .u) k + i - j else i - j) + j * lda,
        .packed_ => if (uplo == .u) j * (j + 1) / 2 + i else j * (2 * n - j + 1) / 2 + i - j,
    };
}

fn inTriangle(layout: Layout, uplo: Uplo, k: usize, i: usize, j: usize) bool {
    if (uplo == .u) {
        if (i > j) return false;
        return layout != .band or j - i <= k;
    }
    if (i < j) return false;
    return layout != .band or i - j <= k;
}

/// `op(A)` as a dense row-major `n` by `n` matrix.
fn dense(comptime T: type, layout: Layout, uplo: Uplo, op: Trans, unit: bool, n: usize, k: usize, lda: usize, a: []const T, out: []T) void {
    const N = Num(T);
    const real_op = canon(T, op);
    for (0..n) |i| for (0..n) |j| {
        // (i,j) of op(A) is (r,c) of A
        const r = if (real_op.transposed()) j else i;
        const c = if (real_op.transposed()) i else j;
        var v: T = if (comptime isComplex(T)) .init(0, 0) else 0;
        if (inTriangle(layout, uplo, k, r, c)) {
            if (r == c and unit) {
                v = if (comptime isComplex(T)) .init(1, 0) else 1;
            } else {
                v = a[index(layout, uplo, n, k, lda, r, c)];
            }
            if (real_op.conjA()) v = N.conj(v, true);
        }
        out[i * n + j] = v;
    };
}

fn storageLen(layout: Layout, n: usize, lda: usize) usize {
    return switch (layout) {
        .full => lda * n,
        .band => lda * n,
        .packed_ => n * (n + 1) / 2,
    };
}

/// Physical index of logical element `i` of a vector of `n` at stride `inc`.
fn phys(i: usize, n: usize, inc: isize) usize {
    const a: usize = @intCast(@abs(inc));
    return if (inc > 0) i * a else (n - 1 - i) * a;
}

fn expectSame(comptime T: type, want: []const T, got: []const T, tol: Float(T)) !void {
    try check.expectCloseSlices(T, want, got, tol);
}

fn checkKind(comptime T: type, comptime kind: Kind, vt: *const Vtable) !void {
    const F = Float(T);
    const N = Num(T);
    const layout = comptime layoutOf(kind);
    const solve = comptime isSolve(kind);
    // scratch from the heap holds NaN, so a missing clear shows
    var poison: guard.Poison = .{ .child = std.testing.allocator };
    const gpa = poison.allocator();
    for (sizes) |n| for ([_]usize{ 0, 1, 3, 1000 }) |k_req| {
        if (layout != .band and k_req != 0) continue;
        const k = if (layout == .band) @min(k_req, n - 1) else 0;
        const lda = if (layout == .band) k + 1 + (n % 3) else n + (n % 3);
        for ([_]Uplo{ .u, .l }) |uplo| for (ops) |op| for (std.enums.values(Diag)) |diag| for (incs) |inc| {
            // static: at max_n a complex f80 frame would be ~330 KB
            const a = guard.operand(T, .a, max_n * (max_n + 3));
            check.fill(T, a[0..storageLen(layout, n, lda)], n * 31 + k + @backingInt(op));
            const alen = storageLen(layout, n, lda);
            // small off-diagonals and a dominant diagonal keep solves tame
            const sc: F = 1 / @sqrt(@as(F, @floatFromInt(n + 1)));
            for (a[0..alen]) |*v| v.* = if (comptime isComplex(T)) .init(v.re * sc, v.im * sc) else v.* * sc;
            for (0..n) |i| {
                const dg = &a[index(layout, uplo, n, k, lda, i, i)];
                dg.* = if (comptime isComplex(T)) .init(3 + @abs(dg.re), dg.im) else 3 + @abs(dg.*);
            }
            const d = guard.operand(T, .want, max_n * max_n);
            dense(T, layout, uplo, op, diag == .u, n, k, lda, a[0..alen], d[0 .. n * n]);

            const span = (n - 1) * @abs(inc) + 1;
            var x: [max_n * 3]T = undefined;
            check.fill(T, x[0..span], n + 17);
            var want = x;
            var xl: [max_n]T = undefined;
            for (0..n) |i| xl[i] = x[phys(i, n, inc)];
            var rl: [max_n]T = undefined;
            const zero: T = if (comptime isComplex(T)) .init(0, 0) else 0;
            if (!solve) {
                for (0..n) |i| {
                    var acc = zero;
                    for (0..n) |j| acc = N.add(acc, N.mul(d[i * n + j], xl[j], false));
                    rl[i] = acc;
                }
                for (0..n) |i| want[phys(i, n, inc)] = rl[i];
            }

            const nn: isize = @intCast(n);
            const kk: isize = @intCast(k);
            const ll: isize = @intCast(lda);
            const ii: isize = @intCast(inc);
            try call(gpa, check.io, T, kind, vt, uplo, op, diag, nn, kk, a.ptr, ll, &x, ii);
            const tol = check.tolerance(F, 4 * n + 4);
            if (!solve) {
                expectSame(T, want[0..span], x[0..span], tol) catch |e| return report(e, T, kind, uplo, op, diag, n, k, inc, vt.dtb_entries);
            } else {
                // op(A) * solution must give back the right-hand side
                for (0..n) |i| {
                    var acc = zero;
                    for (0..n) |j| acc = N.add(acc, N.mul(d[i * n + j], x[phys(j, n, inc)], false));
                    rl[i] = acc;
                }
                for (0..n) |i| {
                    check.expectClose(T, xl[i], rl[i], tol) catch |e| return report(e, T, kind, uplo, op, diag, n, k, inc, vt.dtb_entries);
                }
                // and the gaps between strided elements are untouched
                for (0..span) |p| {
                    var hit = false;
                    for (0..n) |i| hit = hit or phys(i, n, inc) == p;
                    if (!hit) try check.expectClose(T, want[p], x[p], 0);
                }
            }
        };
    };
}

fn report(e: anyerror, comptime T: type, comptime kind: Kind, uplo: Uplo, op: Trans, diag: Diag, n: usize, k: usize, inc: isize, dtb: i32) anyerror {
    check.say("{s} {s} {s} {s} {s} n={d} k={d} inc={d} dtb={d}: ", .{ @tagName(kind), @typeName(T), @tagName(uplo), @tagName(op), @tagName(diag), n, k, inc, dtb });
    return e;
}

const every_kind = std.enums.values(Kind);

/// `list`'s types through `kinds`, at the vtable's block size and a tiny one
/// that reaches the `gemv` paths.
fn Check(comptime list: anytype, comptime kinds: []const Kind) fn (*const Vtable) anyerror!void {
    return struct {
        fn run(vt: *const Vtable) anyerror!void {
            for ([_]i32{ vt.dtb_entries, 3 }) |dtb| {
                var v = vt.*;
                v.dtb_entries = dtb;
                inline for (list) |T| {
                    inline for (kinds) |kind| try checkKind(T, kind, &v);
                }
            }
        }
    }.run;
}

test "triangular routines match dense arithmetic with portable kernels" {
    var vt: Vtable = .{};
    portable.install(&vt);
    try Check(real ++ complex, every_kind)(&vt);
}

test "f80 without slots: triangular routines match dense arithmetic with portable kernels" {
    if (comptime check.f80_unslotted.len == 0) return error.SkipZigTest;
    var vt: Vtable = .{};
    portable.install(&vt);
    try Check(check.f80_unslotted, every_kind)(&vt);
}

// in three, so each fits a rig image (`build/qemu/machines.zig`'s `sparc64-debug`)
test "every linked core's kernels give the same triangular results, real" {
    try check.eachCore(Check(real, every_kind));
}

test "every linked core's kernels give the same triangular results, complex products" {
    try check.eachCore(Check(complex, &.{ .trmv, .tbmv, .tpmv }));
}

test "every linked core's kernels give the same triangular results, complex solves" {
    try check.eachCore(Check(complex, &.{ .trsv, .tbsv, .tpsv }));
}

test "triangular work space beyond the stack comes from the allocator" {
    var vt: Vtable = .{};
    portable.install(&vt);
    const n = 300;
    const gpa = std.testing.allocator;
    const a = try gpa.alloc(f64, n * n);
    defer gpa.free(a);
    const x = try gpa.alloc(f64, 2 * n);
    defer gpa.free(x);
    check.fill(f64, a, 1);
    for (a) |*v| v.* *= 0.05;
    for (0..n) |i| a[i + i * n] = 4;
    check.fill(f64, x, 2);
    const b = try gpa.dupe(f64, x);
    defer gpa.free(b);
    try call(gpa, check.io, f64, .trsv, &vt, .l, .n, .n, n, 0, a.ptr, n, x.ptr, 2);
    for (0..n) |i| {
        var acc: f64 = 0;
        for (0..i + 1) |j| acc += a[i + j * n] * x[2 * j];
        try check.expectClose(f64, b[2 * i], acc, 1e-9);
    }
    var failing = std.testing.FailingAllocator.init(gpa, .{ .fail_index = 0 });
    try std.testing.expectError(error.OutOfMemory, call(failing.allocator(), check.io, f64, .trmv, &vt, .u, .n, .u, n, 0, a.ptr, n, x.ptr, 2));
}

test "triangular argument errors, in the reference's order" {
    var vt: Vtable = .{};
    portable.install(&vt);
    const gpa = std.testing.allocator;
    var a: [16]f64 = @splat(1);
    var x: [4]f64 = @splat(1);
    const e = std.testing.expectError;
    try e(error.InvalidN, call(gpa, check.io, f64, .trmv, &vt, .u, .n, .u, -1, 0, &a, 0, &x, 0));
    try e(error.InvalidLda, call(gpa, check.io, f64, .trmv, &vt, .u, .n, .u, 4, 0, &a, 3, &x, 0));
    try e(error.InvalidIncx, call(gpa, check.io, f64, .trsv, &vt, .u, .n, .u, 4, 0, &a, 4, &x, 0));
    try e(error.InvalidK, call(gpa, check.io, f64, .tbmv, &vt, .u, .n, .u, 4, -1, &a, 4, &x, 1));
    try e(error.InvalidLda, call(gpa, check.io, f64, .tbsv, &vt, .u, .n, .u, 4, 2, &a, 2, &x, 1));
    try e(error.InvalidIncx, call(gpa, check.io, f64, .tpmv, &vt, .l, .t, .u, 4, 0, &a, 1, &x, 0));
    // n == 0 returns before touching anything, even a bad lda
    try call(gpa, check.io, f64, .trmv, &vt, .u, .n, .u, 0, 0, &a, 1, &x, 1);
    // `v` and `m` are Hermitian forms, `o` to `d` conjugate `x`: neither applies
    try e(error.InvalidUplo, call(gpa, check.io, f64, .trmv, &vt, .m, .o, .u, -1, 0, &a, 0, &x, 0));
    try e(error.InvalidTrans, call(gpa, check.io, f64, .tpsv, &vt, .l, .d, .u, -1, 0, &a, 0, &x, 0));
}
