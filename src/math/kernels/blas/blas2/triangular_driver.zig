//! The triangular drivers (`triangular.zig`): `driver/level2/{,z}{trmv,trsv,tbmv,tbsv,tpmv,tpsv}_{U,L}.c`.

const types = @import("../types.zig");
const util = @import("../util.zig");
const Num = @import("ops.zig").Num;
const Ops = @import("ops.zig").Ops;

const BlasLong = types.BlasLong;
const Trans = types.Trans;
const Float = types.Float;
const at = util.at;
const off = util.off;

/// Which reference driver file: `_U` or `_L`.
pub const File = enum { u, l };

/// Shared by every driver: `m` is `n` for the full and packed ones.
pub fn Drive(comptime T: type, comptime file: File, comptime op: Trans, comptime unit: bool) type {
    const F = Float(T);
    const N = Num(T);
    const tr = comptime op.transposed();
    const cj = comptime op.conjA();
    return struct {
        /// `trmv_{U,L}.c`: `b := op(A) * b`.
        pub fn trmv(o: Ops(T), m: BlasLong, a: [*]T, lda: BlasLong, b: [*]T, g: [*]F) void {
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
        pub fn trsv(o: Ops(T), m: BlasLong, a: [*]T, lda: BlasLong, b: [*]T, g: [*]F) void {
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
        pub fn tbmv(o: Ops(T), n: BlasLong, k: BlasLong, a0: [*]T, lda: BlasLong, b: [*]T) void {
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
        pub fn tbsv(o: Ops(T), n: BlasLong, k: BlasLong, a0: [*]T, lda: BlasLong, b: [*]T) void {
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
        pub fn tpmv(o: Ops(T), m: BlasLong, a0: [*]T, b: [*]T) void {
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
        pub fn tpsv(o: Ops(T), m: BlasLong, a0: [*]T, b: [*]T) void {
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
