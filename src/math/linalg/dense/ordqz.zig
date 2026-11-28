//! Reordering a generalized Schur form (MATLAB's `ordqz`), by `tgsen`: a
//! logical select, a region keyword, or clusters (`reorder_select.zig`).

const std = @import("std");
const lapack = @import("../../kernels/lapack/__root__.zig");
const matrix = @import("../matrix.zig");
const vector = @import("../vector.zig");
const ops = @import("../../kernels/lapack/ops.zig");
const util = @import("util.zig");
const schur_form = @import("schur_form.zig");
const reorder_select = @import("reorder_select.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = matrix.Mat;
const isComplex = vector.isComplex;
const Real = vector.Real;
const ElementOf = util.ElementOf;
const constOf = util.constOf;
const copyOf = util.copyOf;
const Qz = schur_form.Qz;
const Select = reorder_select.Select;

pub const OrdqzError = reorder_select.OrderError;

/// `qz`'s form `Q A Z = AA`, `Q B Z = BB` reordered as `select` asks
/// (MATLAB's `[AAS, BBS, QS, ZS] = ordqz(AA, BB, Q, Z, select)`, `keyword`
/// or `clusters`, the eigenvalues `ordeig(AA, BB)`), by `tgsen`; clusters by
/// one `tgsen` per index above the lowest. `BBS`'s diagonal is made
/// nonnegative (real) or real and nonnegative (complex), as `tgsen`'s.
pub fn ordqz(allocator: Allocator, io: Io, aa: anytype, bb: anytype, q: anytype, z: anytype, select: Select) OrdqzError!Qz(ElementOf(@TypeOf(aa))) {
    const T = ElementOf(@TypeOf(aa));
    const ca = constOf(aa);
    const cb = constOf(bb);
    const cq = constOf(q);
    const cz = constOf(z);
    if (ca.rows != ca.cols or cb.rows != cb.cols) return error.NotSquare;
    const n = ca.rows;
    for ([_]@TypeOf(ca){ cb, cq, cz }) |m| if (m.rows != n or m.cols != n) return error.DimensionMismatch;
    try reorder_select.checkLength(select, n);
    const ra = try copyOf(T, allocator, ca);
    errdefer ra.deinit(allocator);
    const rb = try copyOf(T, allocator, cb);
    errdefer rb.deinit(allocator);
    const rq = try copyOf(T, allocator, cq);
    errdefer rq.deinit(allocator);
    const rz = try copyOf(T, allocator, cz);
    errdefer rz.deinit(allocator);
    // LAPACK's `Q` is MATLAB's `Qᴴ`, both ways
    util.adjointInPlace(T, rq);
    const form: Form(T) = .{ .allocator = allocator, .io = io, .n = n, .a = ra, .b = rb, .q = rq, .z = rz };
    try reorder_select.apply(allocator, form, select);
    util.adjointInPlace(T, rq);
    return .{ .aa = ra, .bb = rb, .q = rq, .z = rz };
}

/// The generalized Schur form as `reorder_select.apply` drives it.
fn Form(comptime T: type) type {
    return struct {
        allocator: Allocator,
        io: Io,
        n: usize,
        a: Mat(T),
        b: Mat(T),
        q: Mat(T),
        z: Mat(T),

        const Self = @This();
        const F = Real(T);

        pub fn opensBlock(f: Self, k: usize) bool {
            if (comptime isComplex(T)) return false;
            return k + 1 < f.n and f.a.data[k + 1 + k * f.n] != 0;
        }

        /// `λ = a / b` for a 1 x 1 block (`b = 0`: an infinite or NaN real
        /// part, an infinite modulus unless `a` is 0 too); a 2 x 2 block's
        /// pair from `det(A - λ B) = 0`.
        pub fn eigen(f: Self, k: usize, size: usize) reorder_select.Eigen(F) {
            const n = f.n;
            const a = f.a.data;
            const b = f.b.data;
            if (comptime !isComplex(T)) if (size == 2) return pairEigen(F, a, b, n, k);
            const x = a[k + k * n];
            const y = b[k + k * n];
            if (comptime isComplex(T)) {
                if (y.re == 0 and y.im == 0) return .{ .re = x.re / @as(F, 0), .mod2 = (x.re * x.re + x.im * x.im) / @as(F, 0) };
                const l = x.div(y);
                return .{ .re = l.re, .mod2 = l.re * l.re + l.im * l.im };
            }
            const l = x / y;
            return .{ .re = l, .mod2 = l * l };
        }

        /// A real 2 x 2 block's pair: `det(B) λ² - s λ + det(A) = 0`.
        fn pairEigen(comptime R: type, a: []const R, b: []const R, n: usize, k: usize) reorder_select.Eigen(R) {
            const a11 = a[k + k * n];
            const a12 = a[k + (k + 1) * n];
            const a21 = a[k + 1 + k * n];
            const a22 = a[k + 1 + (k + 1) * n];
            const b11 = b[k + k * n];
            const b12 = b[k + (k + 1) * n];
            const b21 = b[k + 1 + k * n];
            const b22 = b[k + 1 + (k + 1) * n];
            // the roots a conjugate pair: sum `s / det(B)`, product `det(A) / det(B)`
            const db = b11 * b22 - b12 * b21;
            const s = a11 * b22 + a22 * b11 - a12 * b21 - a21 * b12;
            return .{ .re = s / (2 * db), .mod2 = (a11 * a22 - a12 * a21) / db };
        }

        /// One `tgsen` in place.
        pub fn reorder(f: Self, mask: []const bool) OrdqzError!void {
            const n = f.n;
            const w = try f.allocator.alloc(T, 3 * n);
            defer f.allocator.free(w);
            const ev: lapack.Eigenvalues(T) = if (comptime isComplex(T)) w[0..n] else .{ .re = w[0..n], .im = w[n .. 2 * n] };
            _ = lapack.tgsen(f.allocator, f.io, T, .none, .update, .update, mask, n, f.a.data, n, f.b.data, n, ev, w[2 * n ..], f.q.data, n, f.z.data, n) catch |e| return switch (e) {
                error.SwapRejected => error.EigenvaluesTooClose,
                else => ops.narrow(OrdqzError, e),
            };
        }
    };
}
