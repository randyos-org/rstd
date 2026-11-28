//! Reordering a Schur form (MATLAB's `ordschur`), by `trsen`: a logical
//! select, a region keyword, or clusters (`reorder_select.zig`).

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
const Schur = schur_form.Schur;

pub const OrdschurError = reorder_select.OrderError;
pub const Region = reorder_select.Region;
pub const Select = reorder_select.Select;

/// The Schur form `A = U T Uᴴ` reordered as `select` asks (MATLAB's `[US,
/// TS] = ordschur(U, T, select)`, `keyword` or `clusters`), by `trsen`;
/// clusters by one `trsen` per index above the lowest.
pub fn ordschur(allocator: Allocator, io: Io, u: anytype, t: anytype, select: Select) OrdschurError!Schur(ElementOf(@TypeOf(t))) {
    const T = ElementOf(@TypeOf(t));
    const cu = constOf(u);
    const ct = constOf(t);
    if (ct.rows != ct.cols) return error.NotSquare;
    const n = ct.rows;
    if (cu.rows != n or cu.cols != n) return error.DimensionMismatch;
    try reorder_select.checkLength(select, n);
    const rt = try copyOf(T, allocator, ct);
    errdefer rt.deinit(allocator);
    const ru = try copyOf(T, allocator, cu);
    errdefer ru.deinit(allocator);
    const form: Form(T) = .{ .allocator = allocator, .io = io, .n = n, .t = rt, .u = ru };
    try reorder_select.apply(allocator, form, select);
    return .{ .u = ru, .t = rt };
}

/// The Schur form as `reorder_select.apply` drives it.
fn Form(comptime T: type) type {
    return struct {
        allocator: Allocator,
        io: Io,
        n: usize,
        t: Mat(T),
        u: Mat(T),

        const Self = @This();

        pub fn opensBlock(f: Self, k: usize) bool {
            if (comptime isComplex(T)) return false;
            return k + 1 < f.n and f.t.data[k + 1 + k * f.n] != 0;
        }

        /// A 2 x 2 block's `[a b; c d]` pair has real part `(a + d) / 2`
        /// and squared modulus `a d - b c`.
        pub fn eigen(f: Self, k: usize, size: usize) reorder_select.Eigen(Real(T)) {
            const n = f.n;
            const t = f.t.data;
            if (comptime isComplex(T)) {
                const x = t[k + k * n];
                return .{ .re = x.re, .mod2 = x.re * x.re + x.im * x.im };
            }
            if (size == 2) {
                const a = t[k + k * n];
                const b = t[k + (k + 1) * n];
                const c = t[k + 1 + k * n];
                const d = t[k + 1 + (k + 1) * n];
                return .{ .re = (a + d) / 2, .mod2 = a * d - b * c };
            }
            const re = t[k + k * n];
            return .{ .re = re, .mod2 = re * re };
        }

        /// One `trsen` in place.
        pub fn reorder(f: Self, mask: []const bool) OrdschurError!void {
            const n = f.n;
            const w = try f.allocator.alloc(T, if (comptime isComplex(T)) n else 2 * n);
            defer f.allocator.free(w);
            const ev: lapack.Eigenvalues(T) = if (comptime isComplex(T)) w else .{ .re = w[0..n], .im = w[n..] };
            _ = lapack.trsen(f.allocator, f.io, T, .none, .update, mask, n, f.t.data, n, f.u.data, n, ev) catch |e| return switch (e) {
                error.SwapRejected => error.EigenvaluesTooClose,
                else => ops.narrow(OrdschurError, e),
            };
        }
    };
}
