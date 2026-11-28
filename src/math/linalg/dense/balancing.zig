//! Balancing for eigenvalue accuracy (MATLAB's `balance`), by `gebal` and
//! `gebak`.

const std = @import("std");
const lapack = @import("../../kernels/lapack/__root__.zig");
const matrix = @import("../matrix.zig");
const vector = @import("../vector.zig");
const ops = @import("../../kernels/lapack/ops.zig");
const util = @import("util.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = matrix.Mat;
const Real = vector.Real;
const ElementOf = util.ElementOf;
const constOf = util.constOf;
const copyOf = util.copyOf;

/// `balance`'s permutations and scalings (MATLAB's default) or scalings only
/// (`"noperm"`).
pub const BalanceMode = enum { both, noperm };

/// `balance`'s result: `B = T⁻¹ A T` with `T = I(:, p) diag(s)` up to the
/// order MATLAB's `[s, p, B] = balance(A)` gives them: `T(i, p[i]) = s[i]`.
pub fn Balance(comptime T: type) type {
    return struct {
        b: Mat(T),
        /// Powers of 2.
        s: []Real(T),
        p: []usize,

        const Self = @This();

        pub fn deinit(f: Self, allocator: Allocator) void {
            f.b.deinit(allocator);
            allocator.free(f.s);
            allocator.free(f.p);
        }

        /// `T` (MATLAB's `[T, B] = balance(A)`).
        pub fn t(f: Self, allocator: Allocator) Allocator.Error!Mat(T) {
            const n = f.s.len;
            const r = try Mat(T).zeros(allocator, n, n);
            for (f.p, f.s, 0..) |c, x, i| r.data[i + c * n] = ops.scalar(T, x);
            return r;
        }
    };
}

pub const BalanceError = error{ NotSquare, NanInA } || Allocator.Error || Io.Cancelable;

/// `B = T⁻¹ A T` (MATLAB's `balance`), its rows and columns as near in norm
/// as powers of 2 make them, isolated eigenvalues permuted out first unless
/// `mode` is `.noperm`. `error.NanInA` for a NaN in `a`.
pub fn balance(allocator: Allocator, io: Io, a: anytype, mode: BalanceMode) BalanceError!Balance(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    if (ca.rows != ca.cols) return error.NotSquare;
    const n = ca.rows;
    const ld = @max(n, 1);
    const job: lapack.Balancing = if (mode == .both) .both else .scale;
    const b = try copyOf(T, allocator, ca);
    errdefer b.deinit(allocator);
    const scale = try allocator.alloc(Real(T), n);
    defer allocator.free(scale);
    const range = lapack.gebal(io, T, job, n, b.data, ld, scale) catch |e| switch (e) {
        error.NanInA => return error.NanInA,
        error.Canceled => return error.Canceled,
        else => unreachable, // shapes are ours
    };
    // `T = P D`: `gebak` on the identity
    const t = try Mat(T).zeros(allocator, n, n);
    defer t.deinit(allocator);
    for (0..n) |i| t.data[i + i * n] = ops.one(T);
    lapack.gebak(io, T, job, .r, n, range.lo, range.hi, scale, n, t.data, ld) catch |e| switch (e) {
        error.Canceled => return error.Canceled,
        else => unreachable, // `gebal`'s own block and scale
    };
    const s = try allocator.alloc(Real(T), n);
    errdefer allocator.free(s);
    const p = try allocator.alloc(usize, n);
    errdefer allocator.free(p);
    for (0..n) |i| {
        const c = for (0..n) |j| {
            if (!ops.isZero(T, t.data[i + j * n])) break j;
        } else unreachable; // a permutation times a nonsingular diagonal
        p[i] = c;
        s[i] = ops.re(T, t.data[i + c * n]);
    }
    return .{ .b = b, .s = s, .p = p };
}
