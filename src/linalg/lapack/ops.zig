//! What the LAPACK routines share: a column-major matrix view, element
//! arithmetic over a real or `Complex` `T`, and error narrowing. The argument
//! checks, vtable and thread rules are `interface.zig`'s, the BLAS calls
//! `blas_calls.zig`'s, the LAPACK auxiliaries `auxiliary.zig`'s.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
pub const Float = types.Float;
pub const isComplex = types.isComplex;
pub const Trans = types.Trans;

/// What the internal routines return.
pub const Error = Allocator.Error || Io.Cancelable;

/// Tests only: caps the iterative solvers (`bidiag_svd.qrAlgorithm`, both
/// `lahqr`s) at this many iterations, so `NoConvergence` can be reached.
pub var test_iter_limit: ?usize = null;

/// `limit`, or `test_iter_limit` when set.
pub fn iterLimit(limit: usize) usize {
    return test_iter_limit orelse limit;
}

pub const requireFloat = util.requireFloat;
pub const requireComplex = util.requireComplex;

/// A column-major view: element `(i, j)` at `ptr[i + j * ld]`.
pub fn Mat(comptime T: type) type {
    return struct {
        ptr: [*]T,
        rows: usize,
        cols: usize,
        ld: usize,

        const Self = @This();

        pub fn at(m: Self, i: usize, j: usize) *T {
            return &m.ptr[i + j * m.ld];
        }

        pub fn get(m: Self, i: usize, j: usize) T {
            return m.ptr[i + j * m.ld];
        }

        /// The `r x c` block at `(i, j)`.
        pub fn sub(m: Self, i: usize, j: usize, r: usize, c: usize) Self {
            return .{ .ptr = m.ptr + i + j * m.ld, .rows = r, .cols = c, .ld = m.ld };
        }

        pub fn col(m: Self, j: usize) [*]T {
            return m.ptr + j * m.ld;
        }

        pub fn fill(m: Self, v: T) void {
            for (0..m.cols) |j| @memset(m.col(j)[0..m.rows], v);
        }

        /// Zero, with ones on the diagonal.
        pub fn identity(m: Self) void {
            m.fill(zero(T));
            for (0..@min(m.rows, m.cols)) |i| m.at(i, i).* = one(T);
        }

        pub fn copyFrom(m: Self, src: Self) void {
            for (0..m.cols) |j| @memcpy(m.col(j)[0..m.rows], src.col(j)[0..m.rows]);
        }
    };
}

/// An owned `rows x cols` matrix, `ld = max(rows, 1)`.
pub fn alloc(comptime T: type, al: Allocator, rows: usize, cols: usize) Allocator.Error!Mat(T) {
    const ld = @max(rows, 1);
    const buf = try al.alloc(T, ld * cols);
    return .{ .ptr = buf.ptr, .rows = rows, .cols = cols, .ld = ld };
}

pub fn free(comptime T: type, al: Allocator, m: Mat(T)) void {
    al.free(m.ptr[0 .. m.ld * m.cols]);
}

// ---- element arithmetic -----------------------------------------------------------

pub fn scalar(comptime T: type, v: Float(T)) T {
    return if (comptime isComplex(T)) .init(v, 0) else v;
}

pub fn zero(comptime T: type) T {
    return scalar(T, 0);
}

pub fn one(comptime T: type) T {
    return scalar(T, 1);
}

pub fn re(comptime T: type, x: T) Float(T) {
    return if (comptime isComplex(T)) x.re else x;
}

pub fn add(comptime T: type, x: T, y: T) T {
    return if (comptime isComplex(T)) .init(x.re + y.re, x.im + y.im) else x + y;
}

pub fn sub(comptime T: type, x: T, y: T) T {
    return if (comptime isComplex(T)) .init(x.re - y.re, x.im - y.im) else x - y;
}

pub fn neg(comptime T: type, x: T) T {
    return if (comptime isComplex(T)) .init(-x.re, -x.im) else -x;
}

pub fn mul(comptime T: type, x: T, y: T) T {
    if (comptime !isComplex(T)) return x * y;
    return .init(x.re * y.re - x.im * y.im, x.re * y.im + x.im * y.re);
}

pub fn conj(comptime T: type, x: T) T {
    return if (comptime isComplex(T)) .init(x.re, -x.im) else x;
}

pub fn mulReal(comptime T: type, x: T, r: Float(T)) T {
    return if (comptime isComplex(T)) .init(x.re * r, x.im * r) else x * r;
}

/// `conj(x) * y`.
pub fn mulConj(comptime T: type, x: T, y: T) T {
    if (comptime !isComplex(T)) return x * y;
    return .init(x.re * y.re + x.im * y.im, x.re * y.im - x.im * y.re);
}

/// `x / d`, part by part: a real `d`, so each part's quotient is exact.
pub fn divReal(comptime T: type, x: T, d: Float(T)) T {
    return if (comptime isComplex(T)) .init(x.re / d, x.im / d) else x / d;
}

/// `x / y`, complex by Smith's algorithm (`util.cdiv`), as `ZLADIV`.
pub fn div(comptime T: type, x: T, y: T) T {
    if (comptime !isComplex(T)) return x / y;
    const q = util.cdiv(Float(T), x.re, x.im, y.re, y.im);
    return .init(q[0], q[1]);
}

/// `|x|`, by `hypot` for complex.
pub fn abs(comptime T: type, x: T) Float(T) {
    return if (comptime isComplex(T)) std.math.hypot(x.re, x.im) else @abs(x);
}

/// The principal square root. Complex: std's for `f32` and `f64`, which has
/// no other type (issue #202), and the same musl `csqrt` for the rest (T124).
pub fn sqrt(comptime T: type, x: T) T {
    if (comptime !isComplex(T)) return @sqrt(x);
    const F = Float(T);
    if (F == f32 or F == f64) return std.math.complex.sqrt(x);
    // past this, x + hypot(x, y) may overflow
    const threshold = comptime std.math.floatMax(F) / (1 + std.math.sqrt2);
    var xr = x.re;
    var xi = x.im;
    if (xr == 0 and xi == 0) return .init(0, xi);
    if (std.math.isInf(xi)) return .init(std.math.inf(F), xi);
    if (std.math.isNan(xr)) return .init(xr, (xi - xi) / (xi - xi));
    if (std.math.isInf(xr)) {
        return if (std.math.signbit(xr)) .init(@abs(xi - xi), std.math.copysign(xr, xi)) else .init(xr, std.math.copysign(xi - xi, xi));
    }
    const scale = @abs(xr) >= threshold or @abs(xi) >= threshold;
    if (scale) {
        xr *= 0.25;
        xi *= 0.25;
    }
    var r: T = undefined;
    if (xr >= 0) {
        const t = @sqrt((xr + std.math.hypot(xr, xi)) * 0.5);
        r = .init(t, xi / (2 * t));
    } else {
        const t = @sqrt((-xr + std.math.hypot(xr, xi)) * 0.5);
        r = .init(@abs(xi) / (2 * t), std.math.copysign(t, xi));
    }
    return if (scale) .init(2 * r.re, 2 * r.im) else r;
}

/// `|re| + |im|` (faer's `abs1`, the reference's `CABS1`); `|x|` for real.
pub fn abs1(comptime T: type, x: T) Float(T) {
    return if (comptime isComplex(T)) @abs(x.re) + @abs(x.im) else @abs(x);
}

pub fn abs2(comptime T: type, x: T) Float(T) {
    return if (comptime isComplex(T)) x.re * x.re + x.im * x.im else x * x;
}

/// `1 / x`, complex by Smith's algorithm (`util.cdiv`).
pub fn recip(comptime T: type, x: T) T {
    if (comptime !isComplex(T)) return 1 / x;
    const q = util.cdiv(Float(T), 1, 0, x.re, x.im);
    return .init(q[0], q[1]);
}

pub fn isZero(comptime T: type, x: T) bool {
    return if (comptime isComplex(T)) x.re == 0 and x.im == 0 else x == 0;
}

pub fn isFinite(comptime T: type, x: T) bool {
    return if (comptime isComplex(T)) std.math.isFinite(x.re) and std.math.isFinite(x.im) else std.math.isFinite(x);
}

/// The largest `|x(i)|` of a real `x`, 0 when empty.
pub fn absMax(comptime F: type, x: []const F) F {
    var m: F = 0;
    for (x) |v| m = @max(m, @abs(v));
    return m;
}

/// `x *= 1 / max`, by division where that reciprocal overflows (faer
/// multiplies, which a denormal `max` makes infinite; issue #252).
pub fn scaleDown(comptime F: type, x: []F, max: F) void {
    if (max >= std.math.floatMin(F)) {
        const inv = 1 / max;
        for (x) |*v| v.* *= inv;
    } else for (x) |*v| {
        v.* /= max;
    }
}

/// Columns `i` and `j` of `m` swapped.
pub fn swapCols(comptime T: type, m: Mat(T), i: usize, j: usize) void {
    for (0..m.rows) |r| std.mem.swap(T, m.at(r, i), m.at(r, j));
}

/// `x / |x|` (1 for 0), complex `x`, scaled first so `|x|` cannot overflow.
/// Faer multiplies by `1 / norm1`, infinite for a denormal `x`; that divides.
pub fn normalized(comptime T: type, x: T) T {
    if (isZero(T, x)) return one(T);
    const norm1 = @max(@abs(x.re), @abs(x.im));
    const y: T = if (norm1 >= std.math.floatMin(Float(T))) mulReal(T, x, 1 / norm1) else .init(x.re / norm1, x.im / norm1);
    return mulReal(T, y, 1 / abs(T, y));
}

/// `e` narrowed to `E`, for a callee given arguments it accepts: the errors
/// `E` names pass (such as `OutOfMemory`, `Canceled` or `Singular`); any
/// other is `unreachable`.
pub fn narrow(comptime E: type, e: anyerror) E {
    inline for (@typeInfo(E).error_set.?) |x| {
        if (e == @field(anyerror, x.name)) return @errorCast(e);
    }
    unreachable;
}

test {
    _ = Mat(f64);
}

const testing = std.testing;

test "ops sqrt: the complex square root in every float type agrees with f64's (T124)" {
    const C64 = std.math.Complex(f64);
    const cases = [_]C64{ .init(5, 3), .init(-5, 3), .init(-5, -3), .init(0, 2), .init(-4, 0), .init(1e300, 1e300), .init(-1e-300, 1e-310) };
    inline for (.{ f80, f128 }) |F| {
        const C = std.math.Complex(F);
        for (cases) |z| {
            const want = std.math.complex.sqrt(z);
            const got = sqrt(C, .init(z.re, z.im));
            try testing.expectApproxEqRel(want.re, @as(f64, @floatCast(got.re)), 4 * std.math.floatEps(f64));
            try testing.expectApproxEqRel(want.im, @as(f64, @floatCast(got.im)), 4 * std.math.floatEps(f64));
        }
        const inf = std.math.inf(F);
        const r = sqrt(C, .init(-inf, 1));
        try testing.expect(r.re == 0 and r.im == inf);
        try testing.expect(std.math.isNan(sqrt(C, .init(std.math.nan(F), 1)).re));
    }
}
