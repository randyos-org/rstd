//! Strided vector views, the vectors the public BLAS and LAPACK routines
//! take, and the error pieces their exact error sets are built from
//! (design.md, "Planned: a Zig-idiomatic interface"). A view holds its slice,
//! element count and stride: the count is never inferred, so spare elements
//! past the vector cannot change how many a routine works on. Matrices stay
//! a slice, sizes and a leading dimension.

const std = @import("std");

/// A slice too short for the sizes, strides and leading dimensions given.
pub const SliceError = error{SliceTooShort};
/// A stride of zero, which no view takes.
pub const StrideError = error{ZeroStride};
/// Vectors whose element counts should agree and don't.
pub const LengthError = error{LengthMismatch};
/// A leading dimension below the matrix's row count, or below 1.
pub const LeadingDimensionError = error{LeadingDimensionTooSmall};
/// What `init` returns.
pub const InitError = SliceError || StrideError;

/// `n` elements of `T` in `data`, `stride` apart. A negative stride walks
/// `data` backward: element 0 is the last in memory, as the reference BLAS
/// has it. `data.ptr`, `n` and `stride` are the kernels' pointer, size and
/// increment unchanged.
pub fn Strided(comptime T: type) type {
    return View(T, false);
}

/// `Strided(T)` that only reads.
pub fn ConstStrided(comptime T: type) type {
    return View(T, true);
}

fn View(comptime T: type, comptime read_only: bool) type {
    return struct {
        const Self = @This();
        pub const Element = T;
        const Slice = if (read_only) []const T else []T;

        data: Slice,
        n: usize,
        stride: isize,

        /// All of `s`, contiguous.
        pub fn of(s: Slice) Self {
            return .{ .data = s, .n = s.len, .stride = 1 };
        }

        /// `n` elements of `s`, `stride` apart; `s` must hold
        /// `(n - 1) * |stride| + 1` of them.
        pub fn init(s: Slice, n: usize, stride: isize) InitError!Self {
            if (stride == 0) return error.ZeroStride;
            if (n > 0) {
                const span = std.math.mul(usize, n - 1, @abs(stride)) catch return error.SliceTooShort;
                if (span >= s.len) return error.SliceTooShort;
            }
            return .{ .data = s, .n = n, .stride = stride };
        }

        /// Where element `i` is in `data`; writes go through it. Asserts `i < n`
        /// (callers index within `n`).
        pub fn index(v: Self, i: usize) usize {
            std.debug.assert(i < v.n);
            const step: usize = @abs(v.stride);
            return if (v.stride > 0) i * step else (v.n - 1 - i) * step;
        }

        pub fn get(v: Self, i: usize) T {
            return v.data[v.index(i)];
        }

        /// The same elements, read-only.
        pub fn toConst(v: Self) ConstStrided(T) {
            return .{ .data = v.data, .n = v.n, .stride = v.stride };
        }

        /// Whether `v` and `w` have as many elements.
        pub fn sameLength(v: Self, w: anytype) LengthError!void {
            if (v.n != w.n) return error.LengthMismatch;
        }
    };
}

// ---- tests ------------------------------------------------------------------

const testing = std.testing;

test "of: the whole slice, contiguous" {
    var a = [_]f64{ 1, 2, 3 };
    const v = Strided(f64).of(&a);
    try testing.expectEqual(3, v.n);
    try testing.expectEqual(1, v.stride);
    v.data[v.index(1)] = 5;
    try testing.expectEqual(5, a[1]);
    try testing.expectEqual(3, v.toConst().get(2));
}

test "init: positive and negative strides, element 0 last in memory for a negative one" {
    const a = [_]f32{ 0, 1, 2, 3, 4, 5, 6 };
    const up = try ConstStrided(f32).init(&a, 3, 3);
    try testing.expectEqual([3]f32{ 0, 3, 6 }, [3]f32{ up.get(0), up.get(1), up.get(2) });
    const down = try ConstStrided(f32).init(&a, 3, -3);
    try testing.expectEqual([3]f32{ 6, 3, 0 }, [3]f32{ down.get(0), down.get(1), down.get(2) });
}

test "init: the slice must reach the last element" {
    const a = [_]f64{ 0, 1, 2, 3, 4, 5 };
    try testing.expectError(error.SliceTooShort, ConstStrided(f64).init(&a, 3, 3));
    try testing.expectError(error.SliceTooShort, ConstStrided(f64).init(&a, 7, 1));
    try testing.expectError(error.SliceTooShort, ConstStrided(f64).init(&a, std.math.maxInt(usize), 2));
    _ = try ConstStrided(f64).init(&a, 2, 5);
    _ = try ConstStrided(f64).init(&a, 6, -1);
}

test "init: a zero stride is refused, an empty vector fits any slice" {
    const a = [_]f64{1};
    try testing.expectError(error.ZeroStride, ConstStrided(f64).init(&a, 1, 0));
    const empty = try ConstStrided(f64).init(&.{}, 0, 4);
    try testing.expectEqual(0, empty.n);
}

test "sameLength: counts must agree" {
    const a = [_]f64{ 1, 2, 3, 4 };
    const x = ConstStrided(f64).of(&a);
    try x.sameLength(try ConstStrided(f64).init(&a, 4, -1));
    try testing.expectError(error.LengthMismatch, x.sameLength(try ConstStrided(f64).init(&a, 2, 2)));
}
