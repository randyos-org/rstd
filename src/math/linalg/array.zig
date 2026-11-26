//! N-d arrays beside `Mat(T)` (src/math/linalg/design.md, "Planned: N-d
//! arrays"): `StaticArray` (comptime shape) and `Array` with its views,
//! column-major (the first index varies fastest), up to `max_rank`. Stubs
//! for now: the types and `Array`'s storage are real; the view and
//! element-wise bodies return `error.NotImplemented`, which leaves each
//! error set when its body lands, and their tests skip until then.

const std = @import("std");
const matrix = @import("matrix.zig");

const Allocator = std.mem.Allocator;

pub const max_rank = 8;

/// A shape: `dims[0..rank]`, the rest unused.
pub const Shape = struct {
    dims: [max_rank]usize = @splat(1),
    rank: u8,

    pub fn init(dims: []const usize) error{RankTooHigh}!Shape {
        if (dims.len > max_rank) return error.RankTooHigh;
        var s: Shape = .{ .rank = @intCast(dims.len) };
        @memcpy(s.dims[0..dims.len], dims);
        return s;
    }

    pub fn count(s: Shape) usize {
        var n: usize = 1;
        for (s.dims[0..s.rank]) |d| n *= d;
        return n;
    }
};

pub const Range = struct { start: usize, end: usize, step: isize = 1 };
pub const Index = union(enum) { all, at: usize, range: Range };

/// A comptime shape (`.{ 4, 4, 3 }`), a value with no allocator.
pub fn StaticArray(comptime dims: anytype, comptime T: type) type {
    const ds: [dims.len]usize = dims;
    if (ds.len > max_rank) @compileError("rstd.math.linalg: StaticArray of rank past max_rank");
    var n: usize = 1;
    for (ds) |d| n *= d;
    const count = n;
    return struct {
        data: [count]T,

        pub const shape: [ds.len]usize = ds;
    };
}

/// A borrowed array: `data[Σ index[k] strides[k]]`. `ConstArrayView` is its
/// read-only form.
pub fn ArrayView(comptime T: type) type {
    return ViewOf(T, false);
}

pub fn ConstArrayView(comptime T: type) type {
    return ViewOf(T, true);
}

pub const SliceError = error{ IndexOutOfRange, RankMismatch, NotImplemented };
pub const PermuteError = error{ InvalidPermutation, NotImplemented };
pub const ReshapeError = error{ CountMismatch, NotContiguous, NotImplemented };
pub const BroadcastError = error{ RankMismatch, DimensionMismatch };
pub const AsMatrixError = error{ RankMismatch, NotColumnMajor, NotImplemented };

fn ViewOf(comptime T: type, comptime read_only: bool) type {
    return struct {
        data: if (read_only) []const T else []T,
        shape: Shape,
        strides: [max_rank]isize,

        const Self = @This();

        /// A view of `indices` (one per dimension), no copy.
        pub fn slice(v: Self, indices: []const Index) SliceError!Self {
            if (indices.len != v.shape.rank) return error.RankMismatch;
            return error.NotImplemented;
        }

        /// MATLAB's `permute`: dimension `k` of the result is `order[k]`.
        pub fn permute(v: Self, order: []const u8) PermuteError!Self {
            _ = v;
            _ = order;
            return error.NotImplemented;
        }

        pub fn reshape(v: Self, dims: []const usize) ReshapeError!Self {
            _ = v;
            _ = dims;
            return error.NotImplemented;
        }

        /// The view over `shape`, dimensions of 1 expanding by stride 0.
        pub fn broadcastTo(v: Self, shape: Shape) (BroadcastError || error{NotImplemented})!Self {
            _ = try broadcastShape(v.shape, shape);
            return error.NotImplemented;
        }

        /// Rank 2 with a first stride of 1, as a matrix view.
        pub fn asMatrix(v: Self) AsMatrixError!(if (read_only) matrix.ConstView(T) else matrix.View(T)) {
            if (v.shape.rank != 2) return error.RankMismatch;
            return error.NotImplemented;
        }
    };
}

/// An owned array, contiguous and column-major.
pub fn Array(comptime T: type) type {
    return struct {
        data: []T,
        shape: Shape,

        const Self = @This();

        pub fn init(allocator: Allocator, dims: []const usize) (Allocator.Error || error{RankTooHigh})!Self {
            const shape = try Shape.init(dims);
            return .{ .data = try allocator.alloc(T, shape.count()), .shape = shape };
        }

        pub fn deinit(self: *Self, allocator: Allocator) void {
            allocator.free(self.data);
            self.* = undefined;
        }

        /// The whole array as a view: strides by the dimensions before.
        pub fn view(self: Self) ArrayView(T) {
            var strides: [max_rank]isize = @splat(0);
            var s: isize = 1;
            for (0..self.shape.rank) |k| {
                strides[k] = s;
                s *= @intCast(self.shape.dims[k]);
            }
            return .{ .data = self.data, .shape = self.shape, .strides = strides };
        }
    };
}

/// The shape both broadcast to: equal ranks, each pair of dimensions equal
/// or one of them 1 (never padded; design.md's Broadcasting).
pub fn broadcastShape(a: Shape, b: Shape) BroadcastError!Shape {
    if (a.rank != b.rank) return error.RankMismatch;
    var s: Shape = .{ .rank = a.rank };
    for (0..a.rank) |k| {
        const x = a.dims[k];
        const y = b.dims[k];
        s.dims[k] = if (x == y or y == 1) x else if (x == 1) y else return error.DimensionMismatch;
    }
    return s;
}

// ---- tests ----

const testing = std.testing;

/// A stub's tests skip until its body lands.
fn skip(err: anyerror) anyerror {
    return if (err == error.NotImplemented) error.SkipZigTest else err;
}

test "Array, Shape, broadcastShape" {
    var a = try Array(f64).init(testing.allocator, &.{ 2, 3, 4 });
    defer a.deinit(testing.allocator);
    try testing.expectEqual(24, a.data.len);
    const v = a.view();
    try testing.expectEqualSlices(isize, &.{ 1, 2, 6 }, v.strides[0..3]);
    const s = try broadcastShape(try .init(&.{ 2, 1, 4 }), try .init(&.{ 2, 3, 1 }));
    try testing.expectEqualSlices(usize, &.{ 2, 3, 4 }, s.dims[0..3]);
    try testing.expectError(error.RankMismatch, broadcastShape(try .init(&.{2}), try .init(&.{ 2, 1 })));
    try testing.expectError(error.DimensionMismatch, broadcastShape(try .init(&.{2}), try .init(&.{3})));
    try testing.expectEqual(48, @sizeOf(StaticArray(.{ 2, 3 }, f64)));
    _ = v.permute(&.{ 1, 0, 2 }) catch |e| return skip(e);
}
