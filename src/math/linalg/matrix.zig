//! `rstd.math.linalg`'s matrices (design.md, "Planned: the public `rstd.math.linalg`
//! layer"): `StaticMat(R, C, T)` and free functions over it, `T` a float or
//! its `std.math.Complex` (interleaved), inferring the sizes and `T`; vectors
//! are `linalg.vector`'s. Shaderlib re-exports them under GLSL names.

const std = @import("std");
const ops = @import("../kernels/lapack/ops.zig");
const root = @import("root");
const zigstd = @import("../../zigstd.zig");
const mops = @import("matrix/ops.zig");
const msol = @import("matrix/solve.zig");
const mplumb = @import("matrix/plumbing.zig");

pub const add = mops.add;
pub const col = mops.col;
pub const ctranspose = mops.ctranspose;
pub const fromCols = mops.fromCols;
pub const fromRows = mops.fromRows;
pub const identity = mops.identity;
pub const mul = mops.mul;
const mulAnySize = mops.mulAnySize;
pub const outer = mops.outer;
pub const row = mops.row;
pub const scale = mops.scale;
pub const sub = mops.sub;
pub const trace = mops.trace;
pub const transpose = mops.transpose;
pub const transposeVecs = mops.transposeVecs;
pub const det = msol.det;
const detAnySize = msol.detAnySize;
pub const inv = msol.inv;
const invAnySize = msol.invAnySize;
pub const luInPlace = msol.luInPlace;
pub const luSolve = msol.luSolve;
pub const solve = msol.solve;
const solveAnySize = msol.solveAnySize;
pub const eMul = mplumb.eMul;
pub const eMulAdd = mplumb.eMulAdd;

/// `rstd.math.linalg`'s build settings, from a root's `rstd_linalg_options`
/// (design.md, Dispatch, decision 3).
pub const Options = struct {
    /// `L`, per dimension: at or below it the inline kernels, past it the
    /// BLAS and LAPACK. At least 4, for shaderlib's `Mat2` to `Mat4`.
    inline_limit: usize = 12,
};

/// The root's options, or the defaults.
pub const options: Options = if (@hasDecl(root, "rstd_linalg_options")) root.rstd_linalg_options else .{};

pub const inline_limit = options.inline_limit;

comptime {
    if (inline_limit < 4) @compileError(std.fmt.comptimePrint("rstd_linalg_options.inline_limit is {d}, below 4 (shaderlib's Mat4)", .{inline_limit}));
}

/// `mul`, `det`, `inv` and `solve` at every size, past the inline limit
/// too: `rstd.math.linalg`'s inline path (its tests run it to 16) and timing.
pub const any_size = struct {
    pub const mul = mulAnySize;
    pub const det = detAnySize;
    pub const inv = invAnySize;
    pub const solve = solveAnySize;
};

/// `R` rows by `C` columns of `T`, packed column-major: `data[j]` is column
/// `j` (design.md, decisions 1 and 14). `mul`, `transpose`, `scale`, `add`
/// and `sub` of real `T` go through `@Matrix` where the compiler has it
/// (`rstd_matrix`), else through columns in its order, to the same bits
/// (decision 16).
pub fn StaticMat(comptime R: comptime_int, comptime C: comptime_int, comptime T: type) type {
    return struct {
        data: [C][R]T,

        pub const rows = R;
        pub const cols = C;
        pub const Element = T;

        /// A view of `m`'s storage, valid while `m` is.
        pub fn view(m: *@This()) View(T) {
            const flat: *[C * R]T = @ptrCast(&m.data);
            return .{ .data = flat, .rows = R, .cols = C, .ld = R };
        }

        pub fn constView(m: *const @This()) ConstView(T) {
            const flat: *const [C * R]T = @ptrCast(&m.data);
            return .{ .data = flat, .rows = R, .cols = C, .ld = R };
        }
    };
}

/// A borrowed `rows` by `cols` matrix in `data`, column-major with column
/// stride `ld` (`data[i + j * ld]`): all or part of a `Mat`, a `StaticMat`
/// or the caller's storage. `ConstView` is its read-only form.
pub fn View(comptime T: type) type {
    return ViewOf(T, false);
}

pub fn ConstView(comptime T: type) type {
    return ViewOf(T, true);
}

fn ViewOf(comptime T: type, comptime read_only: bool) type {
    return struct {
        data: if (read_only) []const T else []T,
        rows: usize,
        cols: usize,
        ld: usize,

        const Self = @This();
        pub const Element = T;

        /// `error.DimensionMismatch` for `ld < rows`, `error.SliceTooShort`
        /// when `data` ends before the last element.
        pub fn init(data: @FieldType(Self, "data"), rows: usize, cols: usize, ld: usize) error{ DimensionMismatch, SliceTooShort }!Self {
            if (ld < @max(rows, 1)) return error.DimensionMismatch;
            if (rows > 0 and cols > 0 and data.len < (cols - 1) * ld + rows) return error.SliceTooShort;
            return .{ .data = data, .rows = rows, .cols = cols, .ld = ld };
        }

        pub fn at(v: Self, i: usize, j: usize) T {
            return v.data[i + j * v.ld];
        }

        /// Rows `i..i + rows`, columns `j..j + cols`, sharing `v`'s storage.
        pub fn sub(v: Self, i: usize, j: usize, rows: usize, cols: usize) error{DimensionMismatch}!Self {
            if (i + rows > v.rows or j + cols > v.cols) return error.DimensionMismatch;
            const start = if (rows > 0 and cols > 0) i + j * v.ld else 0;
            const end = if (rows > 0 and cols > 0) start + (cols - 1) * v.ld + rows else 0;
            return .{ .data = v.data[start..end], .rows = rows, .cols = cols, .ld = v.ld };
        }

        /// Column `j`, contiguous; empty, and in bounds, with no rows.
        pub fn col(v: Self, j: usize) @FieldType(Self, "data") {
            if (v.rows == 0) return v.data[0..0];
            return v.data[j * v.ld ..][0..v.rows];
        }

        pub fn constView(v: Self) ConstView(T) {
            return .{ .data = v.data, .rows = v.rows, .cols = v.cols, .ld = v.ld };
        }
    };
}

/// A `rows` by `cols` matrix on the heap, packed column-major (`ld = rows`),
/// unmanaged as `std.ArrayList`: freed with `deinit` on the allocator
/// that made it.
pub fn Mat(comptime T: type) type {
    return struct {
        data: []T,
        rows: usize,
        cols: usize,

        const Self = @This();
        pub const Element = T;

        /// Elements undefined.
        pub fn init(allocator: std.mem.Allocator, rows: usize, cols: usize) std.mem.Allocator.Error!Self {
            return .{ .data = try allocator.alloc(T, rows * cols), .rows = rows, .cols = cols };
        }

        pub fn zeros(allocator: std.mem.Allocator, rows: usize, cols: usize) std.mem.Allocator.Error!Self {
            const m = try init(allocator, rows, cols);
            @memset(m.data, ops.zero(T));
            return m;
        }

        pub fn deinit(m: Self, allocator: std.mem.Allocator) void {
            allocator.free(m.data);
        }

        pub fn at(m: Self, i: usize, j: usize) T {
            return m.data[i + j * m.rows];
        }

        pub fn view(m: Self) View(T) {
            return .{ .data = m.data, .rows = m.rows, .cols = m.cols, .ld = @max(m.rows, 1) };
        }

        pub fn constView(m: Self) ConstView(T) {
            return .{ .data = m.data, .rows = m.rows, .cols = m.cols, .ld = @max(m.rows, 1) };
        }
    };
}

pub const M2f32 = StaticMat(2, 2, f32);
pub const M3f32 = StaticMat(3, 3, f32);
pub const M4f32 = StaticMat(4, 4, f32);
pub const M2f64 = StaticMat(2, 2, f64);
pub const M3f64 = StaticMat(3, 3, f64);
pub const M4f64 = StaticMat(4, 4, f64);

/// `mul` into `out`.
pub fn mulInto(out: anytype, a: anytype, b: anytype) void {
    out.* = mul(a, b);
}

/// `transpose` into `out`.
pub fn transposeInto(out: anytype, m: anytype) void {
    out.* = transpose(m);
}

/// `ctranspose` into `out`.
pub fn ctransposeInto(out: anytype, m: anytype) void {
    out.* = ctranspose(m);
}

/// `outer` into `out`.
pub fn outerInto(out: anytype, a: anytype, b: anytype) void {
    out.* = outer(a, b);
}

/// `scale` into `out`.
pub fn scaleInto(out: anytype, m: anytype, s: @TypeOf(m).Element) void {
    out.* = scale(m, s);
}

/// `add` into `out`.
pub fn addInto(out: anytype, a: anytype, b: @TypeOf(a)) void {
    out.* = add(a, b);
}

/// `sub` into `out`.
pub fn subInto(out: anytype, a: anytype, b: @TypeOf(a)) void {
    out.* = sub(a, b);
}

/// `inv` into `out`, untouched on `error.Singular`.
pub fn invInto(out: anytype, m: anytype) error{Singular}!void {
    out.* = try inv(m);
}

/// `solve` into `out`, untouched on `error.Singular`.
pub fn solveInto(out: anytype, a: anytype, b: anytype) error{Singular}!void {
    out.* = try solve(a, b);
}

test {
    zigstd.testRefAllDecls(@This());
    _ = @import("_test/matrix.zig");
}
