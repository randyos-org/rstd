//! The `rstd_matrix` module where the compiler has `@Matrix` (this project's
//! Zig fork). Stock Zig rejects the builtin's name in any file it parses, so
//! rstd's build picks this file or `fallback.zig`; only code that spells
//! `@Matrix` or `@transpose` belongs here, the rest in rstd, branching on
//! `has_matrix`.

pub const has_matrix = true;

pub fn Matrix(comptime rows: comptime_int, comptime cols: comptime_int, comptime T: type) type {
    return @Matrix(rows, cols, T);
}

pub fn Transposed(comptime M: type) type {
    const info = @typeInfo(M).matrix;
    return @Matrix(info.cols, info.rows, info.child);
}

pub inline fn transpose(m: anytype) Transposed(@TypeOf(m)) {
    return @transpose(m);
}

test transpose {
    const m: Matrix(2, 3, f32) = .{ .{ 1, 2 }, .{ 3, 4 }, .{ 5, 6 } };
    const t: [2][3]f32 = transpose(m);
    try @import("std").testing.expectEqual([2][3]f32{ .{ 1, 3, 5 }, .{ 2, 4, 6 } }, t);
}
