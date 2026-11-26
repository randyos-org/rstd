//! A routed call's parts and arguments (`apple_accel.zig`).
//!
//! Accelerate runs single-threaded; a call splits by size alone into parts on
//! rstd's workers, one Accelerate call each (threading.md, Vendor backends).

const std = @import("std");
const types = @import("../types.zig");
const threading = @import("../threading.zig");
const binding = @import("binding.zig");
const routines = @import("routines.zig");

const Io = std.Io;
const isComplex = types.isComplex;
const Float = types.Float;
const Blas = binding.Blas;
const Symbol = routines.Symbol;
const Order = binding.Order;
const Transpose = binding.Transpose;
const Uplo = binding.Uplo;
const Diag = binding.Diag;
const Side = binding.Side;

/// Level 1's parts start on multiples of this many elements, as rstd's own
/// do (`blas1/interface.zig`).
pub const vector_quantum = 1024;
/// Level 2's and 3's on multiples of this many rows or columns, as rstd's
/// level 2 (`blas2/interface.zig`): parts stay off each other's cache lines.
pub const matrix_quantum = 64;

/// `part(ctx, chunk)` for each of `family`'s chunks of `len` units, `work`
/// being its rule's measure: inline for one, else on `io`'s workers.
pub fn split(io: Io, comptime family: threading.Family, comptime T: type, work: u64, len: usize, quantum: usize, ctx: anytype, comptime part: fn (@TypeOf(ctx), threading.Chunk) void) Io.Cancelable!void {
    const p = threading.plan(family, T, work, len, .even, quantum);
    if (p.chunks.count == 1) return part(ctx, p.chunks.get(0));
    return threading.forEach(io, &p, ctx, part);
}

/// `split` for a reduction over `len` elements: the parts' sums, added in
/// chunk order.
pub fn fold(io: Io, comptime family: threading.Family, comptime T: type, comptime W: type, len: usize, ctx: anytype, comptime part: fn (@TypeOf(ctx), threading.Chunk) W) Io.Cancelable!W {
    const p = threading.plan(family, T, len, len, .even, vector_quantum);
    if (p.chunks.count == 1) return part(ctx, p.chunks.get(0));
    const S = struct {
        fn add(a: W, b: W) W {
            return if (comptime isComplex(W)) a.add(b) else a + b;
        }
    };
    return threading.reduce(io, &p, W, ctx, part, S.add);
}

/// What a reduction over `T` sums its parts in: `f64` for `f32`, so the
/// result rounds once, as rstd's own fold (`blas1/interface.zig`).
pub fn Sum(comptime T: type) type {
    const W = if (Float(T) == f32) f64 else Float(T);
    return if (isComplex(T)) std.math.Complex(W) else W;
}

pub fn rounded(comptime T: type, w: Sum(T)) T {
    if (comptime isComplex(T)) return .init(@floatCast(w.re), @floatCast(w.im));
    return @floatCast(w);
}

/// Whether the table has `sym` and, on LP64, every entry of `ints` (each
/// integer argument's magnitude and each operand's span) fits `c_int`; then
/// so does any call over a part of the operands.
pub fn fits(comptime sym: Symbol, b: *const Blas, ints: []const isize) bool {
    switch (b.*) {
        inline else => |*bound, tag| {
            if (@field(bound.fns, @tagName(sym)) == null) return false;
            if (tag == .lp64) for (ints) |e| if (e > std.math.maxInt(c_int)) return false;
            return true;
        },
    }
}

/// `|v|`, as `fits` takes a stride.
pub fn mag(v: isize) isize {
    return @intCast(@abs(v));
}

/// The whole of `len` units, as one chunk.
pub fn whole(len: isize) threading.Chunk {
    return .{ .index = 0, .start = 0, .end = @intCast(len) };
}

/// Elements `c.start..c.end` of the `n` at `p`, stride `inc`, as CBLAS takes
/// a vector: from its lowest address, which is the last one for `inc < 0`.
pub fn sub(p: anytype, n: isize, inc: isize, c: threading.Chunk) @TypeOf(p) {
    const at = if (inc >= 0) c.start else @as(usize, @intCast(n)) - c.end;
    return p + at * @abs(inc);
}

/// Column `j` of a matrix at `p`, leading dimension `ld`.
pub fn column(p: anytype, j: usize, ld: isize) @TypeOf(p) {
    return p + j * @as(usize, @intCast(ld));
}

pub fn units(c: threading.Chunk) isize {
    return @intCast(c.end - c.start);
}

/// A scalar as CBLAS takes it: by value if real, by pointer if complex.
pub inline fn scalar(comptime T: type, v: *const T) if (isComplex(T)) *const anyopaque else T {
    return if (comptime isComplex(T)) v else v.*;
}

/// CBLAS's `TRANS` for `t`; null for the letters it has none for (complex
/// `r`, and `o` to `d`). A real `r` and `c` are `n` and `t`.
pub fn transpose(comptime T: type, t: types.Trans) ?Transpose {
    if (t.conjX()) return null;
    if (comptime !isComplex(T)) return if (t.transposed()) .trans else .no_trans;
    return switch (t) {
        .n => .no_trans,
        .t => .trans,
        .c => .conj_trans,
        else => null,
    };
}

/// CBLAS's `UPLO`; null for `HEMVREV`'s `v` and `m`.
pub fn uplo(u: types.Uplo) ?Uplo {
    return switch (u) {
        .u => .upper,
        .l => .lower,
        .v, .m => null,
    };
}

pub fn diag(d: types.Diag) Diag {
    return if (d == .u) .unit else .non_unit;
}

pub fn side(s: types.Side) Side {
    return if (s == .l) .left else .right;
}

pub const col = Order.col_major;
