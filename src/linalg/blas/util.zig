const std = @import("std");

/// Equivalent of C's `base[offset]` where `offset` may be negative -- BLAS
/// routines allow negative `inc_x`/`inc_y` (iterate a vector in reverse),
/// which does not map onto Zig's `usize`-indexed pointer/slice arithmetic
/// the way it does onto a C pointer (`x[ix]` with negative `ix` is just
/// `*(x + ix)` in C).
pub inline fn at(comptime T: type, base: [*]T, offset: i64) *T {
    const addr: usize = @intFromPtr(base);
    const byte_offset: i64 = offset * @as(i64, @sizeOf(T));
    const new_addr: usize = if (byte_offset >= 0)
        addr + @as(usize, @intCast(byte_offset))
    else
        addr - @as(usize, @intCast(-byte_offset));
    return @ptrFromInt(new_addr);
}

/// The symbol of the kernel `role` (`daxpy`, ...) of `core`. The per-core objects
/// (`build/blaskern.zig`) export it, with hidden visibility so it never leaves a
/// shared library rstd is built into, and `kernels/<arch>.zig`'s `install` looks
/// it up by the same name.
pub fn kernelSymbol(comptime role: []const u8, comptime core: []const u8) []const u8 {
    return "blaskern_" ++ role ++ "_" ++ core;
}

/// Loads a `V` from `p`, which is not assumed aligned to `V`'s natural
/// alignment -- e.g. a kernel walking a caller's buffer by an arbitrary
/// stride, where nothing guarantees the vector-width alignment `@Vector`
/// otherwise wants. (Lives here, not in `simd/`, because the per-core kernel
/// objects that use it are their own module rooted at `_kernel.zig` and
/// cannot reach outside `linalg/blas/`; see `kernels/__root__.zig`.)
pub inline fn load(comptime V: type, p: [*]const @typeInfo(V).vector.child) V {
    return @as(*align(1) const V, @ptrCast(p)).*;
}

/// Stores `v` to `p`, unaligned; see `load`.
pub inline fn store(comptime V: type, p: [*]@typeInfo(V).vector.child, v: V) void {
    @as(*align(1) V, @ptrCast(p)).* = v;
}

/// Whether a vector of `n` elements of `size` bytes at stride `inc` fits in
/// the address space: `((n - 1) * |inc| + 1) * size` bytes, with no overflow.
/// Where it does not, no such vector can exist; the interfaces refuse it so
/// no index computation overflows.
pub fn fits(n: isize, inc: isize, size: usize) bool {
    if (inc == std.math.minInt(isize)) return false;
    if (n <= 1) return true;
    return fitsMatrix(1, n, @intCast(@abs(inc)), size);
}

/// Whether a column-major `rows` by `cols` matrix with leading dimension `ld`
/// fits, as `fits` does for a vector: `(cols - 1) * ld + rows` elements.
pub fn fitsMatrix(rows: isize, cols: isize, ld: isize, size: usize) bool {
    if (rows <= 0 or cols <= 0) return true;
    const span = std.math.mul(isize, cols - 1, ld) catch return false;
    const elems = std.math.add(isize, span, rows) catch return false;
    _ = std.math.mul(isize, elems, std.math.cast(isize, size) orelse return false) catch return false;
    return true;
}
