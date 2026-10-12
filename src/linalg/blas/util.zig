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
