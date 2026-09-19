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
