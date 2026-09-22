//! AArch64 `daxpy` block for the ThunderX: the function `daxpy_kernel_8` of
//! OpenBLAS `kernel/arm64/daxpy_thunderx.c`, which the shared driver in
//! `../kernel.zig` calls for the contiguous part in blocks of 32.
//!
//! The reference is C that is scheduled by hand around prefetches and empty
//! `asm` statements that pin the results in vector registers so the compiler
//! keeps the eight multiply-adds together. The prefetches (`prfm PLDL1STRM`)
//! are `@prefetch`es here, three lines ahead of the block as they are there;
//! the register pinning has no counterpart and the loop is otherwise the same.

/// `y[0..n] += alpha * x[0..n]`, eight elements at a time.
pub fn daxpy_kernel_8(n: isize, x_in: [*]f64, y_in: [*]f64, alpha: [*]const f64) void {
    const a = alpha[0];
    var x = x_in;
    var y = y_in;
    const line = 128 / @sizeOf(f64);

    // the distances the reference starts with
    @prefetch(x + 2 * line, .{ .locality = 0 });
    @prefetch(y + 2 * line, .{ .locality = 0 });
    @prefetch(x + 3 * line, .{ .locality = 0 });
    @prefetch(y + 3 * line, .{ .locality = 0 });
    @prefetch(x + 4 * line, .{ .locality = 0 });
    @prefetch(y + 4 * line, .{ .locality = 0 });

    var i: isize = 0;
    while (i < n) : (i += 8) {
        var y0: [8]f64 = undefined;
        inline for (0..8) |k| y0[k] = @mulAdd(f64, a, x[k], y[k]);
        inline for (0..8) |k| y[k] = y0[k];

        @prefetch(x + 4 * line, .{ .locality = 0 });
        @prefetch(y + 4 * line, .{ .locality = 0 });

        y += 8;
        x += 8;
    }
}
