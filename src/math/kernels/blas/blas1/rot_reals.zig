//! `CSROT`/`ZDROT` through the real rot kernel (`interface/zrot.c`, OpenBLAS
//! 04511ccb1): real `c` and `s` rotate the real and imaginary parts alike,
//! so at unit strides `n` complex elements are `2n` reals. On x86_64 the
//! real kernels of Haswell and newer are about twice the complex SSE ones
//! on long vectors.

const sysinfo = @import("builtin");
const blas = @import("../__root__.zig");
const types = @import("../types.zig");
const util = @import("../util.zig");
const portable = @import("../portable.zig");

const Float = types.Float;

/// `zrot.c`'s `RROT_MIN`: below, the real kernels cost more per call.
pub const n_min = 4096;
/// `zrot.c`'s `RROT_MAX`: keeps `2n` in the reference's `int`.
pub const n_max = 1 << 29;

/// Whether `T` routes at all: `zrot.c`'s `ARCH_X86_64` and not `XDOUBLE`.
pub fn serves(comptime T: type) bool {
    const F = Float(T);
    return sysinfo.cpu.arch == .x86_64 and types.isComplex(T) and (F == f32 or F == f64);
}

/// `zrot.c`'s `RROT_USE(n) && incx == 1 && incy == 1`. The real kernels
/// fuse where the complex ones do not, so the bits differ: the caller asks
/// once, for the whole `n`, never per chunk.
pub fn routed(comptime T: type, n: isize, incx: isize, incy: isize) bool {
    if (comptime !serves(T)) return false;
    return n >= n_min and n <= n_max and incx == 1 and incy == 1;
}

/// `len` unit-stride complex elements of `x` and `y`, rotated as `2 * len`
/// reals; `T` one `serves` accepts.
pub fn run(comptime T: type, len: isize, x: [*]T, y: [*]T, c: Float(T), s: Float(T)) void {
    _ = portable.fetch(Float(T), .rot_k, blas.vtable)(2 * len, util.raw(T, x), 1, util.raw(T, y), 1, c, s);
}
