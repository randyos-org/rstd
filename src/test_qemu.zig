//! Root of the tests that run on emulated machines (see `build/qemutest.zig`):
//! hardware-independent code only, since there is no OS underneath (so nothing
//! that needs an allocator, files or threads). Extend as more of rstd qualifies.

test {
    _ = @import("simd/__root__.zig");
    _ = @import("linalg/blas/blas1/axpy/__root__.zig");
}
