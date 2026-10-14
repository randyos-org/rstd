//! Root of the tests that run on emulated machines (see `build/qemutest.zig`):
//! hardware-independent code only, since there is no OS underneath (so nothing
//! that needs an allocator, files or threads). Extend as more of rstd qualifies.

test {
    _ = @import("simd/__root__.zig");
    // the whole BLAS tree; `-Dqemu-test-filter` picks what runs
    _ = @import("linalg/blas/__root__.zig");
}
