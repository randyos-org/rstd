//! Root of the tests that run on emulated machines (see `build/qemutest.zig`):
//! hardware-independent code only, since there is no OS underneath (so nothing
//! that needs an allocator, files or threads). Extend as more of rstd qualifies.
//! The one rig file in `src/`: a module's tests are those its own files reach,
//! and it cannot import past its directory. The rest is under `test/`.

test {
    _ = @import("simd/__root__.zig");
    // the whole BLAS tree; `-Dqemu-test-filter` picks what runs
    _ = @import("linalg/blas/__root__.zig");
}
