//! Root of the tests that run on emulated machines (see `build/qemutest.zig`):
//! hardware-independent code only, since there is no OS underneath (so nothing
//! that needs files or threads; the runner's heap serves allocations). Extend
//! as more of rstd qualifies.
//! With `_test_wasm.zig`, the rig files in `src/`: a module's tests are those
//! its own files reach, and it cannot import past its directory. The rest is
//! under `test/`.

const sysinfo = @import("builtin");

/// The UEFI rigs' `std.Io` (`test/qemu/uefi.zig`), which this module imports
/// as `qemutest_root`; the runner's files cannot reach `src/`.
pub const uefi_io = if (sysinfo.target.os.tag == .uefi) @import("_os/uefi/io/__root__.zig") else struct {};

test {
    // the whole BLAS tree; `-Dqemu-test-filter` picks what runs
    _ = @import("linalg/blas/__root__.zig");
    _ = @import("linalg/lapack/__root__.zig");
}
