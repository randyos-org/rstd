//! Root of the WASI tests (`build/wasmtest.zig`): the rigs' OS-free tests
//! (`_test_qemu.zig`) and LAPACK, whose allocations WASI serves.

test {
    _ = @import("_test_qemu.zig");
    _ = @import("linalg/lapack/__root__.zig");
}
