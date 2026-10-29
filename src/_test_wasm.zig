//! Root of the WASI tests (`build/wasmtest.zig`): the rigs' OS-free tests
//! (`_test_qemu.zig`).

test {
    _ = @import("_test_qemu.zig");
}
