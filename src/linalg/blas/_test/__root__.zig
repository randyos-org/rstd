//! Test-only code: each level's shared test helpers, the edge-case and
//! denormal tests (`blas1/edge.zig`), and OpenBLAS's own suites ported (`from_openblas/`).

test {
    _ = @import("types.zig");
    _ = @import("check.zig");
    _ = @import("bitwise.zig");
    _ = @import("blas1/edge.zig");
    _ = @import("blas1/extremum.zig");
    _ = @import("blas1/threaded.zig");
    _ = @import("blas2/edge.zig");
    _ = @import("blas2/threaded.zig");
    _ = @import("bf16/gemv.zig");
    _ = @import("blas3/gemmt.zig");
    _ = @import("blas3/threaded.zig");
    _ = @import("from_openblas/__root__.zig");
}
