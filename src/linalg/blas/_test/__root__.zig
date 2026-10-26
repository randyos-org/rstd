//! Test-only code: each level's shared test helpers, the edge-case and
//! denormal tests, and OpenBLAS's own suites ported (`from_openblas/`).

test {
    _ = @import("types.zig");
    _ = @import("blas1/check.zig");
    _ = @import("blas1/edge.zig");
    _ = @import("blas1/denormal.zig");
    _ = @import("blas2/edge.zig");
    _ = @import("from_openblas/__root__.zig");
}
