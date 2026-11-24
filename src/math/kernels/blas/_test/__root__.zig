//! Test-only code: each level's shared test helpers, the edge-case and
//! denormal tests (`blas1/edge.zig`), and OpenBLAS's own suites ported (`from_openblas/`).

test {
    _ = @import("types.zig");
    _ = @import("check.zig");
    _ = @import("bitwise.zig");
    _ = @import("blas1/edge.zig");
    _ = @import("blas1/extremum.zig");
    _ = @import("blas1/interface.zig");
    _ = @import("blas1/interface_any.zig");
    _ = @import("blas1/rot_reals.zig");
    _ = @import("blas1/threaded.zig");
    _ = @import("blas2/banded.zig");
    _ = @import("blas2/edge.zig");
    _ = @import("blas2/gemv.zig");
    _ = @import("blas2/interface.zig");
    _ = @import("blas2/interface_any.zig");
    _ = @import("blas2/threaded.zig");
    _ = @import("blas2/triangular.zig");
    _ = @import("bf16/gemv.zig");
    _ = @import("blas3/driver.zig");
    _ = @import("blas3/gemmt.zig");
    _ = @import("blas3/half.zig");
    _ = @import("blas3/interface.zig");
    _ = @import("blas3/params.zig");
    _ = @import("blas3/small.zig");
    _ = @import("blas3/threaded.zig");
    _ = @import("from_openblas/__root__.zig");
}
