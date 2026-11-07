//! Tests of the `gemm3m` kernels.

// ---- tests -----------------------------------------------------------------

const check = @import("../../_test/blas3/check.zig");
const blas1check = @import("../../_test/blas1/check.zig");
const Vtable = @import("../../Vtable.zig");
const portable = @import("../../portable.zig");

fn cases(vt: *const Vtable) !void {
    try check.routineAll(vt, .gemm3m);
}

test "portable gemm3m kernels match the wide reference through the driver" {
    var vt: Vtable = .{};
    portable.install(&vt);
    try cases(&vt);
}

test "every linked core's gemm3m kernels match the wide reference, where this machine runs them" {
    try blas1check.eachCore(cases);
}

test {
    _ = @import("kernel.zig");
    _ = @import("driver.zig");
}
