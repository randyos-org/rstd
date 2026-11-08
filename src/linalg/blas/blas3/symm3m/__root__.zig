//! Tests of the `symm3m` kernels.

// ---- tests -----------------------------------------------------------------

const check3 = @import("../../_test/blas3/check.zig");
const check = @import("../../_test/check.zig");
const Vtable = @import("../../Vtable.zig");
const portable = @import("../../portable.zig");

fn cases(vt: *const Vtable) !void {
    try check3.routineAll(vt, .symm3m);
}

test "portable symm3m kernels match the wide reference through the driver" {
    var vt: Vtable = .{};
    portable.install(&vt);
    try cases(&vt);
}

test "every linked core's symm3m kernels match the wide reference, where this machine runs them" {
    try check.eachCore(cases);
}
