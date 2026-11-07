//! Tests of the `syrk` and `herk` kernels.

const check = @import("../../_test/blas3/check.zig");
const Vtable = @import("../../Vtable.zig");

fn cases(vt: *const Vtable) !void {
    try check.routineAll(vt, .syrk);
    try check.routineAll(vt, .herk);
}

const blas1check = @import("../../_test/blas1/check.zig");
const portable = @import("../../portable.zig");

test "portable syrk and herk match the wide reference through the driver" {
    var vt: Vtable = .{};
    portable.install(&vt);
    try cases(&vt);
}

test "every linked core's syrk and herk match the wide reference, where this machine runs them" {
    try blas1check.eachCore(cases);
}
