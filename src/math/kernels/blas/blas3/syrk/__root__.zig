//! Tests of the `syrk` and `herk` kernels.

const check3 = @import("../../_test/blas3/check.zig");
const Vtable = @import("../../Vtable.zig");
const check = @import("../../_test/check.zig");
const portable = @import("../../portable.zig");

fn cases(vt: *const Vtable) !void {
    try check3.routineAll(vt, .syrk);
    try check3.routineAll(vt, .herk);
}

test "portable syrk and herk match the wide reference through the driver" {
    var vt: Vtable = .{};
    portable.install(&vt);
    try cases(&vt);
}

test "every linked core's syrk and herk match the wide reference, where this machine runs them" {
    try check.eachCore(cases);
}
