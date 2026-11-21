//! Tests of the `trmm` kernels.

const check3 = @import("../../_test/blas3/check.zig");
const Vtable = @import("../../Vtable.zig");
const check = @import("../../_test/check.zig");
const portable = @import("../../portable.zig");

fn cases(vt: *const Vtable) !void {
    try check3.routineAll(vt, .trmm);
}

test "portable trmm matches the wide reference through the driver" {
    var vt: Vtable = .{};
    portable.install(&vt);
    try cases(&vt);
}

test "every linked core's trmm matches the wide reference, where this machine runs it" {
    try check.eachCore(cases);
}
