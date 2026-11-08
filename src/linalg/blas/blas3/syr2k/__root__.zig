//! Tests of the `syr2k` and `her2k` kernels.

const check3 = @import("../../_test/blas3/check.zig");
const Vtable = @import("../../Vtable.zig");
const check = @import("../../_test/check.zig");
const portable = @import("../../portable.zig");

fn cases(vt: *const Vtable) !void {
    try check3.routineAll(vt, .syr2k);
    try check3.routineAll(vt, .her2k);
}

test "portable syr2k and her2k match the wide reference through the driver" {
    var vt: Vtable = .{};
    portable.install(&vt);
    try cases(&vt);
}

test "every linked core's syr2k and her2k match the wide reference, where this machine runs them" {
    try check.eachCore(cases);
}
