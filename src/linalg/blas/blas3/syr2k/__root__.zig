//! Tests of the `syr2k` and `her2k` kernels.

const check = @import("../../_test/blas3/check.zig");
const Vtable = @import("../../Vtable.zig");

fn cases(vt: *const Vtable) !void {
    try check.routineAll(vt, .syr2k);
    try check.routineAll(vt, .her2k);
}

const blas1check = @import("../../_test/blas1/check.zig");
const portable = @import("../../portable.zig");

test "portable syr2k and her2k match the wide reference through the driver" {
    var vt: Vtable = .{};
    portable.install(&vt);
    try cases(&vt);
}

test "every linked core's syr2k and her2k match the wide reference, where this machine runs them" {
    try blas1check.eachCore(cases);
}
