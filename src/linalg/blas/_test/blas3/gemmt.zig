//! The tests of `gemmt` through the driver (`driver.gemmt`), which runs the
//! gemm kernels over one triangle of `C`; it has no kernels of its own.

const check3 = @import("check.zig");
const check = @import("../check.zig");
const portable = @import("../../portable.zig");
const Vtable = @import("../../Vtable.zig");

fn cases(vt: *const Vtable) !void {
    try check3.routineAll(vt, .gemmt);
}

test "portable gemmt matches the wide reference through the driver" {
    var vt: Vtable = .{};
    portable.install(&vt);
    try cases(&vt);
}

test "every linked core's gemmt matches the wide reference, where this machine runs it" {
    try check.eachCore(cases);
}
