//! The `symm` and `hemm` packing copies (`kernel.zig`), and tests of both
//! routines.

pub const kernel = @import("kernel.zig");
const check3 = @import("../../_test/blas3/check.zig");
const Vtable = @import("../../Vtable.zig");
const check = @import("../../_test/check.zig");
const portable = @import("../../portable.zig");

fn cases(vt: *const Vtable) !void {
    try check3.routineAll(vt, .symm);
    try check3.routineAll(vt, .hemm);
}

test "portable symm and hemm match the wide reference through the driver" {
    var vt: Vtable = .{};
    portable.install(&vt);
    try cases(&vt);
}

test "every linked core's symm and hemm match the wide reference, where this machine runs them" {
    try check.eachCore(cases);
}
