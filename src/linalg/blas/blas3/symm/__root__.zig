//! The `symm` and `hemm` packing copies (`kernel.zig`), and tests of both
//! routines.

pub const kernel = @import("kernel.zig");

const check = @import("../../_test/blas3/check.zig");
const Vtable = @import("../../Vtable.zig");

fn cases(vt: *const Vtable) !void {
    try check.routineAll(vt, .symm);
    try check.routineAll(vt, .hemm);
}

const blas1check = @import("../../_test/blas1/check.zig");
const portable = @import("../../portable.zig");

test "portable symm and hemm match the wide reference through the driver" {
    var vt: Vtable = .{};
    portable.install(&vt);
    try cases(&vt);
}

test "every linked core's symm and hemm match the wide reference, where this machine runs them" {
    try blas1check.eachCore(cases);
}
