//! Fallback for architectures with no detector of their own.
//!
//! The reference has the same notion: a `<ARCH>_GENERIC` kernel set built
//! from portable C, used whenever nothing more specific applies.

const std = @import("std");
const coretypes = @import("_coretypes.zig");
const zigstd = @import("../../../zigstd.zig");

pub const CoreType = coretypes.GenericCoreType;
pub const kernel = @import("../kernels/generic.zig");

pub fn detect() CoreType {
    return .generic;
}

const Vtable = @import("../Vtable.zig");
const axpy = @import("../blas1/axpy/__root__.zig");
pub fn vtableFromCoreType(core_type: CoreType) Vtable {
    var vt: Vtable = .{};
    axpy.install(&vt, core_type);
    return vt;
}

test "generic detect is always available" {
    try zigstd.testExpectEqual(CoreType.generic, detect());
}
