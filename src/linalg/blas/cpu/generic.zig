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
const blas1 = @import("../blas1/__root__.zig");
const blas = @import("../__root__.zig");
const gemm_params = @import("../blas3/gemm/params.zig");
pub fn vtableFromCoreType(core_type: CoreType) Vtable {
    var vt: Vtable = .{};
    blas1.installGeneric(&vt);
    kernel.install(&vt, core_type);
    gemm_params.installCore(&vt, .generic, core_type, .fromConfig(blas.config));
    return vt;
}

test "generic detect is always available" {
    try zigstd.testExpectEqual(CoreType.generic, detect());
}
