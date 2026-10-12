//! 32-bit ARM CPU detection (`_detect/arm.zig`), with its kernel registry
//! and vtable. The detector is its own file so the build can use it
//! without importing the kernels.

const detector = @import("_detect/arm.zig");
const Vtable = @import("../Vtable.zig");
const blas1 = @import("../blas1/__root__.zig");
const blas = @import("../__root__.zig");
const gemm_params = @import("../blas3/gemm/params.zig");

pub const CoreType = detector.CoreType;
pub const kernel = @import("../kernels/arm.zig");
pub const detect = detector.detect;

pub fn vtableFromCoreType(core_type: CoreType) Vtable {
    var vt: Vtable = .{};
    blas1.installGeneric(&vt);
    kernel.install(&vt, core_type);
    gemm_params.installCore(&vt, .arm, core_type, .fromConfig(blas.config));
    return vt;
}

test {
    _ = detector;
}
