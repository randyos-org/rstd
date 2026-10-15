//! AArch64 CPU detection (`_detect/arm64.zig`), with its kernel registry
//! and vtable. The detector is its own file so the build can use it
//! without importing the kernels.

const detector = @import("_detect/arm64.zig");
const Vtable = @import("../Vtable.zig");
const portable = @import("../portable.zig");
const blas = @import("../__root__.zig");
const gemm_params = @import("../blas3/gemm/params.zig");

pub const CoreType = detector.CoreType;
pub const kernel = @import("../kernels/arm64.zig");
pub const detect = detector.detect;

pub fn vtableFromCoreType(core_type: CoreType) Vtable {
    var vt: Vtable = .{};
    portable.install(&vt);
    kernel.install(&vt, core_type);
    portable.alias(&vt);
    gemm_params.installCore(&vt, .arm64, core_type, .fromConfig(blas.config));
    return vt;
}

test {
    _ = detector;
}
