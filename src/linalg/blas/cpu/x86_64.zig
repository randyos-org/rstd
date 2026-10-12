//! x86_64 CPU detection: the detector shared with 32-bit x86 (`_detect/x86.zig`),
//! with x86_64's own kernel registry and vtable.

const detector = @import("_detect/x86.zig");
const Vtable = @import("../Vtable.zig");
const blas1 = @import("../blas1/__root__.zig");
const blas = @import("../__root__.zig");
const gemm_params = @import("../blas3/gemm/params.zig");

pub const CoreType = detector.CoreType;
pub const kernel = @import("../kernels/x86_64.zig");
pub const detect = detector.detect;
/// Whether this machine can run `core`'s kernels (for tests).
pub const canRun = detector.supportsCore;

pub fn vtableFromCoreType(core_type: CoreType) Vtable {
    var vt: Vtable = .{};
    blas1.installGeneric(&vt);
    kernel.install(&vt, core_type);
    gemm_params.installCore(&vt, .x86_64, core_type, .fromConfig(blas.config));
    return vt;
}

test {
    _ = detector;
}
