//! s390x CPU detection (`_detect/s390x.zig`), with its kernel registry
//! and vtable. The detector is its own file so the build can use it
//! without importing the kernels.

const detector = @import("_detect/s390x.zig");
const Vtable = @import("../Vtable.zig");
const portable = @import("../portable.zig");
const blas = @import("../__root__.zig");
const gemm_params = @import("../blas3/gemm/params.zig");

pub const CoreType = detector.CoreType;
pub const canRun = detector.supportsCore;
pub const kernel = @import("../kernels/s390x.zig");
/// Detects the running core, under the configured `LD_HWCAP_MASK`.
pub fn detect() CoreType {
    return detector.detectMasked(blas.config.ld_hwcap_mask orelse ~@as(u64, 0));
}

pub fn vtableFromCoreType(core_type: CoreType) Vtable {
    var vt: Vtable = .{};
    portable.install(&vt);
    // blocking first: a core's level-3 kernels set their own unrolls
    gemm_params.installCore(&vt, .s390x, core_type, .fromConfig(blas.config));
    kernel.install(&vt, core_type);
    portable.finish(&vt);
    return vt;
}

test {
    _ = detector;
}
