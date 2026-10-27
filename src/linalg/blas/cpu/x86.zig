//! 32-bit x86 CPU detection: the detector shared with x86_64 (`_detect/x86.zig`),
//! with 32-bit x86's own kernel registry and vtable.

const detector = @import("_detect/x86.zig");
const Vtable = @import("../Vtable.zig");
const portable = @import("../portable.zig");
const blas = @import("../__root__.zig");
const gemm_params = @import("../blas3/gemm/params.zig");
const gemv = @import("../blas2/gemv/__root__.zig");

pub const CoreType = detector.CoreType;
pub const kernel = @import("../kernels/x86.zig");
pub const detect = detector.detect;
/// Whether this machine can run `core`'s kernels (for tests).
pub const canRun = detector.supportsCore;

pub fn vtableFromCoreType(core_type: CoreType) Vtable {
    var vt: Vtable = .{};
    portable.install(&vt);
    // blocking first: a core's level-3 kernels set their own unrolls
    gemm_params.installCore(&vt, .x86, core_type, .fromConfig(blas.config));
    kernel.install(&vt, core_type);
    // the x87 `q`/`x` gemv_n block by the core's DTB entries
    gemv.installExtendedCore(&vt, core_type);
    portable.finish(&vt);
    return vt;
}

test {
    _ = detector;
}
