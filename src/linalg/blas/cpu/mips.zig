//! MIPS CPU detection (`_detect/mips.zig`), with its kernel registry
//! and vtable. The detector is its own file so the build can use it
//! without importing the kernels.

const detector = @import("_detect/mips.zig");
const Vtable = @import("../Vtable.zig");
const portable = @import("../portable.zig");
const blas = @import("../__root__.zig");
const gemm_params = @import("../blas3/gemm/params.zig");

pub const CoreType = detector.CoreType;
pub const kernel = @import("../kernels/mips.zig");
pub const detect = detector.detect;
/// Whether this machine can run `core`'s kernels (for tests).
pub const canRun = detector.canRun;

pub fn vtableFromCoreType(core_type: CoreType) Vtable {
    var vt: Vtable = .{};
    portable.install(&vt);
    kernel.install(&vt, core_type);
    gemm_params.installCore(&vt, (if (@sizeOf(usize) == 8) .mips64 else .mips32), core_type, .fromConfig(blas.config));
    return vt;
}

test {
    _ = detector;
}
