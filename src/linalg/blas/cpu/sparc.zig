//! SPARC CPU detection (`_detect/sparc.zig`), with its kernel registry
//! and vtable. The detector is its own file so the build can use it
//! without importing the kernels.

const detector = @import("_detect/sparc.zig");
const Vtable = @import("../Vtable.zig");
const portable = @import("../portable.zig");
const blas = @import("../__root__.zig");
const gemm_params = @import("../blas3/gemm/params.zig");

pub const CoreType = detector.CoreType;
pub const kernel = @import("../kernels/sparc.zig");
pub const detect = detector.detect;

pub fn vtableFromCoreType(core_type: CoreType) Vtable {
    var vt: Vtable = .{};
    portable.install(&vt);
    // blocking first: a core's level-3 kernels set their own unrolls
    gemm_params.installCore(&vt, .sparc, core_type, .fromConfig(blas.config));
    kernel.install(&vt, core_type);
    portable.finish(&vt);
    return vt;
}

test {
    _ = detector;
}
