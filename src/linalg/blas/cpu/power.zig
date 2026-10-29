//! POWER CPU detection (`_detect/power.zig`), with its kernel registry
//! and vtable. The detector is its own file so the build can use it
//! without importing the kernels.

const sysinfo = @import("builtin");
const detector = @import("_detect/power.zig");
const Vtable = @import("../Vtable.zig");
const portable = @import("../portable.zig");
const blas = @import("../__root__.zig");
const gemm_params = @import("../blas3/gemm/params.zig");

pub const CoreType = detector.CoreType;
pub const kernel = @import("../kernels/power.zig");
pub const detect = detector.detect;
pub const canRun = detector.supportsCore;

/// The blocking tables of this word size and byte order.
const family: gemm_params.Family = switch (sysinfo.target.cpu.arch) {
    .powerpc64le => .power64,
    .powerpc64 => .power64be,
    .powerpcle => .power32le,
    else => .power32,
};

pub fn vtableFromCoreType(core_type: CoreType) Vtable {
    var vt: Vtable = .{};
    portable.install(&vt);
    // blocking first: a core's level-3 kernels set their own unrolls
    gemm_params.installCore(&vt, family, core_type, .fromConfig(blas.config));
    kernel.install(&vt, core_type);
    portable.finish(&vt);
    return vt;
}

test {
    _ = detector;
}
