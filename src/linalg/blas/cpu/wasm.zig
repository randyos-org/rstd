//! WebAssembly CPU detection (`_detect/wasm.zig`), with its kernel registry
//! and vtable. The detector is its own file so the build can use it
//! without importing the kernels.

const detector = @import("_detect/wasm.zig");
const Vtable = @import("../Vtable.zig");
const portable = @import("../portable.zig");
const blas = @import("../__root__.zig");
const gemm_params = @import("../blas3/gemm/params.zig");

pub const CoreType = detector.CoreType;
pub const kernel = @import("../kernels/wasm.zig");
pub const detect = detector.detect;

pub fn vtableFromCoreType(core_type: CoreType) Vtable {
    var vt: Vtable = .{};
    portable.install(&vt);
    kernel.install(&vt, core_type);
    gemm_params.installCore(&vt, .wasm, core_type, .fromConfig(blas.config));
    return vt;
}

test {
    _ = detector;
}
