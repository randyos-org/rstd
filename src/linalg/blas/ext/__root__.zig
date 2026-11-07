//! BLAS extensions and auxiliary kernels (bf16 conversion, geadd, matcopy, laswp, neg copies).
//!
//! One directory per routine, holding its portable kernels and their tests;
//! `interface.zig` is the checks in front of the kernels.

pub const bf16 = @import("bf16/__root__.zig");
pub const geadd = @import("geadd/__root__.zig");
pub const laswp = @import("laswp/__root__.zig");
pub const neg = @import("neg/__root__.zig");
/// The portable `omatcopy` and `imatcopy` kernels, which share their scaling;
/// column-major only, the row-major `ORDER` `R` ones are left out.
pub const matcopy = @import("matcopy.zig");

/// The extension calls themselves (`geadd`, `omatcopy`, `imatcopy`,
/// `laswp`), over the kernels.
pub const interface = @import("interface.zig");

const Vtable = @import("../Vtable.zig");

/// Points the vtable at the extensions' portable kernels; the bfloat16 ones
/// go in with level 1, and `neg_tcopy` and `laswp_ncopy` with the gemm
/// kernels.
pub fn installGeneric(vt: *Vtable) void {
    geadd.installGeneric(vt);
    matcopy.installGeneric(vt);
}

const zigstd = @import("../../../zigstd.zig");

test {
    zigstd.testRefAllDecls(@This());
}
