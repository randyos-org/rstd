//! BLAS Level 2 (matrix-vector) routines.
//!
//! One directory per kernel-backed routine, holding its vtable wrappers and
//! portable kernels; `interface.zig`, `triangular.zig` and `banded.zig` are the
//! checks and drivers in front of them (the OpenBLAS `interface/` and
//! `driver/level2/` code).

pub const banded = @import("banded.zig");
pub const gemv = @import("gemv/__root__.zig");
pub const ger = @import("ger/__root__.zig");
pub const hemv = @import("hemv/__root__.zig");
pub const interface = @import("interface.zig");
pub const symv = @import("symv/__root__.zig");
pub const triangular = @import("triangular.zig");

const Vtable = @import("../Vtable.zig");

/// Points the vtable at every level-2 routine's portable kernels, where it has
/// them.
pub fn installGeneric(vt: *Vtable) void {
    inline for (.{ gemv, ger, hemv, symv }) |r| {
        if (@hasDecl(r, "installGeneric")) r.installGeneric(vt);
    }
}

const zigstd = @import("../../../zigstd.zig");

test {
    zigstd.testRefAllDecls(@This());
}
