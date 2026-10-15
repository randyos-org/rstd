//! BLAS Level 2 (matrix-vector) routines.
//!
//! One file per routine, holding the vtable wrappers for that routine.

const zigstd = @import("../../../zigstd.zig");

pub const gemv = @import("gemv/__root__.zig");
pub const ger = @import("ger/__root__.zig");
pub const hemv = @import("hemv/__root__.zig");
pub const symv = @import("symv/__root__.zig");

const Vtable = @import("../Vtable.zig");

/// Points the vtable at every level-2 routine's portable kernels, where it has
/// them.
pub fn installGeneric(vt: *Vtable) void {
    inline for (.{ gemv, ger, hemv, symv }) |r| {
        if (@hasDecl(r, "installGeneric")) r.installGeneric(vt);
    }
}

test {
    zigstd.testRefAllDecls(@This());
}
