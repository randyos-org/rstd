//! BLAS Level 1 (vector-vector) routines.
//!
//! One directory per routine, holding its portable kernels and tests.

const Vtable = @import("../Vtable.zig");

pub const asum = @import("asum/__root__.zig");
pub const axpby = @import("axpby/__root__.zig");
pub const axpy = @import("axpy/__root__.zig");
pub const copy = @import("copy/__root__.zig");
pub const dot = @import("dot/__root__.zig");
/// `amax`, `amin`, `iamax`, `iamin`, `imax`, `imin`, `max` and `min`.
pub const extremum = @import("extremum.zig");
pub const nrm2 = @import("nrm2/__root__.zig");
pub const rot = @import("rot/__root__.zig");
pub const rotm = @import("rotm/__root__.zig");
pub const scal = @import("scal/__root__.zig");
pub const sum = @import("sum/__root__.zig");
pub const swap = @import("swap/__root__.zig");

/// The BLAS calls themselves (`copy`, `dot`, `nrm2`, ...), over the kernels.
pub const interface = @import("interface.zig");

/// Every routine above, for `installGeneric`.
const routines = .{
    axpy,
    asum,
    axpby,
    copy,
    dot,
    extremum,
    nrm2,
    rot,
    rotm,
    scal,
    sum,
    swap,
    // the bfloat16 conversions live with the extensions; `sbdot` uses them
    @import("../ext/bf16/__root__.zig"),
};

/// Points the vtable at every level-1 routine's portable kernels, where it has
/// them; a core's own kernels go over these (`cpu.kernel.install`).
pub fn installGeneric(vt: *Vtable) void {
    inline for (routines) |r| {
        if (@hasDecl(r, "installGeneric")) r.installGeneric(vt);
    }
}

const zigstd = @import("../../../zigstd.zig");

test {
    zigstd.testRefAllDecls(@This());
}
