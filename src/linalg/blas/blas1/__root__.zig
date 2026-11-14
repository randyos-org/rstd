//! BLAS Level 1 (vector-vector) routines.
//!
//! One directory per routine, holding its portable kernels and tests.

const std = @import("std");
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

/// The bfloat16 conversions live with the extensions; `sbdot` uses them.
const bf16 = @import("../ext/bf16/__root__.zig");

/// Points the vtable at every level-1 routine's portable kernels, where it has
/// them; a core's own kernels go over these (`cpu.kernel.install`).
pub fn installGeneric(vt: *Vtable) void {
    // the public routine modules, by reflection
    inline for (comptime std.meta.declarations(@This())) |d| {
        const r = @field(@This(), d);
        if (@TypeOf(r) == type) {
            if (@hasDecl(r, "installGeneric")) r.installGeneric(vt);
        }
    }
    bf16.installGeneric(vt);
}

const zigstd = @import("../../../zigstd.zig");

test {
    zigstd.testRefAllDecls(@This());
}
