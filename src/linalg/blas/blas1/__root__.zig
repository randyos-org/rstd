//! BLAS Level 1 (vector-vector) routines.
//!
//! One file per routine, holding the vtable wrappers for that routine.

const zigstd = @import("../../../zigstd.zig");

pub const amax = @import("amax/__root__.zig");
pub const amin = @import("amin/__root__.zig");
pub const asum = @import("asum/__root__.zig");
pub const axpby = @import("axpby/__root__.zig");
pub const axpy = @import("axpy/__root__.zig");
pub const copy = @import("copy/__root__.zig");
pub const dot = @import("dot/__root__.zig");
pub const iamax = @import("iamax/__root__.zig");
pub const iamin = @import("iamin/__root__.zig");
pub const imax = @import("imax/__root__.zig");
pub const imin = @import("imin/__root__.zig");
pub const max = @import("max/__root__.zig");
pub const min = @import("min/__root__.zig");
pub const nrm2 = @import("nrm2/__root__.zig");
pub const rot = @import("rot/__root__.zig");
pub const rotm = @import("rotm/__root__.zig");
pub const scal = @import("scal/__root__.zig");
pub const sum = @import("sum/__root__.zig");
pub const swap = @import("swap/__root__.zig");

test {
    zigstd.testRefAllDecls(@This());
}
