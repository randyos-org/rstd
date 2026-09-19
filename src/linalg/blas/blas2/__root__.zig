//! BLAS Level 2 (matrix-vector) routines.
//!
//! One file per routine, holding the vtable wrappers for that routine.

const zigstd = @import("../../../zigstd.zig");

pub const gemv = @import("gemv/__root__.zig");
pub const ger = @import("ger/__root__.zig");
pub const hemv = @import("hemv/__root__.zig");
pub const symv = @import("symv/__root__.zig");

test {
    zigstd.testRefAllDecls(@This());
}
