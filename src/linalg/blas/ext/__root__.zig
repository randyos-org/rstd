//! BLAS extensions and auxiliary kernels (bf16 conversion, geadd, matcopy, laswp, neg copies).
//!
//! One file per routine, holding the vtable wrappers for that routine.

pub const bf16 = @import("bf16/__root__.zig");
pub const geadd = @import("geadd/__root__.zig");
pub const imatcopy = @import("imatcopy/__root__.zig");
pub const laswp = @import("laswp/__root__.zig");
pub const neg = @import("neg/__root__.zig");
pub const omatcopy = @import("omatcopy/__root__.zig");

const zigstd = @import("../../../zigstd.zig");

test {
    zigstd.testRefAllDecls(@This());
}
