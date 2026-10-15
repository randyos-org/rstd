//! Every level's portable kernels, which the vtable holds wherever a core has
//! none of its own; a core's kernels then go over them (`cpu.kernel.install`).

const Vtable = @import("Vtable.zig");

const levels = .{
    @import("blas1/__root__.zig"),
    @import("blas2/__root__.zig"),
};

pub fn install(vt: *Vtable) void {
    inline for (levels) |l| l.installGeneric(vt);
    alias(vt);
}

/// Points the bfloat16-named slots with no kernel of their own at the `f32`
/// kernels `vt` holds (`ext/bf16/kernel.zig`); again after a core's kernels.
pub const alias = @import("ext/bf16/kernel.zig").alias;
