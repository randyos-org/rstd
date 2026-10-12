//! The kernels whose upstream source is a C driver around a microkernel
//! block, for any routine: its driver, ported, composed with the block an
//! architecture's `uk/__root__.zig` names for a core.
//!
//! An entry says the driver's block size, and either an assembly function of
//! a converted module, or that the block is C intrinsics upstream (written
//! on vectors by the architecture and passed in), or that the driver has no
//! block at all and the scalar loops do everything.

/// What a driver calls for its contiguous part.
pub const Block = enum {
    /// A function of `module`, called as the reference's block kernels are.
    assembly,
    /// C intrinsics upstream, written on vectors by the architecture.
    vector,
    /// No block: the scalar loops do it all.
    none,
};

/// Whether `role` (`caxpy`, `zdot`, ...) is complex.
pub fn isComplex(comptime role: []const u8) bool {
    return role[0] == 'c' or role[0] == 'z';
}

/// The component type of `role`.
pub fn Float(comptime role: []const u8) type {
    return if (role[0] == 'd' or role[0] == 'z') f64 else f32;
}

/// The kernel `role` of the core whose entry is `e`. `Driver` is the routine's
/// ported driver: `Kernel(role)` is the slot type, `Asm(role, f)` a block that
/// calls assembly function `f` as the reference does, and `entry(role, mask, B)`
/// the kernel around block `B` (null for none). `Vector(role)` is the
/// architecture's block on vectors.
pub fn driven(
    comptime Driver: type,
    comptime role: []const u8,
    comptime e: anytype,
    comptime Vector: fn (comptime []const u8) type,
) Driver.Kernel(role) {
    const B: ?type = switch (e.block) {
        .none => null,
        .vector => Vector(role),
        .assembly => Driver.Asm(role, @field(e.module, e.kernel)),
    };
    return Driver.entry(role, e.mask, B);
}

/// For architectures with no block on vectors.
pub fn NoVector(comptime role: []const u8) type {
    _ = role;
    return struct {};
}
