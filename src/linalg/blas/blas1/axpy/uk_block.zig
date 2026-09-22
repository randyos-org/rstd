//! The kernels whose upstream source is a driver around a microkernel block:
//! the shared driver in `kernel.zig` composed with the block an architecture's
//! `uk/__root__.zig` names for a core.
//!
//! An entry says the driver's block size, and either an assembly function of
//! a converted module, or that the block is C intrinsics upstream (written
//! on vectors by the architecture and passed in), or that the driver has no
//! block at all and the scalar loops do everything.

const std = @import("std");
const types = @import("../../types.zig");
const fma = @import("../../fma.zig");
const kernel = @import("kernel.zig");

const BlasLong = types.BlasLong;
const Complex = std.math.Complex;

/// What a driver calls for its contiguous part, as a block `uk/__root__.zig`
/// names (or an architecture's own entries, which use this type).
pub const Block = enum {
    /// A function of `module`, called as the reference's block kernels are.
    assembly,
    /// C intrinsics upstream, written on vectors by the architecture.
    vector,
    /// No block: the scalar loops do it all.
    none,
};

fn isComplex(comptime role: []const u8) bool {
    return role[0] == 'c' or role[0] == 'z';
}

/// `caxpyc`, `zaxpyc`.
fn isConj(comptime role: []const u8) bool {
    return role.len == 6;
}

/// The component type of `role`.
pub fn Float(comptime role: []const u8) type {
    return if (role[0] == 'd' or role[0] == 'z') f64 else f32;
}

/// The kernel type of `role`.
pub fn Kernel(comptime role: []const u8) type {
    return if (isComplex(role)) kernel.ComplexKernel(Float(role)) else kernel.RealKernel(Float(role));
}

/// The block that `e` names, for the real driver. `Vector` is the
/// architecture's block on vectors, for `e.block == .vector`.
fn RealBlock(comptime T: type, comptime e: anytype, comptime Vector: fn (comptime type) type) type {
    if (e.block == .vector) return Vector(T);
    return struct {
        pub fn run(n1: BlasLong, alpha: T, x: [*]T, y: [*]T) void {
            const a = alpha;
            @field(e.module, e.kernel)(n1, x, y, @ptrCast(&a));
        }
    };
}

/// The block that `e` names, for the complex driver: it takes the components
/// as an array.
fn CplxBlock(comptime F: type, comptime e: anytype) type {
    return struct {
        pub fn run(n1: BlasLong, alpha: Complex(F), x: [*]Complex(F), y: [*]Complex(F)) void {
            const a = [2]F{ alpha.re, alpha.im };
            @field(e.module, e.kernel)(n1, @ptrCast(x), @ptrCast(y), @ptrCast(&a));
        }
    };
}

/// The kernel `role` for the core whose entry is `e`.
pub fn driven(comptime role: []const u8, comptime e: anytype, comptime Vector: fn (comptime type) type) Kernel(role) {
    const F = Float(role);
    const mask = e.mask;
    const fuse = fma.fuse;
    if (comptime isComplex(role)) {
        const B = if (e.block == .none) null else CplxBlock(F, e);
        return kernel.complexEntry(F, kernel.Cplx(F, isConj(role), mask, B, fuse));
    }
    const B = if (e.block == .none) null else RealBlock(F, e, Vector);
    return kernel.realEntry(F, kernel.Real(F, mask, B, fuse));
}

/// For architectures with no block on vectors.
pub fn NoVector(comptime T: type) type {
    _ = T;
    return struct {};
}
