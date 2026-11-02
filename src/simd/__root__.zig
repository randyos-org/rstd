//! SIMD widths and reductions. The vector math (`mix`, `mulAdd`, `dot` and
//! the rest) is `linalg.vector`'s.

const std = @import("std");

pub const zigvec = @import("zigvec.zig");
pub const reduce = @import("reduce.zig");

/// Preferred SIMD vector length (element count, not bits) for `T` on the
/// compile target -- e.g. 2 for `f64` on plain SSE2, 8 once AVX-512 is in
/// play. Thin wrapper over `std.simd.suggestVectorLength`, falling back to
/// scalar (1 lane) where it has no suggestion, so every caller gets a
/// usable `@Vector(vectorLen(T), T)` without its own null-handling.
///
/// This is a *compile-time* choice tied to whatever CPU model the build was
/// configured for (Zig's default x86_64 "baseline" model only guarantees
/// SSE2, so this returns 2 for `f64` unless `-Dcpu=` says otherwise) -- not
/// runtime dispatch. A binary built this way hits an illegal instruction on
/// hardware narrower than what it assumes, exactly like a microarch-locked
/// OpenBLAS `TARGET=SKYLAKEX` build. Runtime "detect the CPU and pick a
/// kernel" (`DYNAMIC_ARCH`) is a separate, deliberately deferred problem.
pub fn vectorLen(comptime T: type) comptime_int {
    return comptime std.simd.suggestVectorLength(T) orelse 1;
}

const zigstd = @import("../zigstd.zig");

test {
    zigstd.testRefAllDecls(@This());
}
