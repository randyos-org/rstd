//! OpenBLAS's x86_64 level-1 kernels for HASWELL and ZEN, which
//! `build/bench.zig` compiles from the `-Dopenblas` checkout as
//! `_ob_<core>_<slot>`: the `l1` section's references on x86_64.

const rstd = @import("rstd");
const bb = @import("blas_bench.zig");

const CoreType = rstd.linalg.blas.cpu.CoreType;
const Fn = bb.Fn;

pub const Core = enum { haswell, zen };

/// The OpenBLAS core built for rstd's core `c`, if any.
pub fn coreFor(c: CoreType) ?Core {
    return switch (c) {
        .haswell => .haswell,
        .zen => .zen,
        else => null,
    };
}

/// OpenBLAS's kernel for vtable slot `name` on `core`.
pub fn kernel(comptime name: []const u8, core: Core) Fn(name) {
    return switch (core) {
        inline else => |k| @extern(Fn(name), .{ .name = "_ob_" ++ @tagName(k) ++ "_" ++ name }),
    };
}
