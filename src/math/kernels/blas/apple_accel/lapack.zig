//! Accelerate's LAPACK routines (`apple_accel.zig`).

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../types.zig");
const accel = @import("../apple_accel.zig");
const binding = @import("binding.zig");
const routines = @import("routines.zig");

const isComplex = types.isComplex;
const x86_64 = binding.x86_64;
const routine = routines.routine;
const covers = accel.covers;
const use = accel.use;
const mat = accel.mat;
const calls = &accel.calls;

/// `?potrf_`'s `INFO`, or null where rstd runs it. Accelerate's LAPACK runs
/// on its BLAS, single-threaded as `load` sets it.
pub inline fn potrf(comptime T: type, ul: types.Uplo, n: isize, a: [*]T, lda: isize) ?isize {
    if (comptime !covers(T) or (x86_64 and isComplex(T))) return null else {
        const b = use(T) orelse return null;
        const u: u8 = switch (ul) {
            .u => 'U',
            .l => 'L',
            .v, .m => return null,
        };
        switch (b.*) {
            inline else => |*bound, tag| {
                const f = @field(bound.fns, @tagName(routine(.potrf, T))) orelse return null;
                const I = if (tag == .ilp64) i64 else c_int;
                if (I != i64 and mat(n, n, lda) > std.math.maxInt(I)) return null;
                const nn = std.math.cast(I, n) orelse return null;
                const ld = std.math.cast(I, lda) orelse return null;
                var info: I = 0;
                if (sysinfo.is_test) _ = calls.fetchAdd(1, .monotonic);
                f(&u, &nn, @ptrCast(a), &ld, &info);
                return info;
            },
        }
    }
}
