//! Accelerate's LAPACK routines (`apple_accel.zig`).

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../types.zig");
const accel = @import("../apple_accel.zig");
const binding = @import("binding.zig");
const routines = @import("routines.zig");
const quirks = @import("quirks.zig");

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

/// `?getrf_`'s `INFO`, `ipiv` 0-based, or null where rstd runs it: x86_64,
/// whose Accelerate trails rstd (issue #772), an LP64 table, and an `A` with
/// NaN or an infinity, as Accelerate's `i?amax` does not take the first NaN
/// that rstd's pivot search takes (issue #840).
pub inline fn getrf(comptime T: type, m: isize, n: isize, a: [*]T, lda: isize, ipiv: []usize) ?isize {
    if (comptime !covers(T) or x86_64) return null else {
        const b = use(T) orelse return null;
        const bound = switch (b.*) {
            .ilp64 => |*x| x,
            .lp64 => return null,
        };
        const f = @field(bound.fns, @tagName(routine(.getrf, T))) orelse return null;
        if (quirks.nonFiniteMat(T, m, n, a, lda)) return null;
        const p: [*]i64 = @ptrCast(ipiv.ptr);
        var info: i64 = 0;
        if (sysinfo.is_test) _ = calls.fetchAdd(1, .monotonic);
        const mi: i64 = m;
        const ni: i64 = n;
        const ld: i64 = lda;
        f(&mi, &ni, @ptrCast(a), &ld, p, &info);
        for (ipiv) |*v| v.* -%= 1;
        return info;
    }
}

/// `?getrs_` with `getrf`'s 0-based `ipiv`; false where rstd runs it: as
/// `getrf`, where `trsm` keeps rstd's (`B` not finite, a tiny pivot),
/// and for `z`.
pub fn getrs(allocator: std.mem.Allocator, comptime T: type, trans: types.Trans, n: isize, nrhs: isize, a: [*]const T, lda: isize, ipiv: []const usize, b_: [*]T, ldb: isize) std.mem.Allocator.Error!bool {
    // `z` too: rstd's threaded solve runs faster (on the M6 at 12 threads
    // routed 0.63 to 0.86 of OpenBLAS alone, its own 0.86 to 1.18)
    if (comptime !covers(T) or x86_64 or T == types.Complex64) return false else {
        const b = use(T) orelse return false;
        const bound = switch (b.*) {
            .ilp64 => |*x| x,
            .lp64 => return false,
        };
        const f = @field(bound.fns, @tagName(routine(.getrs, T))) orelse return false;
        const t: u8 = switch (trans) {
            .n => 'N',
            .t => 'T',
            .c => if (comptime isComplex(T)) 'C' else 'T',
            else => return false,
        };
        if (quirks.nonFiniteMat(T, n, nrhs, b_, ldb) or quirks.tinyDiagonal(T, n, a, lda)) return false;
        const p = try allocator.alloc(i64, ipiv.len);
        defer allocator.free(p);
        for (p, ipiv) |*o, v| o.* = @intCast(v + 1);
        var info: i64 = 0;
        if (sysinfo.is_test) _ = calls.fetchAdd(1, .monotonic);
        const ni: i64 = n;
        const ri: i64 = nrhs;
        const la: i64 = lda;
        const lb: i64 = ldb;
        f(&t, &ni, &ri, @ptrCast(a), &la, p.ptr, @ptrCast(b_), &lb, &info);
        return true;
    }
}
