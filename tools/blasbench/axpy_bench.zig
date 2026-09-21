//! `daxpy` implementations raced against OpenBLAS's own assembled kernel.
//!
//! Build (from the rstd root), assembling the reference first:
//!
//!     zig cc -c -x assembler-with-cpp -w -I <inc> -I <openblas> \
//!       -I <openblas>/kernel/x86_64 -DASSEMBLER -DARCH_X86_64 -DOS_WINNT \
//!       -DPENRYN -DHAVE_SSE -DHAVE_SSE2 -DHAVE_SSE3 -DDOUBLE \
//!       -DASMNAME=ob_daxpy_k -DASMFNAME=ob_daxpy_k_ -DNAME=ob_daxpy_k_ \
//!       -DCNAME=ob_daxpy_k <openblas>/kernel/x86_64/axpy_sse2.S -o ref.o
//!     zig build-exe axpy_bench.zig ref.o -OReleaseFast -lc

const std = @import("std");
const bench = @import("bench.zig");

/// OpenBLAS's kernel, assembled from the `.S` and linked in.
extern fn ob_daxpy_k(
    n: isize,
    dummy0: isize,
    dummy1: isize,
    alpha: f64,
    x: [*]f64,
    inc_x: isize,
    y: [*]f64,
    inc_y: isize,
    dummy: [*]f64,
    dummy2: isize,
) callconv(.c) i32;

/// OpenBLAS's Haswell kernel: `daxpy.c` around `daxpy_microk_haswell-2.c`.
extern fn ob_hsw_daxpy_k(
    n: isize,
    dummy0: isize,
    dummy1: isize,
    alpha: f64,
    x: [*]f64,
    inc_x: isize,
    y: [*]f64,
    inc_y: isize,
    dummy: [*]f64,
    dummy2: isize,
) callconv(.c) i32;

const gen = @import("gen_kernel");

/// The generated kernel: OpenBLAS's body, inline in a Zig function.
fn gen_daxpy(n: isize, _: isize, _: isize, alpha: f64, x: [*]f64, inc_x: isize, y: [*]f64, inc_y: isize, _: [*]f64, _: isize) callconv(.c) i32 {
    gen.daxpy(n, alpha, x, inc_x, y, inc_y);
    return 0;
}

fn scalar(n: isize, _: isize, _: isize, alpha: f64, x: [*]f64, _: isize, y: [*]f64, _: isize, _: [*]f64, _: isize) callconv(.c) i32 {
    const len: usize = @intCast(n);
    for (0..len) |i| y[i] += alpha * x[i];
    return 0;
}

/// Two doubles at a time, the width `axpy_sse2.S` works in.
fn vector2(n: isize, _: isize, _: isize, alpha: f64, x: [*]f64, _: isize, y: [*]f64, _: isize, _: [*]f64, _: isize) callconv(.c) i32 {
    return vectorN(2, n, alpha, x, y);
}

/// Four doubles at a time, what AVX gives if the target allows it.
fn vector4(n: isize, _: isize, _: isize, alpha: f64, x: [*]f64, _: isize, y: [*]f64, _: isize, _: [*]f64, _: isize) callconv(.c) i32 {
    return vectorN(4, n, alpha, x, y);
}

/// Eight at a time, unrolled by two vectors.
fn vector8(n: isize, _: isize, _: isize, alpha: f64, x: [*]f64, _: isize, y: [*]f64, _: isize, _: [*]f64, _: isize) callconv(.c) i32 {
    return vectorN(8, n, alpha, x, y);
}

fn vectorN(comptime w: usize, n: isize, alpha: f64, x: [*]f64, y: [*]f64) i32 {
    const V = @Vector(w, f64);
    const av: V = @splat(alpha);
    const len: usize = @intCast(n);
    var i: usize = 0;
    while (i + w <= len) : (i += w) {
        const xv: V = x[i..][0..w].*;
        const yv: V = y[i..][0..w].*;
        y[i..][0..w].* = yv + av * xv;
    }
    while (i < len) : (i += 1) y[i] += alpha * x[i];
    return 0;
}

pub fn main(init: std.process.Init) !void {
    const gpa = init.arena.allocator();

    const impls = [_]bench.Impl{
        .{ .name = "ob sse2.S", .run = &ob_daxpy_k },
        .{ .name = "ob haswell", .run = &ob_hsw_daxpy_k },
        .{ .name = "generated", .run = &gen_daxpy },
        .{ .name = "zig scalar", .run = &scalar },
        .{ .name = "zig vec2", .run = &vector2 },
        .{ .name = "zig vec4", .run = &vector4 },
        .{ .name = "zig vec8", .run = &vector8 },
    };
    // in L1, in L2, in L3, and past it into memory
    const sizes = [_]usize{ 256, 1024, 8 * 1024, 64 * 1024, 512 * 1024, 4 * 1024 * 1024 };
    try bench.sweep(init.io, gpa, &impls, &sizes, bench.axpy_work, 200 * std.time.ns_per_ms);
}
