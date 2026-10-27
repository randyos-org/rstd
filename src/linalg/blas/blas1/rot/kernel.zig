//! The drivers of `kernel/x86_64/srot.c` and `drot.c` (`rot_compute`): unit
//! strides go to the block (`<p>rot_kernel`), others to the C loop. Each
//! element is `x' = c*x + s*y`, `y' = c*y - s*x`; the blocks are intrinsics
//! (`Vector(slot)` of the architecture) or, where the chain names none, the
//! C block on `../simd/intrin.h` here. Nothing here depends on the BLAS root,
//! so the per-core objects (`_kernel.zig`) can use it.
//!
//! The intrinsics fuse by name (`_mm256_fmadd_pd`: `@mulAdd`); the C
//! statements contract as the core's flags say (`fma.madd`), `c*x + s*y`
//! fusing the first product, so on the FMA cores every path gives one value.

const std = @import("std");
const sysinfo = @import("builtin");
const util = @import("../../util.zig");
const types = @import("../../types.zig");
const fma = @import("../../fma.zig");
const uk = @import("../uk_block.zig");

const BlasLong = types.BlasLong;
const Slot = util.Slot;

/// The kernels an architecture's `<arch>/kernel.zig` iterates over: the real
/// ones driven here, and the complex ones only converted assembly fills.
pub const slots = [_]Slot{ .srot_k, .drot_k, .csrot_k, .zdrot_k };

/// The element type of `slot`, one of the real `slots`.
pub fn Elem(comptime slot: Slot) type {
    return switch (slot) {
        .srot_k => f32,
        .drot_k => f64,
        else => @compileError(@tagName(slot) ++ " is not a real rot kernel"),
    };
}

/// Vtable slot type of `slot`.
pub fn Kernel(comptime slot: Slot) type {
    return util.SlotFn(slot);
}

/// One element: `c*x + s*y` into `x`, `c*y - s*x` into `y`, with the products
/// `s*y` and `s*x` rounded and each fused into the other (`fuse`); `y` is
/// stored first, as the C's loop does.
pub inline fn rotate(comptime F: type, comptime fuse: bool, x: *F, y: *F, c: F, s: F) void {
    const xv = x.*;
    const yv = y.*;
    const temp = fma.madd(fuse, F, c, xv, s * yv);
    y.* = fma.madd(fuse, F, c, yv, -(s * xv));
    x.* = temp;
}

/// `count` vectors of `V` from `i` on, rotated as the intrinsics do:
/// `t = s*y; t = fmadd(c, x, t)` into `x`, `t = s*x; t = fmsub(c, y, t)` into `y`.
pub inline fn rotateVectors(comptime V: type, comptime count: comptime_int, x: [*]@typeInfo(V).vector.child, y: [*]@typeInfo(V).vector.child, i: usize, c: V, s: V) void {
    const l = @typeInfo(V).vector.len;
    var xv: [count]V = undefined;
    var yv: [count]V = undefined;
    inline for (0..count) |k| {
        xv[k] = util.load(V, x + i + k * l);
        yv[k] = util.load(V, y + i + k * l);
    }
    inline for (0..count) |k| util.store(V, x + i + k * l, @mulAdd(V, c, xv[k], s * yv[k]));
    inline for (0..count) |k| util.store(V, y + i + k * l, @mulAdd(V, c, yv[k], -(s * xv[k])));
}

/// `V_SIMD` of `../simd/intrin.h` for the core being built: 512 with
/// AVX-512VL, 256 with AVX2, 128 with SSE2.
const v_simd: comptime_int = blk: {
    if (!sysinfo.target.cpu.arch.isX86()) break :blk 0;
    const f = sysinfo.target.cpu.features;
    const x86 = std.Target.x86;
    if (x86.featureSetHas(f, .avx512vl)) break :blk 512;
    if (x86.featureSetHas(f, .avx2)) break :blk 256;
    if (x86.featureSetHas(f, .sse2)) break :blk 128;
    break :blk 0;
};
const has_fma3 = sysinfo.target.cpu.arch.isX86() and std.Target.x86.featureSetHas(sysinfo.target.cpu.features, .fma);

/// The drivers' own `<p>rot_kernel`, for a core whose chain names no block:
/// vectors of `V_SIMD` bits four at a time, then one, then the C tail; or,
/// where the C takes no vectors, its unrolled scalar loop.
fn Fallback(comptime F: type, comptime fuse: bool) type {
    const vector = if (F == f64) v_simd > 256 else v_simd != 0 and (has_fma3 or v_simd > 128);
    return struct {
        pub fn run(n: BlasLong, x: [*]F, y: [*]F, c: F, s: F) void {
            const len: usize = @intCast(n);
            var i: usize = 0;
            if (comptime vector) {
                const l = v_simd / @bitSizeOf(F);
                const V = @Vector(l, F);
                const cv: V = @splat(c);
                const sv: V = @splat(s);
                const x4 = len & ~@as(usize, 4 * l - 1);
                const x1 = len & ~@as(usize, l - 1);
                while (i < x4) : (i += 4 * l) rotateVectors(V, 4, x, y, i, cv, sv);
                while (i < x1) : (i += l) rotateVectors(V, 1, x, y, i, cv, sv);
            } else {
                const n1 = len & ~@as(usize, 7);
                while (i < n1) : (i += 4) {
                    inline for (0..4) |k| rotate(F, fuse, &x[i + k], &y[i + k], c, s);
                }
            }
            while (i < len) : (i += 1) rotate(F, fuse, &x[i], &y[i], c, s);
        }
    };
}

/// `CNAME` (`rot_compute` on one thread) of `slot` around `B`.
fn Driver(comptime slot: Slot, comptime B: type, comptime fuse: bool) type {
    const F = Elem(slot);
    return struct {
        pub fn call(n: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, c: F, s: F) callconv(.c) i32 {
            if (n <= 0) return 0;
            if (inc_x == 1 and inc_y == 1) {
                B.run(n, x, y, c, s);
                return 0;
            }
            var ix: BlasLong = 0;
            var iy: BlasLong = 0;
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                rotate(F, fuse, util.at(F, x, ix), util.at(F, y, iy), c, s);
                ix += inc_x;
                iy += inc_y;
            }
            return 0;
        }
    };
}

/// The kernel `slot` of the core whose index entry is `e` (`uk_block.driven`).
pub fn entryOf(comptime slot: Slot, comptime e: uk.Entry, comptime Vector: fn (comptime Slot) type) Kernel(slot) {
    const B = switch (e.block) {
        .vector => Vector(slot),
        .none => Fallback(Elem(slot), fma.fuse),
        .assembly => @compileError("no rot microkernel is assembly"),
    };
    return &Driver(slot, B, fma.fuse).call;
}
