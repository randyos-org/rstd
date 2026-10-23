const std = @import("std");
const sysinfo = @import("builtin");

/// Equivalent of C's `base[offset]` where `offset` may be negative -- BLAS
/// routines allow negative `inc_x`/`inc_y` (iterate a vector in reverse),
/// which does not map onto Zig's `usize`-indexed pointer/slice arithmetic
/// the way it does onto a C pointer (`x[ix]` with negative `ix` is just
/// `*(x + ix)` in C).
pub inline fn at(comptime T: type, base: [*]T, offset: i64) *T {
    const addr: usize = @intFromPtr(base);
    const byte_offset: i64 = offset * @as(i64, @sizeOf(T));
    const new_addr: usize = if (byte_offset >= 0)
        addr + @as(usize, @intCast(byte_offset))
    else
        addr - @as(usize, @intCast(-byte_offset));
    return @ptrFromInt(new_addr);
}

/// The symbol of the kernel `role` (`daxpy`, ...) of `core`. The per-core objects
/// (`build/blaskern.zig`) export it, with hidden visibility so it never leaves a
/// shared library rstd is built into, and `kernels/<arch>.zig`'s `install` looks
/// it up by the same name.
pub fn kernelSymbol(comptime role: []const u8, comptime core: []const u8) []const u8 {
    return "blaskern_" ++ role ++ "_" ++ core;
}

/// Loads a `V` from `p`, which is not assumed aligned to `V`'s natural
/// alignment -- e.g. a kernel walking a caller's buffer by an arbitrary
/// stride, where nothing guarantees the vector-width alignment `@Vector`
/// otherwise wants. (Lives here, not in `simd/`, because the per-core kernel
/// objects that use it are their own module rooted at `_kernel.zig` and
/// cannot reach outside `linalg/blas/`; see `kernels/__root__.zig`.)
pub inline fn load(comptime V: type, p: [*]const @typeInfo(V).vector.child) V {
    return @as(*align(1) const V, @ptrCast(p)).*;
}

/// Stores `v` to `p`, unaligned; see `load`.
pub inline fn store(comptime V: type, p: [*]@typeInfo(V).vector.child, v: V) void {
    @as(*align(1) V, @ptrCast(p)).* = v;
}

/// Whether `@Vector` arithmetic in `F` keeps denormals on this target. Not
/// `f32` on ARMv7 NEON, which flushes to zero whatever FPSCR says, nor on
/// AltiVec without VSX, which flushes under `VSCR[NJ]`, as Linux and QEMU
/// start it. Neither has `f64` lanes; LLVM splits those into scalar IEEE
/// operations, as it does every vector under a soft-float ABI. Kernels take
/// their scalar path where this is false (`comparisons.md`, 10).
pub fn vectorIeee(comptime F: type) bool {
    if (F != f32 or sysinfo.abi.float() == .soft) return true;
    const features = sysinfo.cpu.features;
    return switch (sysinfo.cpu.arch) {
        .arm, .armeb, .thumb, .thumbeb => !std.Target.arm.featureSetHas(features, .neon),
        .powerpc, .powerpcle, .powerpc64, .powerpc64le => !std.Target.powerpc.featureSetHas(features, .altivec) or
            std.Target.powerpc.featureSetHas(features, .vsx),
        else => true,
    };
}

/// Whether a vector of `n` elements of `size` bytes at stride `inc` fits in
/// the address space: `((n - 1) * |inc| + 1) * size` bytes, with no overflow.
/// Where it does not, no such vector can exist; the interfaces refuse it so
/// no index computation overflows.
pub fn fits(n: isize, inc: isize, size: usize) bool {
    if (inc == std.math.minInt(isize)) return false;
    if (n <= 1) return true;
    return fitsMatrix(1, n, @intCast(@abs(inc)), size);
}

/// Whether a column-major `rows` by `cols` matrix with leading dimension `ld`
/// fits, as `fits` does for a vector: `(cols - 1) * ld + rows` elements.
pub fn fitsMatrix(rows: isize, cols: isize, ld: isize, size: usize) bool {
    if (rows <= 0 or cols <= 0) return true;
    const span = std.math.mul(isize, cols - 1, ld) catch return false;
    const elems = std.math.add(isize, span, rows) catch return false;
    _ = std.math.mul(isize, elems, std.math.cast(isize, size) orelse return false) catch return false;
    return true;
}

/// `nrm2`'s algorithm at 64 bits, for the portable and `@Vector` kernels:
/// Blue's, as Reference-LAPACK 3.10's `BLAS/SRC/dnrm2.f90` has it
/// (Anderson, "Algorithm 978: Safe Scaling in the Level 1 BLAS", ACM TOMS
/// 44(1), 2017): one pass summing the squares of small, medium and big
/// magnitudes apart, each scaled by a power of two so none overflows or
/// underflows; no division. Infinite for an infinity, NaN for a NaN.
pub fn blue(comptime F: type, comptime w: comptime_int, n: isize, x: [*]F, inc_x: isize) F {
    // Fortran's `digits`, `minexponent` and `maxexponent`
    const digits = std.math.floatMantissaBits(F) + 1;
    const minexp = std.math.floatExponentMin(F) + 1;
    const maxexp = std.math.floatExponentMax(F) + 1;
    const tsml = std.math.ldexp(@as(F, 1), @divFloor(minexp, 2));
    const tbig = std.math.ldexp(@as(F, 1), @divFloor(maxexp - digits + 1, 2));
    const ssml = std.math.ldexp(@as(F, 1), -@divFloor(minexp - digits, 2));
    const sbig = std.math.ldexp(@as(F, 1), -@divFloor(maxexp + digits, 2));
    var notbig = true;
    var asml: F = 0;
    var amed: F = 0;
    var abig: F = 0;
    var i: isize = 0;
    while (i < n) : (i += 1) {
        inline for (0..w) |c| {
            const ax = @abs(at(F, x, i * inc_x * w + c).*);
            if (ax > tbig) {
                abig += (ax * sbig) * (ax * sbig);
                notbig = false;
            } else if (ax < tsml) {
                if (notbig) asml += (ax * ssml) * (ax * ssml);
            } else {
                amed += ax * ax;
            }
        }
    }
    // a NaN lands in `amed`
    const mixed = amed > 0 or amed != amed;
    if (abig > 0) {
        if (mixed) abig += (amed * sbig) * sbig;
        return @sqrt(abig) / sbig;
    }
    if (asml > 0) {
        if (!mixed) return @sqrt(asml) / ssml;
        const m = @sqrt(amed);
        const s = @sqrt(asml) / ssml;
        const lo = if (s > m) m else s;
        const hi = if (s > m) s else m;
        return @sqrt(hi * hi * (1 + (lo / hi) * (lo / hi)));
    }
    return @sqrt(amed);
}
