const std = @import("std");
const sysinfo = @import("builtin");
const Vtable = @import("Vtable.zig");

/// A vtable slot's name; too many fields for the default branch quota.
pub const Slot = blk: {
    @setEvalBranchQuota(10_000);
    break :blk std.meta.FieldEnum(Vtable);
};

/// The function type of a kernel slot.
pub fn SlotFn(comptime slot: Slot) type {
    return @typeInfo(@FieldType(Vtable, @tagName(slot))).optional.child;
}

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

/// The symbol of `core`'s kernel for `slot` (`blaskern_daxpy_k_haswell`). The
/// per-core objects (`build/blaskern.zig`) export it, with hidden visibility so it
/// never leaves a shared library rstd is built into, and `kernels/<arch>.zig`'s
/// `install` looks it up by the same name.
pub fn kernelSymbol(comptime slot: Slot, comptime core: anytype) []const u8 {
    return "blaskern_" ++ @tagName(slot) ++ "_" ++ @tagName(core);
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
/// their scalar path where this is false (`comparisons.md`, 10). LLVM's
/// auto-vectorizer keeps off NEON's `f32` but not AltiVec's: the AltiVec
/// cores' objects are built without it (`cpu/_corespec.zig`), and the loops
/// built for the target pass their `f32` results through `scalar`.
pub fn vectorIeee(comptime F: type) bool {
    return vectorIeeeOn(F, sysinfo.target.cpu);
}

/// `vectorIeee` for `cpu`, a core's of this target.
pub fn vectorIeeeOn(comptime F: type, comptime cpu: std.Target.Cpu) bool {
    if (F != f32 or sysinfo.target.abi.float() == .soft) return true;
    return switch (cpu.arch) {
        .arm, .armeb, .thumb, .thumbeb => !std.Target.arm.featureSetHas(cpu.features, .neon),
        .powerpc, .powerpcle, .powerpc64, .powerpc64le => !std.Target.powerpc.featureSetHas(cpu.features, .altivec) or
            std.Target.powerpc.featureSetHas(cpu.features, .vsx),
        else => true,
    };
}

/// `x` (a float, or an array of them), through an empty asm in a scalar
/// register where AltiVec's `f32` lanes flush: LLVM's vectorizers cannot
/// widen a loop or a group of operations through it, so a loop that passes
/// what it stores or accumulates through it keeps out of those lanes. `x`
/// itself everywhere else; ARMv7's vectorizers already keep off NEON's `f32`.
pub inline fn scalar(comptime T: type, x: T) T {
    if (comptime vectorIeee(f32) or !sysinfo.target.cpu.arch.isPowerpc()) return x;
    if (@typeInfo(T) == .array) {
        var y = x;
        inline for (&y) |*e| e.* = scalar(@typeInfo(T).array.child, e.*);
        return y;
    }
    if (T != f32) return x;
    return asm (""
        : [y] "=f" (-> f32),
        : [x] "0" (x),
    );
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

/// `b / a` in parts, by Smith's algorithm (the larger part of `a` divided
/// out), as gfortran divides; for the triangular solves of levels 2 and 3.
/// No reciprocal is formed, so a tiny `a` whose reciprocal overflows still
/// gives a finite quotient wherever one exists.
pub inline fn cdiv(comptime F: type, br: F, bi: F, ar: F, ai: F) [2]F {
    if (@abs(ar) >= @abs(ai)) {
        const r = ai / ar;
        const d = ar + ai * r;
        return .{ (br + bi * r) / d, (bi - br * r) / d };
    }
    const r = ar / ai;
    const d = ai + ar * r;
    return .{ (br * r + bi) / d, (bi * r - br) / d };
}

/// What `extendedSum` adds up per float.
pub const Term = enum { abs, square };

/// x86_64 `asum.S`'s and `nrm2.S`'s x87 sums and their complex forms', before
/// `nrm2`'s square root, which every float wider than `f64` takes: eight
/// floats a trip into accumulators by float `% 4`, then each remaining
/// element's float `c` (of `w`) into accumulator `c`, folded
/// `a3 + ((a2 + a0) + a1)`, nothing fused. At `f80` it gives the assembly's
/// bits; `kernel/x86/`'s are in the same order.
pub fn extendedSum(comptime F: type, comptime w: comptime_int, comptime term: Term, n: isize, x: [*]F, inc_x: isize) F {
    const per = 8 / w;
    var a: [4]F = @splat(0);
    var i: isize = 0;
    const main = n & -per;
    while (i < main) : (i += per) {
        inline for (0..8) |f| a[f % 4] += of(F, term, at(F, x, (i + f / w) * inc_x * w + f % w).*);
    }
    while (i < n) : (i += 1) {
        inline for (0..w) |c| a[c] += of(F, term, at(F, x, i * inc_x * w + c).*);
    }
    return a[3] + ((a[2] + a[0]) + a[1]);
}

inline fn of(comptime F: type, comptime term: Term, v: F) F {
    return switch (term) {
        .abs => @abs(v),
        .square => v * v,
    };
}
