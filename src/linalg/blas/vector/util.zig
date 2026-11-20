//! What the `@Vector` kernels share (`real.zig`, `complex.zig`, `half.zig`):
//! the lane count, the fused multiply-add, `axpby`'s operands, whether
//! `@max` is an instruction and `nrm2`'s root.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../types.zig");
const util = @import("../util.zig");
const fma = @import("../fma.zig");

const BlasLong = types.BlasLong;

/// Lanes of one native vector of `F`: a NEON register on AArch64, where
/// `std.simd` counts SVE as 256 bits, which LLVM splits into NEON pairs.
pub fn lanes(comptime F: type) comptime_int {
    if (comptime sysinfo.target.cpu.arch.isAarch64()) return 16 / @sizeOf(F);
    return std.simd.suggestVectorLength(F) orelse 1;
}

/// `e` as `T`: itself where `T` is a float, else in every lane of `T`.
pub inline fn splatTo(comptime T: type, e: anytype) T {
    return if (@typeInfo(T) == .vector) @splat(e) else e;
}

pub inline fn at(comptime F: type, p: [*]F, i: BlasLong) *F {
    return util.at(F, p, i);
}

pub const abs = util.abs;

/// `a*b + c`, fused wherever an instruction does it (`fma.fusible`), whatever
/// the setting; elsewhere rounded twice, as `@mulAdd` there is a software
/// `fma` call per element.
pub inline fn mulAdd(comptime T: type, a: T, b: T, c: T) T {
    return if (comptime fma.fusible(T)) fma.fusedMadd(T, a, b, c) else a * b + c;
}

/// Which operands `axpby` reads: a zero scalar's is not.
pub const Pick = enum {
    zero,
    ax,
    by,
    both,

    pub fn of(a: bool, b: bool) Pick {
        return if (a) (if (b) .both else .ax) else if (b) .by else .zero;
    }

    /// `alpha*x + beta*y` from `x` and `y` at `xp` and `yp` as `T`s, a
    /// vector or one float; `beta*y` first, `alpha*x` fused in.
    pub inline fn apply(comptime p: Pick, comptime T: type, a: anytype, b: @TypeOf(a), xp: anytype, yp: anytype) T {
        const F = @TypeOf(a);
        const s = struct {
            inline fn get(q: anytype) T {
                return if (T == F) q.* else q[0..@typeInfo(T).vector.len].*;
            }
        };
        return switch (p) {
            .zero => splatTo(T, @as(F, 0)),
            .ax => splatTo(T, a) * s.get(xp),
            .by => splatTo(T, b) * s.get(yp),
            .both => mulAdd(T, splatTo(T, a), s.get(xp), splatTo(T, b) * s.get(yp)),
        };
    }
};

/// Whether `@max` and `@min` of `F` are instructions here: elsewhere each
/// is a call to libm's `fmax` or `fmin`, lane by lane (ARMv7 VFP, SPARC,
/// pre-VSX POWER, MIPS before release 6, MSA too, the x87).
pub fn hasMinMax(comptime F: type) bool {
    const cpu = sysinfo.target.cpu;
    return switch (cpu.arch) {
        .x86, .x86_64 => !util.onX87(F),
        .arm, .armeb, .thumb, .thumbeb => cpu.has(.arm, .fp_armv8),
        .powerpc, .powerpcle, .powerpc64, .powerpc64le => cpu.has(.powerpc, .vsx),
        .mips, .mipsel, .mips64, .mips64el => cpu.hasAny(.mips, &.{ .mips32r6, .mips64r6 }),
        .sparc, .sparc64 => false,
        else => true,
    };
}

/// `nrm2` from `s`, the plain sum of squares of `n` elements of `w` floats.
/// Summed in `f64` for `f32`, it is exact enough as it is. Else `@sqrt(s)`
/// where it neither overflowed nor could have lost more than `eps^2` of
/// itself to underflow, and Blue's algorithm (`util.blue`) where it could have.
pub fn root(comptime F: type, comptime w: comptime_int, s: anytype, n: BlasLong, x: [*]F, incx: BlasLong) F {
    if (comptime @TypeOf(s) != F) return @floatCast(@sqrt(s));
    const lo = std.math.floatMin(F) / std.math.floatEps(F) * @as(F, @floatFromInt(w * n));
    if (s >= lo and s <= std.math.floatMax(F)) return @sqrt(s);
    return util.blue(F, w, n, x, incx);
}
