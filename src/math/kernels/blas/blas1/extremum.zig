//! The largest or smallest element, or its 1-based index: the portable
//! kernels of `amax`, `amin`, `iamax`, `iamin`, `imax`, `imin`, `max` and
//! `min`, which every core's `@Vector` kernels fall back to (`vector.zig`).
//!
//! The first NaN wins, as in Reference-LAPACK master's `i?amax` (`50dc1bf`,
//! `b358788`), BLIS's `amaxv` and arXiv:2207.09281: its index, and NaN for
//! the value. Else the scan replaces element 0 only on a strict `>` (`<`),
//! so ties keep the first. With `abs` an element counts by `|x|`, by
//! `|re| + |im|` for complex (`CABS1`, not the modulus); where that overflows
//! for a finite element, the index compares `|re|/4 + |im|/4`, as
//! `izamax.f90` does, so an infinite part still beats it. `n <= 0` or
//! `inc_x <= 0` gives 0. Imports stay inside `linalg/blas/`.

const std = @import("std");
const types = @import("../types.zig");
const util = @import("../util.zig");
const Vtable = @import("../Vtable.zig");

const BlasLong = types.BlasLong;
const XDouble = types.XDouble;

pub const Find = struct {
    /// Compare `|x|` rather than `x`.
    abs: bool,
    /// The largest rather than the smallest.
    max: bool,
};

/// The kernels over elements of `w` floats (2 for complex).
pub fn Extremum(comptime F: type, comptime w: comptime_int, comptime find: Find) type {
    return struct {
        fn key(x: [*]F, at: BlasLong) F {
            if (w == 2) return util.abs(util.at(F, x, at).*) + util.abs(util.at(F, x, at + 1).*);
            const v = util.at(F, x, at).*;
            return if (find.abs) @abs(v) else v;
        }

        fn better(k: F, best: F) bool {
            return if (find.max) k > best else k < best;
        }

        /// The extreme element or the first NaN: its index from 0 and key;
        /// `n >= 1`, `inc_x >= 1`.
        pub fn scan(n: BlasLong, x: [*]F, inc_x: BlasLong) struct { BlasLong, F } {
            var best: BlasLong = 0;
            var best_key = key(x, 0);
            if (best_key != best_key) return .{ 0, best_key };
            var i: BlasLong = 1;
            while (i < n) : (i += 1) {
                const k = key(x, i * inc_x * w);
                if (k != k) return .{ i, k };
                if (better(k, best_key)) {
                    best = i;
                    best_key = k;
                }
            }
            if (w == 2 and std.math.isInf(best_key)) return .{ quarter(n, x, inc_x), best_key };
            return .{ best, best_key };
        }

        /// The extreme by `|re|/4 + |im|/4`, finite for finite parts.
        fn quarter(n: BlasLong, x: [*]F, inc_x: BlasLong) BlasLong {
            var best: BlasLong = 0;
            var best_key: F = 0;
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                const at = i * inc_x * w;
                const k = util.abs(util.at(F, x, at).*) / 4 + util.abs(util.at(F, x, at + 1).*) / 4;
                if (i == 0 or better(k, best_key)) {
                    best = i;
                    best_key = k;
                }
            }
            return best;
        }

        pub fn value(n: BlasLong, x: [*]F, inc_x: BlasLong) callconv(types.kernelCc(F)) F {
            if (n <= 0 or inc_x <= 0) return 0;
            return scan(n, x, inc_x)[1];
        }

        pub fn index(n: BlasLong, x: [*]F, inc_x: BlasLong) callconv(types.kernelCc(F)) BlasLong {
            if (n <= 0 or inc_x <= 0) return 0;
            return scan(n, x, inc_x)[0] + 1;
        }
    };
}

/// The portable kernels over `T`, any float or its `Complex`, by the slot
/// they would fill (`portable.kernel`); `max` and `min` are real only.
pub fn Portable(comptime T: type) type {
    const F = types.Float(T);
    const w = if (types.isComplex(T)) 2 else 1;
    return struct {
        pub const amax_k = &Extremum(F, w, amax).value;
        pub const amin_k = &Extremum(F, w, amin).value;
        pub const max_k = &Extremum(F, w, max).value;
        pub const min_k = &Extremum(F, w, min).value;
        pub const iamax_k = &Extremum(F, w, amax).index;
        pub const iamin_k = &Extremum(F, w, amin).index;
        pub const imax_k = &Extremum(F, w, max).index;
        pub const imin_k = &Extremum(F, w, min).index;
    };
}

pub const amax: Find = .{ .abs = true, .max = true };
pub const amin: Find = .{ .abs = true, .max = false };
pub const max: Find = .{ .abs = false, .max = true };
pub const min: Find = .{ .abs = false, .max = false };

/// The precisions of `routines`' slots, by position: `s`, `d`, `q`, `c`, `z`,
/// `x`; float and floats an element.
pub const precisions = .{ .{ f32, 1 }, .{ f64, 1 }, .{ XDouble, 1 }, .{ f32, 2 }, .{ f64, 2 }, .{ XDouble, 2 } };

/// Each routine's search, whether it gives the index, and its slots but
/// bfloat16's; `max` and `min` are real only.
pub const routines = .{
    .{ amax, false, .{ .samax_k, .damax_k, .qamax_k, .camax_k, .zamax_k, .xamax_k } },
    .{ amin, false, .{ .samin_k, .damin_k, .qamin_k, .camin_k, .zamin_k, .xamin_k } },
    .{ amax, true, .{ .isamax_k, .idamax_k, .iqamax_k, .icamax_k, .izamax_k, .ixamax_k } },
    .{ amin, true, .{ .isamin_k, .idamin_k, .iqamin_k, .icamin_k, .izamin_k, .ixamin_k } },
    .{ max, true, .{ .ismax_k, .idmax_k, .iqmax_k } },
    .{ min, true, .{ .ismin_k, .idmin_k, .iqmin_k } },
    .{ max, false, .{ .smax_k, .dmax_k, .qmax_k } },
    .{ min, false, .{ .smin_k, .dmin_k, .qmin_k } },
};

/// Points every extremum slot but bfloat16's at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    inline for (routines) |r| inline for (r[2], 0..) |slot, i| {
        const E = Extremum(precisions[i][0], precisions[i][1], r[0]);
        @field(vt, @tagName(slot)) = if (r[1]) &E.index else &E.value;
    };
}
