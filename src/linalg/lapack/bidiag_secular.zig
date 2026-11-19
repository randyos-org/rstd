//! The secular equation of `bidiag_svd.zig`'s divide and conquer, from faer v0.24.4
//! (`faer/src/linalg/svd/bidiag_svd.rs`, MIT).

const std = @import("std");
const mem = @import("../../mem.zig");
const blas_calls = @import("blas_calls.zig");

const Io = std.Io;

/// The secular function's data: the deflated `col0` and `diag` entries.
pub fn Secular(comptime F: type) type {
    return struct {
        col0_perm: []const F,
        diag_perm: []const F,

        const Self = @This();

        /// `secular_eq`: `1 + Σ c² / ((d - shift - mu) (d + shift + mu))`.
        fn eval(e: Self, shift: F, mu: F) F {
            var res: F = 1;
            for (e.col0_perm, e.diag_perm) |c, d| res += (c / ((d - shift) - mu)) * (c / ((d + shift) + mu));
            return res;
        }

        /// `batch_secular_eq`: four at once, its own way.
        fn batch(e: Self, shift: [4]F, mu: [4]F) [4]F {
            var res: [4]F = @splat(1);
            for (e.col0_perm, e.diag_perm) |c, d| {
                for (&res, mu, shift) |*r, m, s| r.* += (c * c) / (((d - s) - m) * ((d + s) + m));
            }
            return res;
        }
    };
}

/// `compute_singular_values`: each root as `shift + mu`, `shift` the nearer
/// pole, for accuracy.
pub fn singularValues(io: Io, comptime F: type, shifts: []F, mus: []F, s: []F, diag: []const F, col0: []const F, sec: Secular(F)) void {
    const n = diag.len;
    var actual_n = n;
    while (actual_n > 1 and col0[actual_n - 1] == 0) actual_n -= 1;
    for (0..n) |k| {
        s[k] = 0;
        shifts[k] = 0;
        mus[k] = 0;
        if (col0[k] == 0 or actual_n == 1) {
            s[k] = if (k == 0) col0[0] else diag[k];
            shifts[k] = s[k];
            mus[k] = 0;
            continue;
        }
        const last_k = k == actual_n - 1;
        const left = diag[k];
        const right = if (last_k) diag[actual_n - 1] + blas_calls.nrm2(io, F, n, @constCast(col0.ptr), 1) else blk: {
            var l = k + 1;
            while (col0[l] == 0) l += 1;
            break :blk diag[l];
        };
        const r = rootFinder(F, sec, left, right, last_k);
        s[k] = r[0] + r[1];
        shifts[k] = r[0];
        mus[k] = r[1];
    }
}

/// `secular_eq_root_finder`: the root in `(left, right)` of the increasing
/// secular function `sec` (its `eval(shift, mu)` and four at once,
/// `batch`), as `.{ shift, mu }`: bracketing and bisection at geometric
/// midpoints, then the secant method, then bisection if it leaves the
/// bracket. The tridiagonal eigensolver shares it (`tridiag_evd.zig`).
pub fn rootFinder(comptime F: type, sec: anytype, left: F, right: F, last: bool) [2]F {
    const two: F = 2;
    const half: F = 0.5;
    const epsilon = std.math.floatEps(F);
    const inf = std.math.inf(F);
    const mid = left + (right - left) * half;
    const f0 = sec.batch(.{ 0, left, left, right }, .{ mid, if (last) right - left else (right - left) * half, (right - left) * half, (left - right) * half });
    var f_mid = f0[0];
    const f_max = f0[1];
    const f_mid_left_shift = f0[2];
    const f_mid_right_shift = f0[3];
    var shift: F, const mu0: F = if (last or f_mid > 0) .{ left, (right - left) * half } else .{ right, (left - right) * half };
    if (f_mid_left_shift <= 0 and f_mid_right_shift > 0) return .{ shift, mu0 };
    if (!last) {
        if (shift == left) {
            if (f_mid_left_shift < 0) {
                shift = right;
                f_mid = f_mid_right_shift;
            }
        } else if (f_mid_right_shift > 0) {
            shift = left;
            f_mid = f_mid_left_shift;
        }
    }
    var left_shifted: F = undefined;
    var f_left: F = undefined;
    var right_shifted: F = undefined;
    var f_right: F = undefined;
    if (shift == left) {
        left_shifted = 0;
        f_left = -inf;
        right_shifted = if (last) right - left else (right - left) * half;
        f_right = if (last) f_max else f_mid;
    } else {
        left_shifted = (left - right) * half;
        f_left = f_mid;
        right_shifted = 0;
        f_right = inf;
    }
    var iteration_count: usize = 0;
    var f_prev = f_mid;
    const half0 = half;
    const half1 = half0 * half0;
    const half2 = half1 * half1;
    const half3 = half2 * half2;
    const base = if (shift == left) right_shifted else left_shifted;
    const mu_values: [4]F = .{ base * half3, base * half2, base * half1, base * half0 };
    const f_values = sec.batch(@splat(shift), mu_values);
    if (shift == left) {
        var i: usize = 0;
        for (mu_values, f_values, 0..) |m, f, idx| if (f < 0) {
            left_shifted = m;
            f_left = f;
            i = idx + 1;
        };
        if (i < 4) {
            right_shifted = mu_values[i];
            f_right = f_values[i];
        }
    } else {
        var i: usize = 0;
        for (mu_values, f_values, 0..) |m, f, idx| if (f > 0) {
            right_shifted = m;
            f_right = f;
            i = idx + 1;
        };
        if (i < 4) {
            left_shifted = mu_values[i];
            f_left = f_values[i];
        }
    }
    while (right_shifted - left_shifted > two * epsilon * @max(@abs(left_shifted), @abs(right_shifted))) {
        const mid_arith = (left_shifted + right_shifted) * half;
        var mid_geo = @sqrt(@abs(left_shifted)) * @sqrt(@abs(right_shifted));
        if (left_shifted < 0) mid_geo = -mid_geo;
        const mid_shifted = if (mid_geo == 0) mid_arith else mid_geo;
        const fm = sec.eval(shift, mid_shifted);
        if (fm == 0) {
            return .{ shift, mid_shifted };
        } else if (fm > 0) {
            right_shifted = mid_shifted;
            f_prev = f_right;
            f_right = fm;
        } else {
            left_shifted = mid_shifted;
            f_prev = f_left;
            f_left = fm;
        }
        if (iteration_count == 4) break;
        iteration_count += 1;
    }
    const args: [4]F = if (left_shifted == 0)
        .{ right_shifted * two, right_shifted, f_prev, f_right }
    else if (right_shifted == 0)
        .{ left_shifted * two, left_shifted, f_prev, f_left }
    else
        .{ left_shifted, right_shifted, f_left, f_right };
    const sc = secant(F, sec, shift, left, right, args[0], args[1], args[2], args[3]);
    var mu_cur = sc.mu_cur;
    if (sc.left) |lc| if (sc.right) |rc| if (lc < rc) {
        if (lc > left_shifted) left_shifted = lc;
        if (rc < right_shifted) right_shifted = rc;
    };
    if (sc.use_bisection) {
        while (right_shifted - left_shifted > two * epsilon * @max(@abs(left_shifted), @abs(right_shifted))) {
            const m = (left_shifted + right_shifted) * half;
            const fm = sec.eval(shift, m);
            if (fm == 0) {
                break;
            } else if (fm > 0) {
                right_shifted = m;
            } else {
                left_shifted = m;
            }
        }
        mu_cur = (left_shifted + right_shifted) * half;
    }
    return .{ shift, mu_cur };
}

fn SecantResult(comptime F: type) type {
    return struct { use_bisection: bool, mu_cur: F, left: ?F, right: ?F };
}

/// `rootFinder`'s `secant` closure.
fn secant(comptime F: type, sec: anytype, shift: F, left: F, right: F, mu_cur_in: F, mu_prev_in: F, f_cur_in: F, f_prev_in: F) SecantResult(F) {
    const epsilon = std.math.floatEps(F);
    var mu_cur = mu_cur_in;
    var mu_prev = mu_prev_in;
    var f_cur = f_cur_in;
    var f_prev = f_prev_in;
    if (@abs(f_prev) < @abs(f_cur)) {
        mem.swap(F, &f_prev, &f_cur);
        mem.swap(F, &mu_prev, &mu_cur);
    }
    var left_candidate: ?F = null;
    var right_candidate: ?F = null;
    var use_bisection = false;
    const same_sign = (f_prev > 0) == (f_cur > 0);
    if (!same_sign) {
        if (mu_cur < mu_prev) {
            left_candidate = mu_cur;
            right_candidate = mu_prev;
        } else {
            left_candidate = mu_prev;
            right_candidate = mu_cur;
        }
    }
    while (f_cur != 0 and @abs(mu_cur - mu_prev) > 8 * epsilon * @max(@abs(mu_cur), @abs(mu_prev)) and @abs(f_cur - f_prev) > epsilon and !use_bisection) {
        const a = (f_cur - f_prev) * (mu_prev * mu_cur) / (mu_prev - mu_cur);
        const b = f_cur - a / mu_cur;
        const mu_zero = -a / b;
        const f_zero = sec.eval(shift, mu_zero);
        if (f_zero < 0) {
            left_candidate = mu_zero;
        } else {
            right_candidate = mu_zero;
        }
        mu_prev = mu_cur;
        f_prev = f_cur;
        mu_cur = mu_zero;
        f_cur = f_zero;
        if (shift == left and (mu_cur < 0 or mu_cur > right - left)) use_bisection = true;
        if (shift == right and (mu_cur > 0 or mu_cur < left - right)) use_bisection = true;
        if (@abs(f_cur) > @abs(f_prev)) {
            var k: F = 1;
            for (0..4) |_| {
                const mu_opposite = -a / (k * f_zero + b);
                const f_opposite = sec.eval(shift, mu_opposite);
                if (f_zero < 0 and f_opposite >= 0) {
                    right_candidate = mu_opposite;
                    break;
                }
                if (f_zero > 0 and f_opposite <= 0) {
                    left_candidate = mu_opposite;
                    break;
                }
                k *= 2;
            }
            use_bisection = true;
        }
    }
    // Faer stops here with whatever the secant reached when `f` stops changing
    // by more than `eps`, which happens far from a root that lies within
    // rounding of a pole (a tiny `col0` entry): then `mu` is wrong by orders
    // of magnitude, and `perturbCol0` with it. Bisect instead (issue #250).
    if (!use_bisection and f_cur != 0 and !(@abs(f_cur - f_prev) > epsilon) and
        @abs(mu_cur - mu_prev) > 8 * epsilon * @max(@abs(mu_cur), @abs(mu_prev))) use_bisection = true;
    return .{ .use_bisection = use_bisection, .mu_cur = mu_cur, .left = left_candidate, .right = right_candidate };
}
