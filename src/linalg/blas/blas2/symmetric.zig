//! The portable symv and hemv kernels and their tests' reference, shared by
//! `symv/` and `hemv/`: `y += alpha*S*x` for an `m` by `m` matrix `S` given by
//! one stored triangle of `A` (column-major).

const std = @import("std");
const types = @import("../types.zig");
const util = @import("../util.zig");
const check = @import("../blas1/check.zig");

const BlasLong = types.BlasLong;

/// How `S` comes from the stored triangle: symmetric (`S = S^T`), Hermitian
/// (`S = S^H`, the diagonal's imaginary part ignored), or the transpose of
/// the Hermitian one (`HEMVREV`: the `V` and `M` kernels).
pub const Kind = enum { symmetric, hermitian, hermitian_rev };

/// `kernel/generic/symv_k.c` and `zhemv_k.c`, for the columns `offset`
/// selects: the last `offset` of them for the upper triangle, the first
/// `offset` for the lower (`offset == m` is the whole product; the threaded
/// drivers split the columns).
pub fn Kernel(comptime T: type, comptime lower: bool, comptime kind: Kind) type {
    const complex = comptime check.isComplex(T);
    const F = check.FloatOf(T);
    return struct {
        fn get(a: [*]F, i: BlasLong) T {
            if (complex) return .init(util.at(F, a, 2 * i).*, util.at(F, a, 2 * i + 1).*);
            return util.at(F, a, i).*;
        }
        fn add(y: [*]F, i: BlasLong, v: T) void {
            if (complex) {
                util.at(F, y, 2 * i).* += v.re;
                util.at(F, y, 2 * i + 1).* += v.im;
            } else util.at(F, y, i).* += v;
        }
        fn mul(p: T, q: T) T {
            return if (complex) p.mul(q) else p * q;
        }
        /// The stored `a` as row `i`'s entry, and as its mirror's.
        fn pair(a: T) [2]T {
            return switch (kind) {
                .symmetric => .{ a, a },
                .hermitian => .{ a, a.conjugate() },
                .hermitian_rev => .{ a.conjugate(), a },
            };
        }
        fn diag(a: T) T {
            return if (kind == .symmetric) a else .init(a.re, 0);
        }
        fn run(m: BlasLong, offset: BlasLong, alpha: T, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong) void {
            var j: BlasLong = if (lower) 0 else m - offset;
            const end: BlasLong = if (lower) offset else m;
            while (j < end) : (j += 1) {
                const xj = mul(alpha, get(x, j * inc_x));
                var acc: T = if (complex) .init(0, 0) else 0;
                var i: BlasLong = if (lower) j + 1 else 0;
                const stop: BlasLong = if (lower) m else j;
                while (i < stop) : (i += 1) {
                    const e = pair(get(a, i + j * lda));
                    add(y, i * inc_y, mul(e[0], xj));
                    acc = if (complex) acc.add(mul(e[1], get(x, i * inc_x))) else acc + e[1] * get(x, i * inc_x);
                }
                add(y, j * inc_y, mul(alpha, acc));
                add(y, j * inc_y, mul(diag(get(a, j + j * lda)), xj));
            }
        }
        pub fn real(m: BlasLong, offset: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(.c) i32 {
            run(m, offset, alpha, a, lda, x, inc_x, y, inc_y);
            return 0;
        }
        pub fn cplx(m: BlasLong, offset: BlasLong, alpha_r: F, alpha_i: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(.c) i32 {
            run(m, offset, .init(alpha_r, alpha_i), a, lda, x, inc_x, y, inc_y);
            return 0;
        }
    };
}

/// Orders across the kernels' blocks; small enough for the rigs' stacks.
const sizes = [_]usize{ 1, 2, 3, 4, 7, 8, 9, 16, 17 };
const strides = [_][2]usize{ .{ 1, 1 }, .{ 2, 3 } };

/// `kern` against plain arithmetic over the full matrix, for the whole
/// product and for part of it, split as the threaded driver
/// (`driver/level2/symv_thread.c`) splits it: on a multiple of 4 columns,
/// which the assembly kernels rely on.
pub fn checkKernel(comptime T: type, comptime lower: bool, comptime kind: Kind, kern: anytype, name: []const u8) !void {
    const F = check.FloatOf(T);
    const complex = comptime check.isComplex(T);
    const alpha: T = if (complex) .init(1.25, -0.5) else 1.25;
    for (sizes) |m| for (strides) |inc| for ([_]usize{ m, if (m > 4) (if (lower) 4 else m - 4) else m }) |offset| {
        const lda = m + 1;
        var a: [18 * 17]T = undefined;
        var x: [17 * 3]T = undefined;
        var y: [17 * 3]T = undefined;
        check.fill(T, &a, m * 5 + offset);
        check.fill(T, &x, m + 3);
        check.fill(T, &y, m + 7);
        var want = y;
        const lo = if (lower) 0 else m - offset;
        const hi = if (lower) offset else m;
        for (lo..hi) |j| for (0..m) |i| {
            const stored = if (lower) i >= j else i <= j;
            if (!stored) continue;
            const e = a[i + j * lda];
            // `S[i][j]` from the stored entry, and `S[j][i]` if not on the diagonal
            const sij: T, const sji: T = switch (kind) {
                .symmetric => .{ e, e },
                .hermitian => .{ e, if (complex) e.conjugate() else e },
                .hermitian_rev => .{ if (complex) e.conjugate() else e, e },
            };
            if (i == j) {
                const d: T = if (kind == .symmetric) e else if (complex) .init(e.re, 0) else e;
                want[j * inc[1]] = if (complex) want[j * inc[1]].add(alpha.mul(d).mul(x[j * inc[0]])) else want[j * inc[1]] + alpha * d * x[j * inc[0]];
            } else if (complex) {
                want[i * inc[1]] = want[i * inc[1]].add(alpha.mul(sij).mul(x[j * inc[0]]));
                want[j * inc[1]] = want[j * inc[1]].add(alpha.mul(sji).mul(x[i * inc[0]]));
            } else {
                want[i * inc[1]] += alpha * sij * x[j * inc[0]];
                want[j * inc[1]] += alpha * sji * x[i * inc[0]];
            }
        };
        if (complex) {
            _ = kern(@intCast(m), @intCast(offset), alpha.re, alpha.im, @ptrCast(&a), @intCast(lda), @ptrCast(&x), @intCast(inc[0]), @ptrCast(&y), @intCast(inc[1]), check.work(F));
        } else {
            _ = kern(@intCast(m), @intCast(offset), alpha, &a, @intCast(lda), &x, @intCast(inc[0]), &y, @intCast(inc[1]), check.work(F));
        }
        check.expectCloseSlices(T, &want, &y, check.tolerance(F, 8 * m)) catch |e| {
            check.say("{s} m={d} offset={d}: ", .{ name, m, offset });
            return check.case(e, "symv", T, m, inc);
        };
    };
}
