//! `xTREXC`: moving one block of a Schur form, Reference-LAPACK 3.12.0 (`lapack-netlib/SRC/dtrexc.f`,
//! `ztrexc.f`'s rotations, OpenBLAS v0.3.34-317; modified BSD). Faer's `schur_move` translates the same `dtrexc.f`.

const std = @import("std");
const blas = @import("../blas/__root__.zig");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");
const schur = @import("schur.zig");
const laexc_swap = @import("laexc.zig");

const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const Error = ops.Error;
const Compq = schur.Compq;
const laexc = laexc_swap.laexc;

/// What `trexc` returns: a slice too short or a leading dimension too small
/// (`q` checked only with `compq` `v`), `RowOutOfRange` for an `ifst` or
/// `ilst` past `T`'s rows, `SwapRejected` for the reference's `INFO = 1`
/// (real only): two adjacent blocks were too close to swap stably, `T`
/// part reordered and `ilst` the moving block's current row.
pub const TrexcError = interface.ShapeError || error{ RowOutOfRange, SwapRejected } || Io.Cancelable;

/// `?TREXC` (`dtrexc.f`, `ztrexc.f`): moves the diagonal block of the
/// `n x n` Schur form `T` at row `ifst` (0-based) to row `ilst` by
/// orthogonal (unitary) similarity, `T := Zᴴ T Z`, and with `compq` `v`
/// `Q := Q Z`. Real: `T` is quasi-triangular with standardized 2 x 2
/// blocks; `ifst` moves to its block's first row and `ilst` to the moved
/// block's final first row, as the reference's. Complex: both unchanged.
pub fn trexc(io: Io, comptime T: type, compq: Compq, n: usize, t: []T, ldt: usize, q: []T, ldq: usize, ifst: *usize, ilst: *usize) TrexcError!void {
    comptime ops.requireFloat(T, "trexc");
    const wantq = compq == .v;
    try util.holds(t.len, n, n, ldt);
    if (wantq) try util.holds(q.len, n, n, ldq);
    if (n > 0 and (ifst.* >= n or ilst.* >= n)) return error.RowOutOfRange;
    if (n <= 1) return;
    const qq: ?[*]T = if (wantq) q.ptr else null;
    if (comptime isComplex(T)) {
        moveComplex(T, n, t.ptr, ldt, qq, ldq, ifst.* + 1, ilst.* + 1);
        return;
    }
    // the reference's 1-based rows
    var f = ifst.* + 1;
    var l = ilst.* + 1;
    defer {
        ifst.* = f - 1;
        ilst.* = l - 1;
    }
    const done = moveReal(io, T, n, t.ptr, ldt, qq, ldq, &f, &l) catch |e| switch (e) {
        error.OutOfMemory => unreachable,
        error.Canceled => return error.Canceled,
    };
    if (!done) return error.SwapRejected;
}

/// `ztrexc.f` after its checks (1-based `ifst`, `ilst`).
pub fn moveComplex(comptime T: type, n: usize, t: [*]T, ldt: usize, q: ?[*]T, ldq: usize, ifst: usize, ilst: usize) void {
    if (ifst == ilst) return;
    // 0-based swaps of k and k + 1, in the reference's order
    var k: usize = if (ifst < ilst) ifst - 1 else ifst - 2;
    while (true) {
        const t11 = t[k + k * ldt];
        const t22 = t[k + 1 + (k + 1) * ldt];
        var r = t[k + (k + 1) * ldt];
        var cs: Float(T) = undefined;
        var sn: T = undefined;
        blas.zrotg(T, &r, ops.sub(T, t22, t11), &cs, &sn);
        if (k + 2 < n) blas_calls.zrot(T, n - k - 2, t + k + (k + 2) * ldt, ldt, t + k + 1 + (k + 2) * ldt, ldt, cs, sn);
        blas_calls.zrot(T, k, t + k * ldt, 1, t + (k + 1) * ldt, 1, cs, ops.conj(T, sn));
        t[k + k * ldt] = t22;
        t[k + 1 + (k + 1) * ldt] = t11;
        if (q) |qq| blas_calls.zrot(T, n, qq + k * ldq, 1, qq + (k + 1) * ldq, 1, cs, ops.conj(T, sn));
        if (ifst < ilst) {
            if (k + 2 == ilst) break;
            k += 1;
        } else {
            if (k + 1 == ilst) break;
            k -= 1;
        }
    }
}

/// `dtrexc.f` after its checks: `ifst` and `ilst` (1-based) as the
/// reference leaves them; false for `INFO = 1`.
pub fn moveReal(io: Io, comptime F: type, n: usize, t: [*]F, ldt: usize, q: ?[*]F, ldq: usize, ifst: *usize, ilst: *usize) Error!bool {
    // `T(i, j)`, 1-based
    const at = struct {
        fn f(p: [*]F, ld: usize, i: usize, j: usize) F {
            return p[(i - 1) + (j - 1) * ld];
        }
    }.f;
    if (ifst.* > 1 and at(t, ldt, ifst.*, ifst.* - 1) != 0) ifst.* -= 1;
    var nbf: usize = 1;
    if (ifst.* < n and at(t, ldt, ifst.* + 1, ifst.*) != 0) nbf = 2;
    if (ilst.* > 1 and at(t, ldt, ilst.*, ilst.* - 1) != 0) ilst.* -= 1;
    var nbl: usize = 1;
    if (ilst.* < n and at(t, ldt, ilst.* + 1, ilst.*) != 0) nbl = 2;
    if (ifst.* == ilst.*) return true;
    const ex = struct {
        fn f(i: Io, nn: usize, tt: [*]F, ld: usize, qq: ?[*]F, lq: usize, j1: usize, n1: usize, n2: usize) Error!bool {
            return laexc(i, F, nn, tt, ld, qq, lq, j1 - 1, n1, n2);
        }
    }.f;
    var here = ifst.*;
    if (ifst.* < ilst.*) {
        // forward
        if (nbf == 2 and nbl == 1) ilst.* -= 1;
        if (nbf == 1 and nbl == 2) ilst.* += 1;
        while (true) {
            if (nbf == 1 or nbf == 2) {
                // the block is whole
                var nbnext: usize = 1;
                if (here + nbf + 1 <= n and at(t, ldt, here + nbf + 1, here + nbf) != 0) nbnext = 2;
                if (!try ex(io, n, t, ldt, q, ldq, here, nbf, nbnext)) {
                    ilst.* = here;
                    return false;
                }
                here += nbnext;
                if (nbf == 2 and at(t, ldt, here + 1, here) == 0) nbf = 3;
            } else {
                // a 2 x 2 block split in two 1 x 1s: move each
                var nbnext: usize = 1;
                if (here + 3 <= n and at(t, ldt, here + 3, here + 2) != 0) nbnext = 2;
                if (!try ex(io, n, t, ldt, q, ldq, here + 1, 1, nbnext)) {
                    ilst.* = here;
                    return false;
                }
                if (nbnext == 1) {
                    _ = try ex(io, n, t, ldt, q, ldq, here, 1, nbnext);
                    here += 1;
                } else {
                    if (at(t, ldt, here + 2, here + 1) == 0) nbnext = 1;
                    if (nbnext == 2) {
                        if (!try ex(io, n, t, ldt, q, ldq, here, 1, nbnext)) {
                            ilst.* = here;
                            return false;
                        }
                        here += 2;
                    } else {
                        _ = try ex(io, n, t, ldt, q, ldq, here, 1, 1);
                        _ = try ex(io, n, t, ldt, q, ldq, here + 1, 1, 1);
                        here += 2;
                    }
                }
            }
            if (!(here < ilst.*)) break;
        }
    } else {
        // backward
        while (true) {
            if (nbf == 1 or nbf == 2) {
                var nbnext: usize = 1;
                if (here >= 3 and at(t, ldt, here - 1, here - 2) != 0) nbnext = 2;
                if (!try ex(io, n, t, ldt, q, ldq, here - nbnext, nbnext, nbf)) {
                    ilst.* = here;
                    return false;
                }
                here -= nbnext;
                if (nbf == 2 and at(t, ldt, here + 1, here) == 0) nbf = 3;
            } else {
                var nbnext: usize = 1;
                if (here >= 3 and at(t, ldt, here - 1, here - 2) != 0) nbnext = 2;
                if (!try ex(io, n, t, ldt, q, ldq, here - nbnext, nbnext, 1)) {
                    ilst.* = here;
                    return false;
                }
                if (nbnext == 1) {
                    _ = try ex(io, n, t, ldt, q, ldq, here, nbnext, 1);
                    here -= 1;
                } else {
                    if (at(t, ldt, here, here - 1) == 0) nbnext = 1;
                    if (nbnext == 2) {
                        if (!try ex(io, n, t, ldt, q, ldq, here - 1, 2, 1)) {
                            ilst.* = here;
                            return false;
                        }
                        here -= 2;
                    } else {
                        _ = try ex(io, n, t, ldt, q, ldq, here, 1, 1);
                        _ = try ex(io, n, t, ldt, q, ldq, here - 1, 1, 1);
                        here -= 2;
                    }
                }
            }
            if (!(here > ilst.*)) break;
        }
    }
    ilst.* = here;
    return true;
}
