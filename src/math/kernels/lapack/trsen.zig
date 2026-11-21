//! `xTRSEN`: reordering a Schur form by a selection, with the cluster's condition numbers,
//! Reference-LAPACK 3.12.0 (`lapack-netlib/SRC/dtrsen.f`, `ztrsen.f`, OpenBLAS v0.3.34-317; modified BSD).

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const norm = @import("norm.zig");
const trsyl = @import("trsyl.zig");
const schur = @import("schur.zig");
const trexc_move = @import("trexc.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const Error = ops.Error;
const ReorderVectors = schur.ReorderVectors;
const Eigenvalues = schur.Eigenvalues;
const eigenvaluesHold = schur.eigenvaluesHold;
const moveComplex = trexc_move.moveComplex;
const moveReal = trexc_move.moveReal;

/// LAPACK's `JOB`: no condition numbers (`N`), the cluster's (`E`, `S`),
/// the invariant subspace's (`V`, `SEP`), or both (`B`).
pub const ConditionJob = enum { none, eigenvalues, subspace, both };

/// What `trsen` returns: a slice too short (`t`, `select` or `w` under
/// `n`, `q` with `compq` `.update`) or a leading dimension too small,
/// `SwapRejected` for the reference's `INFO = 1` (real only): a swap was
/// rejected, `T` part reordered, the eigenvalues written.
pub const TrsenError = interface.ShapeError || error{SwapRejected} || Allocator.Error || Io.Cancelable;

/// What `trsen` found: the reordered cluster's size (the reference's `m`)
/// and, as `job` asks, the reciprocal condition numbers of the cluster
/// (`eigenvalue_condition`, the reference's `S`) and of the invariant
/// subspace (`subspace_separation`, its `SEP`, an estimate), else null.
pub fn ReorderedSchur(comptime T: type) type {
    return struct { selected: usize, eigenvalue_condition: ?Float(T) = null, subspace_separation: ?Float(T) = null };
}

/// `?TRSEN` (`dtrsen.f`, `ztrsen.f`): reorders the `n x n` Schur form `T`
/// so the eigenvalues `select` picks (`n`; real: either of a pair picks
/// both) lead it, `T := Zᴴ T Z`, `Q := Q Z` with `compq` `.update`. Writes the
/// reordered eigenvalues to `w` and returns their count with the condition
/// numbers `job` asks for.
pub fn trsen(allocator: Allocator, io: Io, comptime T: type, job: ConditionJob, compq: ReorderVectors, select: []const bool, n: usize, t_in: []T, ldt: usize, q: []T, ldq: usize, w: Eigenvalues(T)) TrsenError!ReorderedSchur(T) {
    comptime ops.requireFloat(T, "trsen");
    const wants = job == .eigenvalues or job == .both;
    const wantsp = job == .subspace or job == .both;
    const wantq = compq == .update;
    try util.holds(t_in.len, n, n, ldt);
    if (wantq) try util.holds(q.len, n, n, ldq);
    if (select.len < n) return error.SliceTooShort;
    try eigenvaluesHold(T, w, n);
    const un = n;
    const lt = ldt;
    const t = t_in.ptr;
    const qq: ?[*]T = if (wantq) q.ptr else null;
    // the cluster's size, a real pair counted whole
    var cnt: usize = 0;
    {
        var k: usize = 0;
        while (k < un) : (k += 1) {
            if (!isComplex(T) and k + 1 < un and !ops.isZero(T, t[k + 1 + k * lt])) {
                if (select[k] or select[k + 1]) cnt += 2;
                k += 1;
            } else if (select[k]) cnt += 1;
        }
    }
    var r: ReorderedSchur(T) = .{ .selected = cnt };
    defer eigenvalues(T, un, t, lt, w);
    if (cnt == un or cnt == 0) {
        if (wants) r.eigenvalue_condition = 1;
        if (wantsp) r.subspace_separation = try norm.matNorm(allocator, T, .one, un, un, t, lt);
        return r;
    }
    // move the selected blocks up, in order
    var ks: usize = 0;
    var k: usize = 0;
    while (k < un) : (k += 1) {
        const pair = !isComplex(T) and k + 1 < un and !ops.isZero(T, t[k + 1 + k * lt]);
        const swap = select[k] or (pair and select[k + 1]);
        if (swap) {
            ks += 1;
            if (k + 1 != ks) {
                if (comptime isComplex(T)) {
                    moveComplex(T, un, t, lt, qq, ldq, k + 1, ks);
                } else {
                    var kk = k + 1;
                    if (!try moveReal(io, T, un, t, lt, qq, ldq, &kk, &ks)) return error.SwapRejected;
                }
            }
            if (pair) ks += 1;
        }
        if (pair) k += 1;
    }
    const n1 = cnt;
    const n2 = un - cnt;
    const nn = n1 * n2;
    if (!wants and !wantsp) return r;
    const work = try allocator.alloc(T, 2 * nn);
    defer allocator.free(work);
    const ctx: Sep(T) = .{ .al = allocator, .io = io, .n1 = n1, .n2 = n2, .t = t, .ldt = lt, .x = work.ptr };
    if (wants) {
        // s from the solution of T11 X - X T22 = scale T12
        for (0..n2) |j| @memcpy(work[j * n1 ..][0..n1], t[(n1 + j) * lt ..][0..n1]);
        const scale = try ctx.solve(false, work.ptr);
        const rnorm = try norm.matNorm(allocator, T, .frobenius, n1, n2, work.ptr, n1);
        r.eigenvalue_condition = if (rnorm == 0) 1 else scale / (@sqrt(scale * scale / rnorm + rnorm) * @sqrt(rnorm));
    }
    if (wantsp) {
        var last = Sep(T).Last{ .ctx = ctx };
        const est = try norm.lacn2(allocator, T, Error, nn, work[nn..], work[0..nn], &last, Sep(T).Last.apply);
        r.subspace_separation = last.scale / est;
    }
    return r;
}

/// `trsen`'s Sylvester operator `X ↦ T11 X - X T22` on `n1 x n2` `X`.
fn Sep(comptime T: type) type {
    return struct {
        al: Allocator,
        io: Io,
        n1: usize,
        n2: usize,
        t: [*]T,
        ldt: usize,
        x: [*]T,

        const Self = @This();

        /// Solves for `X` in `x` (the operator's inverse, or its adjoint's),
        /// returning `scale`; a perturbed solve counts, as the reference
        /// ignores its `IERR`.
        fn solve(s: Self, adjoint: bool, x: [*]T) Error!Float(T) {
            const tr: types.Trans = if (adjoint) (if (comptime isComplex(T)) .c else .t) else .n;
            const n1 = s.n1;
            const n2 = s.n2;
            const ld = s.ldt;
            const t11 = s.t[0 .. (n1 - 1) * ld + n1];
            const t22 = (s.t + n1 + n1 * ld)[0 .. (n2 - 1) * ld + n2];
            const r = trsyl.trsyl(s.al, s.io, T, tr, tr, -1, n1, n2, t11, ld, t22, ld, x[0 .. n1 * n2], n1) catch |e| switch (e) {
                error.OutOfMemory => return error.OutOfMemory,
                error.Canceled => return error.Canceled,
                else => unreachable,
            };
            return r.scale;
        }

        /// `lacn2`'s callback, keeping the last solve's `scale`.
        const Last = struct {
            ctx: Self,
            scale: Float(T) = 1,

            fn apply(l: *Last, adjoint: bool, x: [*]T) Error!void {
                l.scale = try l.ctx.solve(adjoint, x);
            }
        };
    };
}

/// `T`'s eigenvalues into `w`: the diagonal, and real pairs `a ± i b` with
/// `b = sqrt(|T(k, k+1)|) sqrt(|T(k+1, k)|)`.
fn eigenvalues(comptime T: type, n: usize, t: [*]const T, ldt: usize, w: Eigenvalues(T)) void {
    if (comptime isComplex(T)) {
        for (0..n) |k| w[k] = t[k + k * ldt];
        return;
    }
    for (0..n) |k| {
        w.re[k] = t[k + k * ldt];
        w.im[k] = 0;
    }
    for (0..n -| 1) |k| if (t[k + 1 + k * ldt] != 0) {
        w.im[k] = @sqrt(@abs(t[k + (k + 1) * ldt])) * @sqrt(@abs(t[k + 1 + k * ldt]));
        w.im[k + 1] = -w.im[k];
    };
}
