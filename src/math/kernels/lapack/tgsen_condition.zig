//! `xTGSEN`'s condition estimates (`IJOB` 1 to 5) for a reordered
//! generalized Schur form, Reference-LAPACK 3.12.0
//! (`lapack-netlib/SRC/dtgsen.f`, `ztgsen.f`, OpenBLAS v0.3.34-317;
//! modified BSD): the projections' norms from a `tgsyl` solve, the
//! separations by `tgsyl`'s Frobenius-norm estimate or `lacn2`'s one-norm.

const std = @import("std");
const types = @import("../blas/types.zig");
const ops = @import("ops.zig");
const norm = @import("norm.zig");
const tgsyl_mod = @import("tgsyl.zig");
const tgexc_mod = @import("tgexc.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const Error = ops.Error;
const Pencil = tgexc_mod.Pencil;

/// LAPACK's `IJOB`: nothing (0), the projections (1), the separations by
/// `tgsyl`'s Frobenius-norm estimate (2) or `lacn2`'s one-norm (3), or the
/// projections with either (4, 5).
pub const Job = enum { none, projections, separations_frobenius, separations_one_norm, both_frobenius, both_one_norm };

/// What `tgsen` found: the reordered cluster's size (the reference's `M`)
/// and, as `job` asks, the reciprocal norms of the projections onto the
/// left and right deflating subspaces (`PL`, `PR`), and the separations
/// `Dif_u`, `Dif_l` of the cluster and the rest (`DIF(1)`, `DIF(2)`,
/// estimates), else null.
pub fn Reordered(comptime F: type) type {
    return struct { selected: usize, projections: ?[2]F = null, separations: ?[2]F = null };
}

fn wantsProjections(job: Job) bool {
    return job == .projections or job == .both_frobenius or job == .both_one_norm;
}

fn wantsFrobenius(job: Job) bool {
    return job == .separations_frobenius or job == .both_frobenius;
}

fn wantsOneNorm(job: Job) bool {
    return job == .separations_one_norm or job == .both_one_norm;
}

/// Nothing moved (`M = 0` or `N`): `PL = PR = 1`, both separations
/// `‖(A, B)‖_F` (`lassq` over each column of `A`, then `B`'s).
pub fn trivial(comptime T: type, p: Pencil(T), job: Job, r: *Reordered(Float(T))) void {
    if (wantsProjections(job)) r.projections = .{ 1, 1 };
    if (!wantsFrobenius(job) and !wantsOneNorm(job)) return;
    var scale: Float(T) = 0;
    var sum: Float(T) = 1;
    for (0..p.n) |i| {
        norm.lassq(T, .of((p.a + i * p.lda)[0..p.n]), &scale, &sum);
        norm.lassq(T, .of((p.b + i * p.ldb)[0..p.n]), &scale, &sum);
    }
    const d = scale * @sqrt(sum);
    r.separations = .{ d, d };
}

/// The estimates `job` asks for, the cluster the leading `m` rows.
pub fn estimate(allocator: Allocator, io: Io, comptime T: type, p: Pencil(T), m: usize, job: Job, r: *Reordered(Float(T))) Error!void {
    const F = Float(T);
    const n1 = m;
    const n2 = p.n - m;
    const nn = n1 * n2;
    const work = try allocator.alloc(T, 2 * nn);
    defer allocator.free(work);
    @memset(work, ops.zero(T));
    const g: Gsep(T) = .{ .al = allocator, .io = io, .p = p, .n1 = n1, .n2 = n2 };
    if (wantsProjections(job)) {
        // `(A11, B11) R - L (A22, B22) = scale (A12, B12)`
        for (0..n2) |j| {
            @memcpy(work[j * n1 ..][0..n1], (p.a + (n1 + j) * p.lda)[0..n1]);
            @memcpy(work[nn + j * n1 ..][0..n1], (p.b + (n1 + j) * p.ldb)[0..n1]);
        }
        const sol = try g.solve(.n, 0, false, work.ptr);
        r.projections = .{ projection(F, sol.scale, frobOf(T, work[0..nn])), projection(F, sol.scale, frobOf(T, work[nn..])) };
    }
    if (wantsFrobenius(job)) {
        const s1 = try g.solve(.n, 3, false, work.ptr);
        const s2 = try g.solve(.n, 3, true, work.ptr);
        r.separations = .{ s1.separation orelse unreachable, s2.separation orelse unreachable };
    } else if (wantsOneNorm(job)) {
        const v = try allocator.alloc(T, 2 * nn);
        defer allocator.free(v);
        var d: [2]F = undefined;
        for ([_]bool{ false, true }, 0..) |swapped, i| {
            var last: Last(T) = .{ .g = g, .swapped = swapped };
            const est = norm.lacn2(allocator, T, Error, 2 * nn, v, work, &last, Last(T).apply) catch |e| switch (e) {
                error.OutOfMemory => return error.OutOfMemory,
                error.Canceled => return error.Canceled,
                error.SliceTooShort => unreachable,
            };
            d[i] = last.scale / est;
        }
        r.separations = d;
    }
}

/// `PL` or `PR` from `‖R‖_F` (`pl`) and `scale`: `1 / sqrt(1 + ‖R / scale‖²)`,
/// as `dtgsen.f` writes it.
fn projection(comptime F: type, scale: F, pl: F) F {
    if (pl == 0) return 1;
    return scale / (@sqrt(scale * scale / pl + pl) * @sqrt(pl));
}

/// `‖x‖₂` by `lassq`, in memory order.
fn frobOf(comptime T: type, x: []const T) Float(T) {
    var scale: Float(T) = 0;
    var sum: Float(T) = 1;
    norm.lassq(T, .of(x), &scale, &sum);
    return scale * @sqrt(sum);
}

/// The generalized Sylvester operator of the cluster `(A11, B11)` (`n1`)
/// and the rest `(A22, B22)` (`n2`), or with `swapped` the rest's and the
/// cluster's; its `R` in `x[0 .. n1 n2]`, `L` after (leading dimension the
/// first block's order).
fn Gsep(comptime T: type) type {
    return struct {
        al: Allocator,
        io: Io,
        p: Pencil(T),
        n1: usize,
        n2: usize,

        const Self = @This();

        fn solve(s: Self, trans: types.Trans, ijob: usize, swapped: bool, x: [*]T) Error!tgsyl_mod.GeneralizedSylvesterSolution(Float(T)) {
            const p = s.p;
            const k1 = if (swapped) s.n2 else s.n1;
            const k2 = if (swapped) s.n1 else s.n2;
            const o1: usize = if (swapped) s.n1 else 0;
            const o2: usize = if (swapped) 0 else s.n1;
            const a1 = (p.a + o1 + o1 * p.lda)[0 .. (k1 - 1) * p.lda + k1];
            const a2 = (p.a + o2 + o2 * p.lda)[0 .. (k2 - 1) * p.lda + k2];
            const b1 = (p.b + o1 + o1 * p.ldb)[0 .. (k1 - 1) * p.ldb + k1];
            const b2 = (p.b + o2 + o2 * p.ldb)[0 .. (k2 - 1) * p.ldb + k2];
            const nn = k1 * k2;
            return tgsyl_mod.tgsyl(s.al, s.io, T, trans, ijob, k1, k2, a1, p.lda, a2, p.lda, x[0..nn], k1, b1, p.ldb, b2, p.ldb, x[nn .. 2 * nn], k1) catch |e| switch (e) {
                error.OutOfMemory => return error.OutOfMemory,
                error.Canceled => return error.Canceled,
                else => unreachable,
            };
        }
    };
}

/// `lacn2`'s callback, the operator's inverse (its adjoint's for
/// `adjoint`), keeping the last solve's `scale`; a perturbed solve counts,
/// as the reference ignores its `IERR`.
fn Last(comptime T: type) type {
    return struct {
        g: Gsep(T),
        swapped: bool,
        scale: Float(T) = 1,

        fn apply(l: *@This(), adjoint: bool, x: [*]T) Error!void {
            const tr: types.Trans = if (adjoint) (if (comptime isComplex(T)) .c else .t) else .n;
            l.scale = (try l.g.solve(tr, 0, l.swapped, x)).scale;
        }
    };
}
