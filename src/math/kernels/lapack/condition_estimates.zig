//! Reciprocal condition numbers from a triangular matrix or a Cholesky
//! factor, `?TRCON` and `?POCON`, Reference-LAPACK 3.12.0
//! (`lapack-netlib/SRC/dtrcon.f`, `ztrcon.f`, `dpocon.f`, `zpocon.f`,
//! OpenBLAS v0.3.34-317; modified BSD): `lacn2` over `latrs` solves, as
//! `gecon.zig`.

const std = @import("std");
const types = @import("../blas/types.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");
const norms = @import("norm.zig");
const structured = @import("norm_structured.zig");
const lacn2_estimate = @import("lacn2.zig");
const latrs_solve = @import("latrs.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const Trans = types.Trans;
const Uplo = types.Uplo;
const Diag = types.Diag;
const Error = ops.Error;
const NormType = norms.NormType;
const lacn2 = lacn2_estimate.lacn2;
const rscl = latrs_solve.rscl;
const scaledSolve = latrs_solve.scaledSolve;

/// What `trcon` returns: `UnsupportedNorm` for other than `.one` or
/// `.infinity`, a bad `uplo` or shape, and the allocator's and
/// cancellation.
pub const TrconError = error{UnsupportedNorm} || interface.TriangleError || Allocator.Error || Io.Cancelable;

/// What `pocon` returns: `NormOutOfRange` for a negative `anorm` (the
/// reference's `INFO = -5`), a bad `uplo` or shape, and the allocator's and
/// cancellation.
pub const PoconError = error{NormOutOfRange} || interface.TriangleError || Allocator.Error || Io.Cancelable;

/// `?TRCON` (`dtrcon.f`, `ztrcon.f`): an estimate of `1 / (‖A‖ ‖A⁻¹‖)` in
/// the one or infinity norm for the `n x n` triangular `A` (`uplo`, `diag`);
/// `‖A‖` by `lantr`, `‖A⁻¹‖` by `lacn2`. 1 for `n = 0`; 0 for a zero (or
/// NaN) `‖A‖`, a solve that would overflow, or a zero estimate.
pub fn trcon(allocator: Allocator, io: Io, comptime T: type, norm: NormType, uplo: Uplo, diag: Diag, n: usize, a: []const T, lda: usize) TrconError!Float(T) {
    comptime ops.requireFloat(T, "trcon");
    const F = Float(T);
    const onenrm = switch (norm) {
        .one => true,
        .infinity => false,
        .max_abs, .frobenius => return error.UnsupportedNorm,
    };
    try interface.holdsTriangle(uplo, n, a.len, lda);
    if (n == 0) return 1;
    const anorm = structured.lantr(allocator, io, T, norm, uplo, diag, n, n, a, lda) catch |e| switch (e) {
        error.OutOfMemory => return error.OutOfMemory,
        // checked above
        error.UnsupportedUplo, error.SliceTooShort, error.LeadingDimensionTooSmall => unreachable,
    };
    if (!(anorm > 0)) return 0;
    const smlnum = std.math.floatMin(F) * @as(F, @floatFromInt(n));
    var ctx: Solves(T) = .{ .al = allocator, .io = io, .kind = .{ .triangular = .{ .onenrm = onenrm, .unit = diag == .u } }, .upper = uplo == .u, .n = n, .a = a.ptr, .lda = lda, .smlnum = smlnum, .cnorm = undefined };
    const ainvnm = try estimate(allocator, T, &ctx);
    if (ainvnm == 0) return 0;
    return (1 / anorm) / ainvnm;
}

/// `?POCON` (`dpocon.f`, `zpocon.f`): an estimate of `1 / (‖A‖₁ ‖A⁻¹‖₁)`
/// for the positive definite `n x n` `A` from `potrf`'s factor in `a`'s
/// `uplo` triangle and `anorm`, `‖A‖₁` before factoring. 1 for `n = 0`; 0
/// for `anorm = 0`, a solve that would overflow, or a zero estimate.
pub fn pocon(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: []const T, lda: usize, anorm: Float(T)) PoconError!Float(T) {
    comptime ops.requireFloat(T, "pocon");
    try interface.holdsTriangle(uplo, n, a.len, lda);
    if (anorm < 0) return error.NormOutOfRange;
    if (n == 0) return 1;
    if (anorm == 0) return 0;
    var ctx: Solves(T) = .{ .al = allocator, .io = io, .kind = .cholesky, .upper = uplo == .u, .n = n, .a = a.ptr, .lda = lda, .smlnum = std.math.floatMin(Float(T)), .cnorm = undefined };
    const ainvnm = try estimate(allocator, T, &ctx);
    if (ainvnm == 0) return 0;
    return (1 / ainvnm) / anorm;
}

/// `lacn2`'s estimate of `‖A⁻¹‖`, 0 where the reference's `GO TO 20` gives
/// up on an overflow.
fn estimate(al: Allocator, comptime T: type, ctx: *Solves(T)) (Allocator.Error || Io.Cancelable)!Float(T) {
    const n = ctx.n;
    const work = try al.alloc(T, 2 * n);
    defer al.free(work);
    const cnorm = try al.alloc(Float(T), n);
    defer al.free(cnorm);
    ctx.cnorm = cnorm.ptr;
    return lacn2(al, T, Solves(T).Stop, n, work[n..], work[0..n], ctx, Solves(T).apply) catch |e| switch (e) {
        error.Overflow => 0,
        error.OutOfMemory => error.OutOfMemory,
        error.Canceled => error.Canceled,
        // sized above
        error.SliceTooShort => unreachable,
    };
}

/// `A⁻¹` (its adjoint's) for `lacn2`: one `latrs` solve with the triangle
/// (`trcon`, its transpose as `KASE` asks), or two with the Cholesky
/// factor (`pocon`); a scale that would overflow ends the estimate.
fn Solves(comptime T: type) type {
    return struct {
        al: Allocator,
        io: Io,
        kind: union(enum) { triangular: struct { onenrm: bool, unit: bool }, cholesky },
        upper: bool,
        n: usize,
        a: [*]const T,
        lda: usize,
        smlnum: Float(T),
        cnorm: [*]Float(T),
        normin: bool = false,

        const Self = @This();
        const Stop = Error || error{Overflow};

        fn apply(s: *Self, adjoint: bool, x: [*]T) Stop!void {
            const op: Trans = if (comptime isComplex(T)) .c else .t;
            var scale: Float(T) = undefined;
            switch (s.kind) {
                .triangular => |t| {
                    // `KASE = KASE1` (`adjoint != onenrm`) solves with `A` itself
                    const tr: Trans = if (adjoint != t.onenrm) .n else op;
                    scale = try scaledSolve(s.al, s.io, T, s.upper, tr, t.unit, s.normin, s.n, s.a, s.lda, x, s.cnorm);
                    s.normin = true;
                },
                .cholesky => {
                    // `Uᴴ` then `U`, or `L` then `Lᴴ`
                    const first: Trans = if (s.upper) op else .n;
                    const second: Trans = if (s.upper) .n else op;
                    const sl = try scaledSolve(s.al, s.io, T, s.upper, first, false, s.normin, s.n, s.a, s.lda, x, s.cnorm);
                    s.normin = true;
                    const su = try scaledSolve(s.al, s.io, T, s.upper, second, false, true, s.n, s.a, s.lda, x, s.cnorm);
                    scale = sl * su;
                },
            }
            if (scale != 1) {
                const ix = blas_calls.iamax(T, s.n, x, 1);
                if (scale < ops.abs1(T, x[ix]) * s.smlnum or scale == 0) return error.Overflow;
                try rscl(s.io, T, scale, .of(x[0..s.n]));
            }
        }
    };
}
