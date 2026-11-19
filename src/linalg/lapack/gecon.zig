//! `xGECON`: the reciprocal condition number from `getrf`'s factor, Reference-LAPACK 3.12.0
//! (`lapack-netlib/SRC/dgecon.f`, `zgecon.f`, OpenBLAS v0.3.34-317; modified BSD).

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");
const norms = @import("norm.zig");
const lacn2_estimate = @import("lacn2.zig");
const latrs_solve = @import("latrs.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const Trans = types.Trans;
const Error = ops.Error;
const Norm = norms.Norm;
const lacn2 = lacn2_estimate.lacn2;
const rscl = latrs_solve.rscl;
const scaledSolve = latrs_solve.scaledSolve;

/// What `gecon` returns: `UnsupportedNorm` for other than `o` or `i`; a
/// slice too short or a leading dimension too small; `NormOutOfRange` for a
/// negative, infinite or NaN `anorm` (the reference's `INFO = -5`, the last
/// two after its checks); `EstimateFailed` for the reference's `INFO = 1`,
/// an estimate that is NaN, infinite, or of a zero inverse norm.
pub const GeconError = error{ UnsupportedNorm, NormOutOfRange, EstimateFailed } || interface.ShapeError || Allocator.Error || Io.Cancelable;

/// `?GECON` (`dgecon.f`, `zgecon.f`): an estimate of the reciprocal
/// condition number `1 / (‖A‖ ‖A⁻¹‖)` in the one (`o`) or infinity (`i`)
/// norm, from `getrf`'s `L U` in `A` (`n x n`) and `anorm`, `A`'s norm
/// before factoring. `‖A⁻¹‖` comes from `lacn2` over solves with `L` and
/// `U` (`latrs`). Returns `rcond`: 1 for `n = 0`, 0 for `anorm = 0` or a
/// solve that would overflow, the estimate otherwise. Takes `5 n` elements
/// of workspace.
pub fn gecon(allocator: Allocator, io: Io, comptime T: type, norm: Norm, n: usize, a: []const T, lda: usize, anorm: Float(T)) GeconError!Float(T) {
    comptime ops.requireFloat(T, "gecon");
    const F = Float(T);
    const onenrm = switch (norm) {
        .o => true,
        .i => false,
        .m, .f => return error.UnsupportedNorm,
    };
    try util.holds(a.len, n, n, lda);
    if (anorm < 0) return error.NormOutOfRange;
    if (n == 0) return 1;
    if (anorm == 0) return 0;
    if (std.math.isNan(anorm) or anorm > std.math.floatMax(F)) return error.NormOutOfRange;
    const un = n;
    const work = try allocator.alloc(T, 2 * un);
    defer allocator.free(work);
    const cnorm = try allocator.alloc(F, 2 * un);
    defer allocator.free(cnorm);
    var ctx: Gecon(T) = .{ .al = allocator, .io = io, .onenrm = onenrm, .n = un, .a = a.ptr, .lda = lda, .cnorm = cnorm.ptr };
    const ainvnm = lacn2(allocator, T, Gecon(T).Stop, un, work[un..], work[0..un], &ctx, Gecon(T).apply) catch |e| switch (e) {
        error.Overflow => return 0,
        error.OutOfMemory, error.Canceled => |x| return x,
        // sized above
        error.SliceTooShort => unreachable,
    };
    if (ainvnm == 0) return error.EstimateFailed;
    const rcond = (1 / ainvnm) / anorm;
    if (std.math.isNan(rcond) or rcond > std.math.floatMax(F)) return error.EstimateFailed;
    return rcond;
}

/// `gecon`'s `B = A⁻¹` for `lacn2`: two `latrs` solves, `L` then `U` (or
/// their transposes in reverse), and the reference's early exit (`GO TO
/// 20`) as `error.Overflow`.
fn Gecon(comptime T: type) type {
    return struct {
        al: Allocator,
        io: Io,
        onenrm: bool,
        n: usize,
        a: [*]const T,
        lda: usize,
        cnorm: [*]Float(T),
        normin: bool = false,

        const Self = @This();
        const Stop = Error || error{Overflow};

        fn apply(s: *Self, adjoint: bool, x: [*]T) Stop!void {
            const F = Float(T);
            const op: Trans = if (comptime isComplex(T)) .c else .t;
            const lower_norms = s.cnorm;
            const upper_norms = s.cnorm + s.n;
            var sl: F = undefined;
            var su: F = undefined;
            if (adjoint != s.onenrm) {
                sl = try scaledSolve(s.al, s.io, T, false, .n, true, s.normin, s.n, s.a, s.lda, x, lower_norms);
                su = try scaledSolve(s.al, s.io, T, true, .n, false, s.normin, s.n, s.a, s.lda, x, upper_norms);
            } else {
                su = try scaledSolve(s.al, s.io, T, true, op, false, s.normin, s.n, s.a, s.lda, x, upper_norms);
                sl = try scaledSolve(s.al, s.io, T, false, op, true, s.normin, s.n, s.a, s.lda, x, lower_norms);
            }
            const scale = sl * su;
            s.normin = true;
            if (scale != 1) {
                const ix = blas_calls.iamax(T, s.n, x, 1);
                if (scale < ops.abs1(T, x[ix]) * std.math.floatMin(F) or scale == 0) return error.Overflow;
                try rscl(s.io, T, scale, .of(x[0..s.n]));
            }
        }
    };
}
