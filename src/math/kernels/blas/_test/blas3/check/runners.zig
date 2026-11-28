//! The runners: each routine through its driver or its public interface.

const std = @import("std");
const types = @import("../../../types.zig");
const Vtable = @import("../../../Vtable.zig");
const blas = @import("../../../__root__.zig");
const l3 = @import("../../../blas3/core.zig");
const driver = @import("../../../blas3/driver.zig");
const driver3m = @import("../../../blas3/driver3m.zig");
const guard = @import("../../guard.zig");
const c3 = @import("../check.zig");

const isComplex = types.isComplex;
const testing = std.testing;
const Case = c3.Case;
const Routine = c3.Routine;
const Trans = c3.Trans;
const io = c3.io;

// ---- runners --------------------------------------------------------------------

/// Calls the driver of `r` on a fixed vtable through a poisoned allocator, so
/// a packing buffer read before it is written shows as NaN. The routine is
/// fixed at compile time, so a test builds only its own driver.
pub fn Driver(comptime r: Routine) type {
    return struct {
        vt: *const Vtable,

        pub fn run(self: @This(), comptime T: type, case: Case, a: [*]T, lda: usize, b: [*]T, ldb: usize, alpha: T, beta: T, c: [*]T, ldc: usize) !void {
            var poison: guard.Poison = .{ .child = testing.allocator };
            return call(T, poison.allocator(), io, self.vt, r, case, a, lda, b, ldb, alpha, beta, c, ldc, 1);
        }
    };
}

/// The one dispatcher of the level-3 tests: driver `r` on `vt` for `case`'s
/// options at up to `want` workers.
pub fn call(comptime T: type, al: std.mem.Allocator, run_io: std.Io, vt: *const Vtable, comptime r: Routine, case: Case, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, alpha: T, beta: T, c: [*]T, ldc: usize, want: u32) !void {
    const herm = comptime isComplex(T);
    switch (r) {
        .gemm => try driver.gemm(al, run_io, T, vt, case.ta.transposed(), case.tb.transposed(), .of(case.ta.conjugated(), case.tb.conjugated()), case.m, case.n, case.k, alpha, a, lda, b, ldb, beta, c, ldc, want),
        .gemmt => try driver.gemmt(al, run_io, T, vt, case.upper, case.ta.transposed(), case.tb.transposed(), .of(case.ta.conjugated(), case.tb.conjugated()), case.n, case.k, alpha, a, lda, b, ldb, beta, c, ldc, want),
        .symm => try driver.symm(al, run_io, T, vt, false, case.left, case.upper, case.m, case.n, alpha, a, lda, b, ldb, beta, c, ldc, want),
        .hemm => if (herm) try driver.symm(al, run_io, T, vt, true, case.left, case.upper, case.m, case.n, alpha, a, lda, b, ldb, beta, c, ldc, want) else unreachable,
        .syrk => try driver.syrk(al, run_io, T, vt, false, case.upper, case.ta.transposed(), case.n, case.k, alpha, a, lda, null, 0, beta, c, ldc, want),
        .syr2k => try driver.syrk(al, run_io, T, vt, false, case.upper, case.ta.transposed(), case.n, case.k, alpha, a, lda, b, ldb, beta, c, ldc, want),
        .herk => if (herm) try driver.syrk(al, run_io, T, vt, true, case.upper, case.ta.transposed(), case.n, case.k, alpha, a, lda, null, 0, beta, c, ldc, want) else unreachable,
        .her2k => if (herm) try driver.syrk(al, run_io, T, vt, true, case.upper, case.ta.transposed(), case.n, case.k, alpha, a, lda, b, ldb, beta, c, ldc, want) else unreachable,
        .trmm, .trsm => if (@hasDecl(driver, "trmm")) {
            const args = .{ al, run_io, T, vt, case.left, case.upper, case.ta.transposed(), case.ta.conjugated(), case.unit, case.m, case.n, alpha, a, lda, c, ldc, want };
            if (r == .trmm) try @call(.auto, driver.trmm, args) else try @call(.auto, driver.trsm, args);
        } else return error.SkipZigTest,
        .gemm3m => if (herm) try driver3m.gemm3m(al, run_io, T, vt, case.ta.transposed(), case.tb.transposed(), .of(case.ta.conjugated(), case.tb.conjugated()), case.m, case.n, case.k, alpha, a, lda, b, ldb, beta, c, ldc, want) else unreachable,
        .symm3m, .hemm3m => if (herm) try driver3m.symm3m(al, run_io, T, vt, r == .hemm3m, case.left, case.upper, case.m, case.n, alpha, a, lda, b, ldb, beta, c, ldc, want) else unreachable,
    }
}

/// `call` for a routine known only at run time.
pub fn callAny(comptime T: type, al: std.mem.Allocator, run_io: std.Io, vt: *const Vtable, r: Routine, case: Case, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, alpha: T, beta: T, c: [*]T, ldc: usize, want: u32) !void {
    switch (r) {
        inline else => |cr| try call(T, al, run_io, vt, cr, case, a, lda, b, ldb, alpha, beta, c, ldc, want),
    }
}

/// Calls the public routines (`l3.gemm` and the rest) through a poisoned
/// allocator; real `r` and `c` stand in for `n` and `t` half the time. As
/// `Driver`, for one routine `r`.
pub fn Iface(comptime r: Routine) type {
    return struct {
        pub fn run(_: @This(), comptime T: type, case: Case, a: [*]T, lda: usize, b: [*]T, ldb: usize, alpha: T, beta: T, c: [*]T, ldc: usize) !void {
            var poison: guard.Poison = .{ .child = testing.allocator };
            const al = poison.allocator();
            const odd = (case.m + case.n + case.k) % 2 == 1;
            const ta = trans(T, case.ta, odd, case.rank());
            const tb = trans(T, case.tb, !odd, false);
            const m: isize = @intCast(case.m);
            const n: isize = @intCast(case.n);
            const k: isize = @intCast(case.k);
            const ia: isize = @intCast(lda);
            const ib: isize = @intCast(ldb);
            const ic: isize = @intCast(ldc);
            const side: blas.Side = if (case.left) .l else .r;
            const uplo: blas.Uplo = if (case.upper) .u else .l;
            const diag: blas.Diag = if (case.unit) .u else .n;
            const herm = comptime isComplex(T);
            switch (r) {
                .gemm => try l3.gemm(al, io, T, ta, tb, m, n, k, alpha, a, ia, b, ib, beta, c, ic),
                .gemmt => try l3.gemmt(al, io, T, uplo, ta, tb, n, k, alpha, a, ia, b, ib, beta, c, ic),
                .symm => try l3.symm(al, io, T, side, uplo, m, n, alpha, a, ia, b, ib, beta, c, ic),
                .hemm => if (herm) try l3.hemm(al, io, T, side, uplo, m, n, alpha, a, ia, b, ib, beta, c, ic) else unreachable,
                .syrk => try l3.syrk(al, io, T, uplo, ta, n, k, alpha, a, ia, beta, c, ic),
                .syr2k => try l3.syr2k(al, io, T, uplo, ta, n, k, alpha, a, ia, b, ib, beta, c, ic),
                .herk => if (herm) try l3.herk(al, io, T, uplo, ta, n, k, alpha.re, a, ia, beta.re, c, ic) else unreachable,
                .her2k => if (herm) try l3.her2k(al, io, T, uplo, ta, n, k, alpha, a, ia, b, ib, beta.re, c, ic) else unreachable,
                .trmm => try l3.trmm(al, io, T, side, uplo, ta, diag, m, n, alpha, a, ia, c, ic),
                .trsm => try l3.trsm(al, io, T, side, uplo, ta, diag, m, n, alpha, a, ia, c, ic),
                .gemm3m => if (herm) try l3.gemm3m(al, io, T, ta, tb, m, n, k, alpha, a, ia, b, ib, beta, c, ic) else unreachable,
                .symm3m => if (herm) try l3.symm3m(al, io, T, side, uplo, m, n, alpha, a, ia, b, ib, beta, c, ic) else unreachable,
                .hemm3m => if (herm) try l3.hemm3m(al, io, T, side, uplo, m, n, alpha, a, ia, b, ib, beta, c, ic) else unreachable,
            }
        }

        /// `t` as `blas.Trans`; real `r` and `c` for `n` and `t` where `alt`.
        pub fn trans(comptime T: type, t: Trans, alt: bool, rank: bool) blas.Trans {
            if (comptime isComplex(T)) return switch (t) {
                inline else => |v| @field(blas.Trans, @tagName(v)),
            };
            // the rank-k updates take no `r`
            return if (t.transposed()) (if (alt) .c else .t) else (if (alt and !rank) .r else .n);
        }
    };
}
