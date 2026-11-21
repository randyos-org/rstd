//! Tests of the SME2 direct real `gemm` (`../driver.zig`, `../kernel.zig`),
//! where this CPU runs it (FEAT_SME2, 512-bit streaming vectors) and the
//! `vortexm4` object is linked; skipped elsewhere.

const std = @import("std");
const sysinfo = @import("builtin");
const blas = @import("../../../../__root__.zig");
const Vtable = @import("../../../../Vtable.zig");
const check = @import("../../../../_test/check.zig");
const bitwise = @import("../../../../_test/bitwise.zig");
const driver = @import("../driver.zig");
const job = @import("../job.zig");

const Allocator = std.mem.Allocator;

/// `vortexm4`'s vtable where its SME2 kernels run here.
fn sme2Vtable() ?Vtable {
    if (comptime sysinfo.target.cpu.arch != .aarch64) return null;
    const vt = blas.cpu.vtableFromCoreType(.vortexm4);
    if (vt.sme_sgemm_kernel == null or vt.sme_dgemm_kernel == null or !driver.usable()) return null;
    return vt;
}

/// `op(X)`'s `(i, l)` of a column-major `X` at `ld`.
fn at(comptime T: type, x: []const T, t: bool, i: usize, l: usize, ld: usize) T {
    return if (t) x[l + i * ld] else x[i + l * ld];
}

fn operands(comptime T: type, al: Allocator, t: bool, rows: usize, cols: usize) ![]T {
    // stored `cols x rows` when transposed, one spare row each
    return al.alloc(T, (if (t) cols else rows) * (if (t) rows else cols) + (if (t) rows else cols) + 1);
}

/// One `gemm` of integer operands (exact in either type) against the exact
/// sum, every transposition and the special scalars.
fn exact(comptime T: type, vt: *const Vtable, m: usize, n: usize, k: usize) !void {
    const al = std.testing.allocator;
    for ([_][2]bool{ .{ false, false }, .{ true, false }, .{ false, true }, .{ true, true } }) |tr| {
        const lda = (if (tr[0]) k else m) + 1;
        const ldb = (if (tr[1]) n else k) + 1;
        const ldc = m + 3;
        const a = try al.alloc(T, lda * (if (tr[0]) m else k));
        defer al.free(a);
        const b = try al.alloc(T, ldb * (if (tr[1]) k else n));
        defer al.free(b);
        const c = try al.alloc(T, ldc * n);
        defer al.free(c);
        const c0 = try al.alloc(T, ldc * n);
        defer al.free(c0);
        check.fillSeq(T, a, m + 2 * k);
        check.fillSeq(T, b, n + 3 * k);
        check.fillSeq(T, c0, m + n);
        for ([_][2]T{ .{ 1, 0 }, .{ -1, 1 }, .{ 0.75, -0.5 }, .{ 1, 1 } }) |s| {
            @memcpy(c, c0);
            if (!try driver.gemm(al, std.Io.failing, T, vt, tr[0], tr[1], m, n, k, s[0], a.ptr, lda, b.ptr, ldb, s[1], c.ptr, ldc, 1)) return error.TestUnexpectedResult;
            for (0..n) |j| for (0..ldc) |i| {
                const got = c[i + j * ldc];
                if (i >= m) {
                    // padding untouched
                    if (!bitwise.sameBits(T, c[i + j * ldc ..][0..1], c0[i + j * ldc ..][0..1])) return error.TestUnexpectedResult;
                    continue;
                }
                var sum: f64 = 0;
                for (0..k) |l| sum += @as(f64, at(T, a, tr[0], i, l, lda)) * at(T, b, tr[1], l, j, ldb);
                const want: f64 = s[0] * sum + s[1] * @as(f64, c0[i + j * ldc]);
                if (got != want) {
                    std.debug.print("sme2 {s} ta={} tb={} m={d} n={d} k={d} alpha={d} beta={d} at ({d},{d}): {d} not {d}\n", .{ @typeName(T), tr[0], tr[1], m, n, k, s[0], s[1], i, j, got, want });
                    return error.TestUnexpectedResult;
                }
            };
        }
    }
}

test "the SME2 gemm is exact on integer operands at every tile edge" {
    var vt = sme2Vtable() orelse return error.SkipZigTest;
    inline for (.{ f32, f64 }) |T| {
        const vl = job.Shape(T).vl;
        // rows and columns of `C` around the vector, the half and full rows
        // of tiles, and the edge kernel's stacks; depths past the 4-step loop
        const sizes = [_]usize{ 1, vl - 1, vl + 1, 2 * vl + 3, 33, 63, 64, 65, 97, 130, 200 };
        for (sizes, 0..) |m, mi| for (sizes, 0..) |n, ni| for ([_]usize{ 1, 3, 9, 70 }, 0..) |k, ki| {
            if ((mi + ni + ki) % 3 != 0) continue;
            try exact(T, &vt, m, n, k);
        };
        // deeper than one `kc` slice and wider than one `nc` block
        try exact(T, &vt, 70, 1100, 1300);
        try exact(T, &vt, 1100, 70, 1300);
    }
}

/// The SME2 job of `C^T` over `blk`, through the kernel slot.
fn runJob(comptime T: type, vt: *const Vtable, al: Allocator, jb0: job.Job(T), blk: job.Blocking) !void {
    var jb = jb0;
    jb.blk = blk;
    const w = job.Work.of(T, blk);
    const buf = try al.alignedAlloc(T, .fromByteUnits(128), w.len);
    defer al.free(buf);
    @memset(buf, std.math.nan(T));
    jb.ac = buf.ptr;
    jb.bc = buf.ptr + w.a_len;
    (if (T == f32) vt.sme_sgemm_kernel else vt.sme_dgemm_kernel).?(&jb);
}

test "the SME2 gemm keeps its bits for any blocking and any number of parts" {
    var vt = sme2Vtable() orelse return error.SkipZigTest;
    const al = std.testing.allocator;
    const saved = driver.test_units;
    defer driver.test_units = saved;
    inline for (.{ f32, f64 }) |T| {
        const m = 301;
        const n = 263;
        const k = 197;
        const a = try al.alloc(T, m * k);
        defer al.free(a);
        const b = try al.alloc(T, k * n);
        defer al.free(b);
        const c0 = try al.alloc(T, m * n);
        defer al.free(c0);
        check.fillMixed(T, a, 11);
        check.fillMixed(T, b, 12);
        check.fillMixed(T, c0, 13);
        for ([_][2]bool{ .{ false, false }, .{ true, false }, .{ false, true }, .{ true, true } }) |tr| {
            const lda: usize = if (tr[0]) k else m;
            const ldb: usize = if (tr[1]) n else k;
            for ([_]T{ 0, 1, -0.375 }) |beta| {
                const want = try al.dupe(T, c0);
                defer al.free(want);
                const alpha: T = -1.25;
                driver.test_units = 1;
                _ = try driver.gemm(al, std.Io.failing, T, &vt, tr[0], tr[1], m, n, k, alpha, a.ptr, lda, b.ptr, ldb, beta, want.ptr, m, 1);
                // the same job on blockings of its smallest blocks
                const whole: job.Job(T) = .{ .m = n, .n = m, .k = k, .alpha = alpha, .beta = beta, .a = b.ptr, .lda = ldb, .a_cols = tr[1], .b = a.ptr, .ldb = lda, .b_cols = tr[0], .c = undefined, .ldc = m, .blk = undefined, .ac = undefined, .bc = undefined };
                const S = job.Shape(T);
                for ([_]job.Blocking{ .{ .mc = S.vl, .nc = S.nr, .kc = 16 }, .{ .mc = 3 * S.vl, .nc = 2 * S.nr, .kc = 1 }, .{ .mc = 304, .nc = 320, .kc = 200 } }) |blk| {
                    const got = try al.dupe(T, c0);
                    defer al.free(got);
                    var jb = whole;
                    jb.c = got.ptr;
                    try runJob(T, &vt, al, jb, blk);
                    try std.testing.expect(bitwise.sameBits(T, want, got));
                }
                // and on 2 to 4 parts
                for ([_]usize{ 2, 3, 4 }) |u| {
                    driver.test_units = u;
                    const got = try al.dupe(T, c0);
                    defer al.free(got);
                    _ = try driver.gemm(al, std.Io.failing, T, &vt, tr[0], tr[1], m, n, k, alpha, a.ptr, lda, b.ptr, ldb, beta, got.ptr, m, 4);
                    try std.testing.expect(bitwise.sameBits(T, want, got));
                }
            }
        }
    }
}

test "a gemm too small for streaming mode takes the plain loops" {
    var vt = sme2Vtable() orelse return error.SkipZigTest;
    inline for (.{ f32, f64 }) |T| {
        try exact(T, &vt, 7, 5, 3);
        try exact(T, &vt, 1, 1, 1);
        try exact(T, &vt, 14, 14, 15);
    }
}
