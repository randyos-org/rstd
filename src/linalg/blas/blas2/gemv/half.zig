//! Half-precision `GEMV` kernels: `y = alpha*op(A)*x + beta*y`, `A` column-major
//! `m` by `n`. `shgemv` is `f16` in and `f32` out, `sbgemv` is `bf16` in and `f32`
//! out, `bgemv` is `bf16` throughout.
//!
//! From `kernel/generic/gemv_n.c` and `gemv_t.c` (OpenBLAS, `BGEMM`/`HFLOAT16`
//! builds): each result is accumulated in `f32` and rounded once on the store.
//! `beta == 0` never reads `y`, so a NaN there does not propagate. `alpha == 0`
//! is the interface's concern, not the kernel's. No arch asm implements these
//! slots (the x86_64, arm64 and power sources are not ported); this is the only
//! kernel. `bf16` widening and narrowing match `ext/bf16` (denormals flush,
//! round to nearest even).

const types = @import("../../types.zig");
const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");
const bf16 = @import("../../ext/bf16/__root__.zig");

const BlasLong = types.BlasLong;
const HF16 = types.HF16;
const BF16 = types.BF16;

fn wide(comptime T: type, v: T) f32 {
    return if (T == BF16) bf16.fromBf16(v) else @floatCast(v);
}

fn narrow(comptime T: type, v: f32) T {
    return if (T == BF16) bf16.toBf16(v) else @floatCast(v);
}

/// `E` is the type of `A` and `x`, `Y` that of `y`, `alpha` and `beta`.
fn Gemv(comptime E: type, comptime Y: type, comptime trans: bool) type {
    return struct {
        fn k(m: BlasLong, n: BlasLong, alpha: Y, a: [*]E, lda: BlasLong, x: [*]E, inc_x: BlasLong, beta: Y, y: [*]Y, inc_y: BlasLong) callconv(.c) i32 {
            const rows = if (trans) n else m; // y's length
            const inner = if (trans) m else n; // x's length
            const alpha_f = wide(Y, alpha);
            const beta_f = wide(Y, beta);
            var r: BlasLong = 0;
            while (r < rows) : (r += 1) {
                var acc: f32 = 0;
                var q: BlasLong = 0;
                while (q < inner) : (q += 1) {
                    const e = util.at(E, a, if (trans) q + r * lda else r + q * lda).*;
                    acc += wide(E, e) * wide(E, util.at(E, x, q * inc_x).*);
                }
                const yp = util.at(Y, y, r * inc_y);
                yp.* = narrow(Y, if (beta_f == 0) alpha_f * acc else alpha_f * acc + beta_f * wide(Y, yp.*));
            }
            return 0;
        }
    };
}

/// Slot prefix, element type, and output type.
const families = .{
    .{ "shgemv", HF16, f32 },
    .{ "bgemv", BF16, BF16 },
    .{ "sbgemv", BF16, f32 },
};

pub fn installGeneric(vt: *Vtable) void {
    inline for (families) |f| {
        @field(vt, f[0] ++ "_n") = &Gemv(f[1], f[2], false).k;
        @field(vt, f[0] ++ "_t") = &Gemv(f[1], f[2], true).k;
    }
}

// ---- tests -----------------------------------------------------------------

const std = @import("std");
const check = @import("../../blas1/check.zig");

const sizes = [_][2]usize{ .{ 1, 1 }, .{ 2, 3 }, .{ 3, 2 }, .{ 4, 4 }, .{ 7, 5 }, .{ 8, 8 }, .{ 9, 16 }, .{ 16, 9 }, .{ 17, 17 }, .{ 33, 4 }, .{ 4, 33 } };
const strides = [_][2]usize{ .{ 1, 1 }, .{ 2, 3 }, .{ 3, 2 } };
const alphas = [_]f32{ 1.25, -0.5 };
const betas = [_]f32{ 0, 1, 0.5 };
const cap = 33 * 3 + 2; // a vector's span at stride 3

fn fillAs(comptime T: type, buf: []T, seed: u64) void {
    var f: [400]f32 = undefined;
    check.fill(f32, f[0..buf.len], seed);
    for (buf, f[0..buf.len]) |*d, s| d.* = narrow(T, s);
}

fn checkHalf(comptime E: type, comptime Y: type, comptime trans: bool, kern: anytype) !void {
    for (sizes) |mn| for (strides) |inc| for (alphas) |alpha| for (betas) |beta| {
        const m = mn[0];
        const n = mn[1];
        const lda = m + 1;
        var a: [400]E = undefined; // lda * n at most 18 * 17
        var x: [cap]E = undefined;
        var y: [cap]Y = undefined;
        fillAs(E, &a, m * 7 + n);
        fillAs(E, &x, m + 3);
        fillAs(Y, &y, n + 5);
        const rows = if (trans) n else m;
        const inner = if (trans) m else n;
        // With beta 0, y must never be read: poison it.
        if (beta == 0) for (&y) |*v| {
            v.* = narrow(Y, std.math.nan(f32));
        };
        const orig = y;
        _ = kern(@intCast(m), @intCast(n), narrow(Y, alpha), &a, @intCast(lda), &x, @intCast(inc[0]), narrow(Y, beta), &y, @intCast(inc[1]));
        var touched: [cap]bool = @splat(false);
        for (0..rows) |r| {
            var acc: f64 = 0;
            for (0..inner) |q| {
                const e = if (trans) a[q + r * lda] else a[r + q * lda];
                acc += @as(f64, wide(E, e)) * wide(E, x[q * inc[0]]);
            }
            const yi = r * inc[1];
            touched[yi] = true;
            const prior: f64 = if (beta == 0) 0 else wide(Y, orig[yi]);
            const want = @as(f64, alpha) * acc + @as(f64, beta) * prior;
            const got: f64 = wide(Y, y[yi]);
            const tol: f64 = if (Y == BF16)
                0.005 * @abs(want) + 1e-4
            else
                check.tolerance(f32, 4 * inner + 4);
            if (!(@abs(want - got) <= tol)) {
                check.say("got {e}, want {e} within {e}: m={d} n={d} trans={} alpha={d} beta={d} row={d}: ", .{ got, want, tol, m, n, trans, alpha, beta, r });
                return check.case(error.TestExpectedApproxEqAbs, "half gemv", E, rows, inc);
            }
        }
        for (0..cap) |i| {
            if (touched[i]) continue;
            if (!std.mem.eql(u8, std.mem.asBytes(&y[i]), std.mem.asBytes(&orig[i]))) {
                check.say("untouched y[{d}] changed: m={d} n={d} trans={}: ", .{ i, m, n, trans });
                return check.case(error.TestUnexpectedResult, "half gemv", E, rows, inc);
            }
        }
    };
}

fn checkAll(vt: *const Vtable) !void {
    inline for (families) |f| {
        try checkHalf(f[1], f[2], false, @field(vt, f[0] ++ "_n").?);
        try checkHalf(f[1], f[2], true, @field(vt, f[0] ++ "_t").?);
    }
}

test "portable half-precision gemv kernels match plain arithmetic" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
}

test "every linked core's half-precision gemv kernels match plain arithmetic, where this machine runs them" {
    try check.eachCore(checkAll);
}
