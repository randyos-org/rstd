//! Tests of `../../blas3/gemm/small.zig`.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../../types.zig");
const fma = @import("../../fma.zig");
const Vtable = @import("../../Vtable.zig");
const portable = @import("../../portable.zig");
const gk = @import("../../blas3/gemm/kernel.zig");
const params = @import("../../blas3/gemm/params.zig");
const small_kernel = @import("../../blas3/gemm/small_kernel.zig");
const sroot = @import("../../blas3/gemm/small.zig");

const BlasLong = types.BlasLong;
const isComplex = types.isComplex;
const HalfPermit = sroot.HalfPermit;
const Permit = sroot.Permit;
const allows = sroot.allows;
const allowsHalf = sroot.allowsHalf;
const halfPermitFor = sroot.halfPermitFor;
const permitFor = sroot.permitFor;

// ---- tests -----------------------------------------------------------------

const check3 = @import("check.zig");
const check = @import("../check.zig");

const fillSeq = check.fillSeq;

/// Every kernel against plain arithmetic, around every vector and column
/// group edge; exact on small integers. `C` padding must survive.
fn checkKernels(comptime T: type, comptime geo: gk.Geometry) !void {
    const K = small_kernel.Kernels(T, geo);
    const complex = comptime isComplex(T);
    const F = types.Float(T);
    const cs: usize = if (complex) 2 else 1;
    const al = std.testing.allocator;
    const mr: usize = geo.mr();
    const sizes = [_]usize{ 1, 2, geo.vl + 1, mr, mr + geo.vl + 1, geo.nr, geo.nr + 2 };
    const ld = std.mem.max(usize, &sizes) + 2;
    const buf = try al.alloc(F, 5 * cs * ld * ld);
    defer al.free(buf);
    const a = buf[0 .. cs * ld * ld];
    const b = buf[cs * ld * ld ..][0 .. cs * ld * ld];
    const c = buf[2 * cs * ld * ld ..][0 .. cs * ld * ld];
    const want = buf[3 * cs * ld * ld ..][0 .. cs * ld * ld];
    const got = buf[4 * cs * ld * ld ..][0 .. cs * ld * ld];
    fillSeq(F, a, 1);
    fillSeq(F, b, 2);
    fillSeq(F, c, 3);
    const ops: []const u8 = if (complex) "ntrc" else "nt";
    const lm: BlasLong = @intCast(ld);
    inline for (comptime ops) |oa| inline for (comptime ops) |ob| inline for (.{ false, true }) |b0| {
        const ta = oa == 't' or oa == 'c';
        const tb = ob == 't' or ob == 'c';
        const sa: F = if (oa == 'r' or oa == 'c') -1 else 1;
        const sb: F = if (ob == 'r' or ob == 'c') -1 else 1;
        const name = "gemm_small_kernel_" ++ (if (b0) "b0_" else "") ++ &[_]u8{ oa, ob };
        const f = @field(K, name);
        for (sizes) |m| for (sizes) |n| for ([_]usize{ 1, 5 }) |k| {
            @memcpy(want, c);
            @memcpy(got, c);
            for (0..n) |j| for (0..m) |i| {
                var sr: F = 0;
                var si: F = 0;
                for (0..k) |l| {
                    const ia = cs * (if (ta) l + i * ld else i + l * ld);
                    const ib = cs * (if (tb) j + l * ld else l + j * ld);
                    const x = a[ia];
                    const y = if (complex) sa * a[ia + 1] else 0;
                    const u = b[ib];
                    const v = if (complex) sb * b[ib + 1] else 0;
                    sr += x * u - y * v;
                    si += x * v + y * u;
                }
                // alpha 2 (2 - i), beta 0 or 3 (3 + 2i)
                const p = cs * (i + j * ld);
                if (complex) {
                    const cr = if (b0) 0 else 3 * c[p] - 2 * c[p + 1];
                    const ci = if (b0) 0 else 3 * c[p + 1] + 2 * c[p];
                    want[p] = cr + 2 * sr + si;
                    want[p + 1] = ci + 2 * si - sr;
                } else want[p] = (if (b0) 0 else 3 * c[p]) + 2 * sr;
            };
            if (complex) {
                if (b0) _ = f(@intCast(m), @intCast(n), @intCast(k), a.ptr, lm, 2, -1, b.ptr, lm, got.ptr, lm) else _ = f(@intCast(m), @intCast(n), @intCast(k), a.ptr, lm, 2, -1, b.ptr, lm, 3, 2, got.ptr, lm);
            } else {
                if (b0) _ = f(@intCast(m), @intCast(n), @intCast(k), a.ptr, lm, 2, b.ptr, lm, got.ptr, lm) else _ = f(@intCast(m), @intCast(n), @intCast(k), a.ptr, lm, 2, b.ptr, lm, 3, got.ptr, lm);
            }
            if (!std.mem.eql(F, want, got)) {
                std.debug.print("{s} {s} m={d} n={d} k={d}\n", .{ @typeName(T), name, m, n, k });
                return error.TestUnexpectedResult;
            }
        };
    };
}

test "small kernels match plain arithmetic around every vector and column edge" {
    inline for (.{ f32, f64 }) |T| try checkKernels(T, gk.geometry(T));
    try checkKernels(f64, .{ .vl = 1, .mv = 3, .nr = 5 });
    try checkKernels(f32, .{ .vl = 4, .mv = 2, .nr = 3 });
    inline for (.{ types.Complex32, types.Complex64 }) |T| try checkKernels(T, gk.geometry(T));
}

test "small permits: the reference's thresholds" {
    // Skylake-X: past 100^3 never; `TN` not for a big `C` or a short depth
    try std.testing.expect(allows(.skylakex, true, false, false, 100, 100, 100));
    try std.testing.expect(!allows(.skylakex, true, false, false, 100, 100, 101));
    try std.testing.expect(!allows(.skylakex, true, true, false, 40, 31, 64));
    try std.testing.expect(allows(.skylakex, true, true, false, 40, 30, 32));
    try std.testing.expect(!allows(.skylakex, false, true, true, 4, 4, 3));
    try std.testing.expect(allows(.skylakex, true, true, true, 4, 4, 3));
    try std.testing.expect(allows(.sve, false, true, true, 64, 64, 64));
    try std.testing.expect(!allows(.sve, true, false, false, 64, 64, 65));
    try std.testing.expect(!allows(.power10, true, false, false, 55, 54, 54));
    try std.testing.expect(allows(.power10, false, false, false, 64, 64, 64));
    try std.testing.expect(!allows(.la464, false, true, false, 15, 14, 14));
    try std.testing.expect(allows(.la264, false, true, false, 80, 80, 80));
    try std.testing.expect(!allows(.generic, false, false, false, 1, 1, 1));
    try std.testing.expectEqual(Permit.skylakex, permitFor(.x86_64, .sapphirerapids));
    try std.testing.expectEqual(Permit.generic, permitFor(.x86, .skylakex));
    try std.testing.expectEqual(Permit.sve, permitFor(.arm64, .neoversev1));
    try std.testing.expectEqual(Permit.generic, permitFor(.arm64, .neoversev2));
    // Cooper Lake's `sbgemm`: 8, 16 or 32 columns, up to 256^3
    try std.testing.expect(allowsHalf(.cooperlake, 256, 32, 2048));
    try std.testing.expect(!allowsHalf(.cooperlake, 256, 32, 2049));
    try std.testing.expect(!allowsHalf(.cooperlake, 4, 9, 4));
    try std.testing.expectEqual(HalfPermit.cooperlake, halfPermitFor(.x86_64, .cooperlake));
    try std.testing.expectEqual(HalfPermit.generic, halfPermitFor(.x86_64, .sapphirerapids));
    // a core's blocking brings its permit
    var vt: Vtable = .{};
    params.installCore(&vt, .x86_64, .skylakex, .{});
    try std.testing.expectEqual(@as(i32, 1), vt.dgemm_small_matrix_permit.?(0, 0, 10, 10, 10, 1, 0));
    try std.testing.expectEqual(@as(i32, 0), vt.zgemm_small_matrix_permit.?(0, 0, 1, 1, 1, 1, 0, 0, 0));
    try std.testing.expectEqual(@as(i32, 0), vt.sbgemm_small_matrix_permit.?(0, 0, 8, 8, 8, 1, 0));
    params.installCore(&vt, .x86_64, .haswell, .{});
    try std.testing.expectEqual(@as(i32, 0), vt.dgemm_small_matrix_permit.?(0, 0, 1, 1, 1, 1, 0));
}

fn always(comptime F: type) fn (i32, i32, BlasLong, BlasLong, BlasLong, F, F) callconv(.c) i32 {
    return struct {
        fn f(_: i32, _: i32, _: BlasLong, _: BlasLong, _: BlasLong, _: F, _: F) callconv(.c) i32 {
            return 1;
        }
    }.f;
}

fn alwaysComplex(comptime F: type) fn (i32, i32, BlasLong, BlasLong, BlasLong, F, F, F, F) callconv(.c) i32 {
    return struct {
        fn f(_: i32, _: i32, _: BlasLong, _: BlasLong, _: BlasLong, _: F, _: F, _: F, _: F) callconv(.c) i32 {
            return 1;
        }
    }.f;
}

/// `cases`' vtable: static, as a vtable is large for a rig's stack.
var permissive: Vtable = undefined;

/// `gemm`'s cases with every problem on the small path.
fn cases(vt: *const Vtable) !void {
    permissive = vt.*;
    permissive.sgemm_small_matrix_permit = &always(f32);
    permissive.dgemm_small_matrix_permit = &always(f64);
    permissive.cgemm_small_matrix_permit = &alwaysComplex(f32);
    permissive.zgemm_small_matrix_permit = &alwaysComplex(f64);
    try check3.routineAll(&permissive, .gemm);
}

test "portable small kernels match the wide reference through gemm" {
    var vt: Vtable = .{};
    portable.install(&vt);
    try cases(&vt);
}

test "every linked core's small kernels match the wide reference, where this machine runs them" {
    try check.eachCore(cases);
}

// ---- Skylake-X's kernels, on any x86_64 --------------------------------------------

/// `small_x86_64.zig`'s kernels, which LLVM splits into narrower vectors
/// on a machine without AVX-512: the same lanes, so the same bits.
pub const skylakex = @import("../../blas3/gemm/small_x86_64.zig");

/// Skylake-X's element `(i, j)` of an `m x n x k` problem of `op` (`"nn"`
/// ...), as upstream computes it: the depth in order, fused, except in the
/// dot blocks (`L` lanes, the masked tail, `REDUCE_4` or
/// `_mm512_reduce_add`); then the store. Fused as the kernel's `fma.madd`
/// does: not without an FMA instruction (a baseline x86_64 build, T154).
fn skylakexModel(comptime T: type, comptime op: []const u8, b0: bool, m: usize, n: usize, k: usize, a: []const T, ld: usize, b: []const T, alpha: T, beta: T, c: T, i: usize, j: usize) T {
    const L = skylakex.Skylakex(T).lanes;
    const ta = op[0] == 't';
    const tb = op[1] == 't';
    const at = struct {
        fn f(x: []const T, t: bool, ldx: usize, r: usize, q: usize) T {
            return if (t) x[q + r * ldx] else x[r + q * ldx];
        }
    }.f;
    // the dot rows: all of `tn`'s; `nn`'s past its vectors, for at most
    // half a vector of them over a depth of two vectors or more
    const mm = m % L;
    const start: ?usize = if (ta and !tb) 0 else if (!ta and !tb and mm > 0 and mm <= L / 2 and k >= 2 * L) m - mm else null;
    const fused = struct {
        fn f(x: T, y: T, z: T) T {
            return fma.madd(true, T, x, y, z);
        }
    }.f;
    if (start == null or i < start.?) {
        var acc: T = 0;
        for (0..k) |l| acc = fused(at(a, ta, ld, i, l), at(b, tb, ld, l, j), acc);
        const r = acc * alpha;
        return if (b0) r else fused(c, beta, r);
    }
    var lanes: [64]T = @splat(0);
    var l0: usize = 0;
    while (l0 < k) : (l0 += L) for (0..L) |x| {
        const l = l0 + x;
        lanes[x] = if (l < k) fused(at(a, ta, ld, i, l), at(b, tb, ld, l, j), lanes[x]) else fused(0, 0, lanes[x]);
    };
    const rows4 = (m - start.?) & ~@as(usize, 3);
    const grouped = i - start.? < rows4 or j < (n & ~@as(usize, 3));
    var w: usize = L;
    while (w > 1) : (w /= 2) for (0..w / 2) |x| {
        lanes[x] = if (grouped) lanes[2 * x] + lanes[2 * x + 1] else lanes[x] + lanes[x + w / 2];
    };
    if (grouped) {
        const s = lanes[0] * alpha;
        return if (b0) s else fused(c, beta, s);
    }
    return if (b0) alpha * lanes[0] else fma.madd(fma.fuse, T, alpha, lanes[0], beta * c);
}

/// Every Skylake-X kernel over sizes that reach each of its blocks: equal
/// bits to `skylakexModel` on uneven data, and to plain arithmetic on small
/// integers; `C`'s padding untouched.
fn checkSkylakex(comptime T: type) !void {
    const K = skylakex.Skylakex(T);
    const L = K.lanes;
    const al = std.testing.allocator;
    const ms = [_]usize{ 1, 3, 4, 7, L / 2 + 1, 3 * L / 4, L + 3, 3 * L + L / 2, 7 * L + 3 };
    const ns = [_]usize{ 1, 2, 3, 7, 13 };
    const ms_wide = [_]usize{ 1, 3, 4, 9, 3 * L / 4 - 1 };
    const ns_wide = [_]usize{ 6 * L + L / 2 + 1, 2 * L + 1 };
    const ks = [_]usize{ 1, 2 * L - 1, 2 * L, 2 * L + L / 2 + 1 };
    const ld = 7 * L + 3 + 2;
    const kmax = 2 * L + L / 2 + 1;
    const size = ld * @max(ld, kmax);
    const buf = try al.alloc(T, 5 * size);
    defer al.free(buf);
    const a = buf[0..size];
    const b = buf[size..][0..size];
    const c = buf[2 * size ..][0..size];
    const want = buf[3 * size ..][0..size];
    const got = buf[4 * size ..][0..size];
    var prng = std.Random.DefaultPrng.init(0x5ca1ab1e);
    const rnd = prng.random();
    inline for (.{ false, true }) |exact| {
        for ([_][]T{ a, b, c }) |x| for (x) |*v| {
            v.* = if (exact) @floatFromInt(rnd.intRangeAtMost(i32, -4, 4)) else std.math.ldexp(rnd.float(T) * 2 - 1, rnd.intRangeAtMost(i32, -6, 6));
        };
        const alpha: T = if (exact) 2 else 0.7;
        const beta: T = if (exact) -3 else 1.3;
        const fns = .{
            .{ "nn", K.gemm_small_kernel_nn, K.gemm_small_kernel_b0_nn },
            .{ "nt", K.gemm_small_kernel_nt, K.gemm_small_kernel_b0_nt },
            .{ "tn", K.gemm_small_kernel_tn, K.gemm_small_kernel_b0_tn },
            .{ "tt", K.gemm_small_kernel_tt, K.gemm_small_kernel_b0_tt },
        };
        inline for (fns) |e| inline for (.{ false, true }) |b0| {
            const op = e[0];
            const f = if (b0) e[2] else e[1];
            for ([_][]const usize{ &ms, &ms_wide }, [_][]const usize{ &ns, &ns_wide }) |mlist, nlist| for (mlist) |m| for (nlist) |n| for (ks) |k| {
                @memcpy(want, c);
                @memcpy(got, c);
                for (0..n) |j| for (0..m) |i| {
                    const p = i + j * ld;
                    if (exact) {
                        var s: T = 0;
                        for (0..k) |l| s += (if (op[0] == 't') a[l + i * ld] else a[i + l * ld]) * (if (op[1] == 't') b[j + l * ld] else b[l + j * ld]);
                        want[p] = alpha * s + (if (b0) 0 else beta * c[p]);
                    } else want[p] = skylakexModel(T, op, b0, m, n, k, a, ld, b, alpha, beta, c[p], i, j);
                };
                const lm: BlasLong = @intCast(ld);
                if (b0) _ = f(@intCast(m), @intCast(n), @intCast(k), a.ptr, lm, alpha, b.ptr, lm, got.ptr, lm) else _ = f(@intCast(m), @intCast(n), @intCast(k), a.ptr, lm, alpha, b.ptr, lm, beta, got.ptr, lm);
                if (!std.mem.eql(u8, std.mem.sliceAsBytes(want), std.mem.sliceAsBytes(got))) {
                    std.debug.print("{s} {s}{s} exact={} m={d} n={d} k={d}\n", .{ @typeName(T), if (b0) "b0_" else "", op, exact, m, n, k });
                    return error.TestUnexpectedResult;
                }
            };
        };
    }
}

test "Skylake-X small kernels: upstream's blocks, lanes and reductions, on any x86_64" {
    if (comptime sysinfo.target.cpu.arch != .x86_64) return error.SkipZigTest;
    try checkSkylakex(f32);
    try checkSkylakex(f64);
}

test "Skylake-X small kernels match the wide reference through gemm" {
    if (comptime sysinfo.target.cpu.arch != .x86_64) return error.SkipZigTest;
    var vt: Vtable = .{};
    portable.install(&vt);
    inline for (.{ .{ portable.Gemm.s, f32 }, .{ portable.Gemm.d, f64 } }) |p| {
        const K = skylakex.Kernels(p[1], gk.geometry(p[1]));
        const S = portable.GemmSlot;
        const kernels = .{
            .{ S.small_kernel_nn, &K.gemm_small_kernel_nn },
            .{ S.small_kernel_nt, &K.gemm_small_kernel_nt },
            .{ S.small_kernel_tn, &K.gemm_small_kernel_tn },
            .{ S.small_kernel_tt, &K.gemm_small_kernel_tt },
            .{ S.small_kernel_b0_nn, &K.gemm_small_kernel_b0_nn },
            .{ S.small_kernel_b0_nt, &K.gemm_small_kernel_b0_nt },
            .{ S.small_kernel_b0_tn, &K.gemm_small_kernel_b0_tn },
            .{ S.small_kernel_b0_tt, &K.gemm_small_kernel_b0_tt },
        };
        inline for (kernels) |e| portable.gemmSlotPtr(&vt, p[0], e[0]).* = e[1];
    }
    try cases(&vt);
}
