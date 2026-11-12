//! `GEMV` kernels, `y += alpha * op(A) * x` on column-major `A`: the portable
//! ones (`*_generic`; the extended types' are `extended.zig`'s x87 sums) and
//! tests of every core's against plain arithmetic. The half-precision gemvs'
//! kernels are `ext/bf16/kernel.zig`'s `Gemv`, the same on every core.

const blas = @import("../../__root__.zig");
const types = @import("../../types.zig");
const portable = @import("../../portable.zig");
const std = @import("std");
const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");
const sysinfo = @import("builtin");
const bf16_kernel = @import("../../ext/bf16/kernel.zig");
const vector = @import("../../vector.zig");
const extended = @import("extended.zig");

const BlasLong = types.BlasLong;
const XDouble = types.XDouble;

// ---- portable kernels ------------------------------------------------------

/// `kernel/arm/gemv_n.c` and `gemv_t.c`: `y += alpha*A*x` (`n`, `y` of `m`)
/// or `y += alpha*A^T*x` (`t`, `y` of `n`). The interface has already scaled
/// `y` by `beta`. Wider than `f64`, the x87 `qgemv_n.S`/`qgemv_t.S` sums
/// instead (`extendedGemv`).
fn Gemv(comptime F: type) type {
    return struct {
        fn n(m: BlasLong, cols: BlasLong, _: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(types.kernelCc(F)) i32 {
            if (comptime types.isExtended(F)) {
                extendedGemv(F, false, extendedColBlock(false, .unknown), m, cols, alpha, a, lda, x, inc_x, y, inc_y);
                return 0;
            }
            var j: BlasLong = 0;
            while (j < cols) : (j += 1) {
                const tv = alpha * util.at(F, x, j * inc_x).*;
                var i: BlasLong = 0;
                while (i < m) : (i += 1) util.at(F, y, i * inc_y).* += tv * util.at(F, a, i + j * lda).*;
            }
            return 0;
        }
        fn t(m: BlasLong, cols: BlasLong, _: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(types.kernelCc(F)) i32 {
            if (comptime types.isExtended(F)) {
                extendedGemv(F, true, extendedRowBlock(false), m, cols, alpha, a, lda, x, inc_x, y, inc_y);
                return 0;
            }
            var j: BlasLong = 0;
            while (j < cols) : (j += 1) {
                const A = types.Acc(F);
                var acc: A = 0;
                var i: BlasLong = 0;
                while (i < m) : (i += 1) acc += @as(A, util.at(F, a, i + j * lda).*) * util.at(F, x, i * inc_x).*;
                util.at(F, y, j * inc_y).* += @floatCast(alpha * acc);
            }
            return 0;
        }
    };
}

/// `kernel/arm/zgemv_n.c` and `zgemv_t.c`: `y += alpha*op(A)*x` with `A`
/// conjugated for `conj` (`CONJ`) and `x` for `xconj` (`XCONJ`): `r` is
/// `conj(A)*x`, `o` is `A*conj(x)`, `c` is `A^H*x`, and so on. Wider than
/// `f64`, the x87 `xgemv_n.S`/`xgemv_t.S` sums instead (`extendedZgemv`).
fn Zgemv(comptime F: type, comptime trans: bool, comptime conj: bool, comptime xconj: bool) type {
    const W = std.math.Complex(types.Acc(F));
    return struct {
        fn wide(v: std.math.Complex(F)) W {
            return .init(v.re, v.im);
        }
        fn get(p: [*]F, i: BlasLong, c: bool) std.math.Complex(F) {
            const re = util.at(F, p, 2 * i).*;
            const im = util.at(F, p, 2 * i + 1).*;
            return .init(re, if (c) -im else im);
        }
        fn k(m: BlasLong, cols: BlasLong, _: BlasLong, alpha_r: F, alpha_i: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(types.kernelCc(F)) i32 {
            if (comptime types.isExtended(F)) {
                const p = comptime if (trans) extendedRowBlock(true) else extendedColBlock(true, .unknown);
                extendedZgemv(F, trans, conj, xconj, p, m, cols, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y);
                return 0;
            }
            const alpha: std.math.Complex(F) = .init(alpha_r, alpha_i);
            var j: BlasLong = 0;
            while (j < cols) : (j += 1) {
                if (trans) {
                    var acc: W = .init(0, 0);
                    var i: BlasLong = 0;
                    while (i < m) : (i += 1) acc = acc.add(wide(get(a, i + j * lda, conj)).mul(wide(get(x, i * inc_x, xconj))));
                    const v = wide(alpha).mul(acc);
                    util.at(F, y, 2 * j * inc_y).* += @floatCast(v.re);
                    util.at(F, y, 2 * j * inc_y + 1).* += @floatCast(v.im);
                } else {
                    const tv = alpha.mul(get(x, j * inc_x, xconj));
                    var i: BlasLong = 0;
                    while (i < m) : (i += 1) {
                        const v = tv.mul(get(a, i + j * lda, conj));
                        util.at(F, y, 2 * i * inc_y).* += v.re;
                        util.at(F, y, 2 * i * inc_y + 1).* += v.im;
                    }
                }
            }
            return 0;
        }
    };
}

// ---- extended types: `extended.zig` ----------------------------------------

const X86CoreType = extended.X86CoreType;
const extendedRowBlock = extended.extendedRowBlock;
const extendedColBlock = extended.extendedColBlock;
const extendedGemv = extended.extendedGemv;
const extendedZgemv = extended.extendedZgemv;
pub const installExtendedCore = extended.installExtendedCore;

pub fn Portable(comptime T: type) type {
    const F = types.Float(T);
    if (!types.isComplex(T)) return struct {
        pub const gemv_n = &Gemv(F).n;
        pub const gemv_t = &Gemv(F).t;
    };
    return struct {
        pub const gemv_n = &Zgemv(F, false, false, false).k;
        pub const gemv_t = &Zgemv(F, true, false, false).k;
        pub const gemv_r = &Zgemv(F, false, true, false).k;
        pub const gemv_c = &Zgemv(F, true, true, false).k;
        pub const gemv_o = &Zgemv(F, false, false, true).k;
        pub const gemv_u = &Zgemv(F, true, false, true).k;
        pub const gemv_s = &Zgemv(F, false, true, true).k;
        pub const gemv_d = &Zgemv(F, true, true, true).k;
    };
}

/// The complex variants by op: transposed, `A` conjugated, `x`
/// conjugated, as `types.Trans` says.
const variants = blk: {
    const all = std.enums.values(types.Trans);
    var v: [all.len]struct { portable.Op, bool, bool, bool } = undefined;
    for (all, &v) |t, *e| e.* = .{ portable.gemvOp(t), t.transposed(), t.conjA(), t.conjX() };
    break :blk v;
};

/// Points every gemv slot at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    bf16_kernel.installGenericGemv(vt);
    vt.sgemv_n = &Gemv(f32).n;
    vt.sgemv_t = &Gemv(f32).t;
    vt.dgemv_n = &Gemv(f64).n;
    vt.dgemv_t = &Gemv(f64).t;
    vt.qgemv_n = &Gemv(XDouble).n;
    vt.qgemv_t = &Gemv(XDouble).t;
    inline for (.{ .{ "c", f32 }, .{ "z", f64 }, .{ "x", XDouble } }) |p| {
        inline for (variants) |v| {
            portable.slotPtr(vt, .from(p[0]), v[0]).* = &Zgemv(p[1], v[1], v[2], v[3]).k;
        }
    }
}

// ---- tests -----------------------------------------------------------------

const check = @import("../../_test/check.zig");
const bitwise = @import("../../_test/bitwise.zig");
const guard = @import("../../_test/guard.zig");
const sym = @import("../symmetric.zig");
const zigstd = @import("../../../../zigstd.zig");
const fma = @import("../../fma.zig");

const C = std.math.Complex;
/// Shapes across the kernels' blocks, then bigger ones through their
/// multi-block loops (storage is static: `guard.operand`).
const sizes = [_][2]usize{
    .{ 1, 1 },   .{ 2, 3 },    .{ 3, 2 },    .{ 4, 4 },   .{ 7, 5 },   .{ 8, 8 }, .{ 9, 16 }, .{ 16, 9 }, .{ 17, 17 }, .{ 33, 4 }, .{ 4, 33 },
    .{ 65, 65 }, .{ 100, 37 }, .{ 37, 100 }, .{ 129, 3 }, .{ 3, 129 },
};
const strides = [_][2]usize{ .{ 1, 1 }, .{ 2, 3 } };

/// `kern` against plain arithmetic; `T` is the element type and the flags
/// are the variant's. Its buffer is the interface's, exactly (`sym.gemvLen`),
/// 32-byte aligned as a stack one is, with nothing written outside it.
fn checkGemv(comptime T: type, comptime trans: bool, comptime conj: bool, comptime xconj: bool, kern: anytype) !void {
    try checkGemvAt(T, trans, conj, xconj, kern, &sizes, &strides, null);
}

/// Past the x87 kernels' 4096-row blocks along the sum (`t` forms) and past
/// every column block (`n` forms, `DTB_ENTRIES` and `extendedColBlock`); on the
/// rigs only the shape that crosses the form's own block, at unit stride.
const big_sizes = [_][2]usize{ .{ 4100, 5 }, .{ 5, 4100 } };

fn checkGemvBig(comptime T: type, comptime trans: bool, comptime conj: bool, comptime xconj: bool, kern: anytype) !void {
    const rig = check.rig();
    const shapes: []const [2]usize = if (!rig) &big_sizes else if (trans) big_sizes[0..1] else big_sizes[1..2];
    const incs: []const [2]usize = if (rig) strides[0..1] else &strides;
    try checkGemvAt(T, trans, conj, xconj, kern, shapes, incs, std.testing.allocator);
}

/// `checkGemv` over `shapes` and `incs`, its operands from `al` where given,
/// else static.
fn checkGemvAt(comptime T: type, comptime trans: bool, comptime conj: bool, comptime xconj: bool, kern: anytype, shapes: []const [2]usize, incs: []const [2]usize, al: ?std.mem.Allocator) !void {
    const F = check.FloatOf(T);
    const complex = comptime check.isComplex(T);
    const cs = if (complex) 2 else 1;
    const alpha: T = if (complex) .init(1.25, -0.5) else 1.25;
    for (shapes) |mn| for (incs) |inc| {
        const m = mn[0];
        const n = mn[1];
        const lda = m + 1;
        const rows = if (trans) n else m; // y's length
        const inner = if (trans) m else n; // x's length
        const ops = try Operands(T).get(al, lda * n, guard.span(inner, inc[0]), guard.span(rows, inc[1]));
        defer ops.free(al);
        const a = ops.a;
        const x = ops.x;
        const y = ops.y;
        const want = ops.want;
        check.fill(T, a, m * 7 + n);
        check.fill(T, x, m + 3);
        check.fill(T, y, n + 5);
        @memcpy(want, y);
        for (0..rows) |r| {
            var acc: T = if (complex) .init(0, 0) else 0;
            for (0..inner) |q| {
                const e = if (trans) a[q + r * lda] else a[r + q * lda];
                const xv = x[q * inc[0]];
                if (complex) {
                    const ae = if (conj) e.conjugate() else e;
                    const xe = if (xconj) xv.conjugate() else xv;
                    acc = acc.add(ae.mul(xe));
                } else acc += e * xv;
            }
            if (complex) want[r * inc[1]] = want[r * inc[1]].add(alpha.mul(acc)) else want[r * inc[1]] += alpha * acc;
        }
        const bytes = sym.gemvLen(F, cs, m, n) * @sizeOf(F);
        // NaN: no kernel lets what the buffer held reach `y` (`buffers.md`)
        const buf = guard.open(F, bytes, 32, .nan);
        if (complex) {
            _ = kern(@intCast(m), @intCast(n), 0, alpha.re, alpha.im, @ptrCast(a.ptr), @intCast(lda), @ptrCast(x.ptr), @intCast(inc[0]), @ptrCast(y.ptr), @intCast(inc[1]), buf);
        } else {
            _ = kern(@intCast(m), @intCast(n), 0, alpha, a.ptr, @intCast(lda), x.ptr, @intCast(inc[0]), y.ptr, @intCast(inc[1]), buf);
        }
        check.expectCloseSlices(T, want, y, check.tolerance(F, 4 * inner)) catch |e| {
            check.say("m={d} n={d} trans={} conj={} xconj={}: ", .{ m, n, trans, conj, xconj });
            return check.case(e, "gemv", T, rows, inc);
        };
        guard.intact(bytes, 32) catch |e| {
            check.say("m={d} n={d} trans={} conj={} xconj={}: ", .{ m, n, trans, conj, xconj });
            return check.case(e, "gemv buffer", T, rows, inc);
        };
    };
}

fn checkAll(vt: *const Vtable) !void {
    try checkGemv(f32, false, false, false, vt.sgemv_n.?);
    try checkGemv(f32, true, false, false, vt.sgemv_t.?);
    try checkGemv(f64, false, false, false, vt.dgemv_n.?);
    try checkGemv(f64, true, false, false, vt.dgemv_t.?);
    inline for (.{ .{ "c", f32 }, .{ "z", f64 } }) |p| {
        inline for (variants) |v| {
            try checkGemv(C(p[1]), v[1], v[2], v[3], portable.slotPtr(vt, .from(p[0]), v[0]).*.?);
        }
    }
}

/// `checkGemv`'s operands: static, or from an allocator for big shapes.
fn Operands(comptime T: type) type {
    return struct {
        a: []T,
        x: []T,
        y: []T,
        want: []T,

        fn get(al: ?std.mem.Allocator, na: usize, nx: usize, ny: usize) !@This() {
            const g = al orelse return .{ .a = guard.operand(T, .a, na), .x = guard.operand(T, .x, nx), .y = guard.operand(T, .y, ny), .want = guard.operand(T, .want, ny) };
            const a = try g.alloc(T, na);
            errdefer g.free(a);
            const x = try g.alloc(T, nx);
            errdefer g.free(x);
            const y = try g.alloc(T, ny);
            errdefer g.free(y);
            return .{ .a = a, .x = x, .y = y, .want = try g.alloc(T, ny) };
        }

        fn free(self: @This(), al: ?std.mem.Allocator) void {
            const g = al orelse return;
            g.free(self.want);
            g.free(self.y);
            g.free(self.x);
            g.free(self.a);
        }
    };
}

/// `checkAll` past the blocks.
fn checkAllBig(vt: *const Vtable) !void {
    try checkGemvBig(f32, false, false, false, vt.sgemv_n.?);
    try checkGemvBig(f32, true, false, false, vt.sgemv_t.?);
    try checkGemvBig(f64, false, false, false, vt.dgemv_n.?);
    try checkGemvBig(f64, true, false, false, vt.dgemv_t.?);
    inline for (.{ .{ "c", f32 }, .{ "z", f64 } }) |p| {
        inline for (variants) |v| {
            try checkGemvBig(C(p[1]), v[1], v[2], v[3], portable.slotPtr(vt, .from(p[0]), v[0]).*.?);
        }
    }
}

test "portable gemv kernels match plain arithmetic" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
    try checkGemv(XDouble, false, false, false, vt.qgemv_n.?);
    try checkGemv(XDouble, true, false, false, vt.qgemv_t.?);
    inline for (variants) |v| try checkGemv(C(XDouble), v[1], v[2], v[3], portable.slotPtr(vt, .x, v[0]).*.?);
}

test "portable gemv kernels at the extended type match plain arithmetic" {
    const E = check.Ext;
    try checkGemv(E, false, false, false, &Gemv(E).n);
    try checkGemv(E, true, false, false, &Gemv(E).t);
    inline for (variants) |v| try checkGemv(C(E), v[1], v[2], v[3], &Zgemv(E, v[1], v[2], v[3]).k);
}

test "every linked core's gemv kernels match plain arithmetic, where this machine runs them" {
    try check.eachCore(checkAll);
}

test "portable gemv kernels past the x87 row block and every column block" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAllBig(&vt);
    try checkGemvBig(XDouble, false, false, false, vt.qgemv_n.?);
    try checkGemvBig(XDouble, true, false, false, vt.qgemv_t.?);
    inline for (variants) |v| try checkGemvBig(C(XDouble), v[1], v[2], v[3], portable.slotPtr(vt, .x, v[0]).*.?);
}

test "every linked core's gemv kernels past the x87 row block and every column block" {
    try check.eachCore(checkAllBig);
}

/// Room for `shifted`'s `A`, `x` and `y` at the largest shift, 64-byte
/// aligned so each shift is the same misalignment on every run.
var shift_mem: [3][(100 + 3) * 13 + 4]C(f64) align(64) = undefined;
const shift_shapes = [_][2]usize{ .{ 7, 5 }, .{ 33, 9 }, .{ 100, 13 } };

/// `kern` on the same `A`, `x` and `y`, each shifted by 0 to 3 floats, any
/// alignment a float may have: the same bits in `y` each time, as
/// `threading.md`'s reproducibility needs.
fn shifted(comptime T: type, comptime trans: bool, comptime conj: bool, comptime xconj: bool, kern: anytype) !void {
    const F = check.FloatOf(T);
    const complex = comptime check.isComplex(T);
    const w = if (complex) 2 else 1;
    const alpha: T = if (complex) .init(1.25, -0.5) else 1.25;
    var vals: [3][(100 + 3) * 13]T = undefined;
    var first: [100]T = undefined;
    for (shift_shapes) |mn| {
        const m = mn[0];
        const n = mn[1];
        const lda = m + 3;
        const rows = if (trans) n else m;
        const inner = if (trans) m else n;
        const lens = [3]usize{ lda * n, inner, rows };
        for (0..3) |k| check.fill(T, vals[k][0..lens[k]], m + 5 * k);
        for (0..64) |o| {
            var p: [3][*]F = undefined;
            for (0..3) |k| {
                p[k] = @as([*]F, @ptrCast(&shift_mem[k])) + (o >> @intCast(2 * k)) % 4;
                @memcpy(p[k][0 .. w * lens[k]], @as([*]const F, @ptrCast(&vals[k]))[0 .. w * lens[k]]);
            }
            const bytes = sym.gemvLen(F, w, m, n) * @sizeOf(F);
            const buf = guard.open(F, bytes, 32, .nan);
            if (complex) {
                _ = kern(@intCast(m), @intCast(n), 0, alpha.re, alpha.im, p[0], @intCast(lda), p[1], 1, p[2], 1, buf);
            } else {
                _ = kern(@intCast(m), @intCast(n), 0, alpha, p[0], @intCast(lda), p[1], 1, p[2], 1, buf);
            }
            const got = @as([*]align(1) T, @ptrCast(p[2]))[0..rows];
            if (o == 0) {
                @memcpy(first[0..rows], got);
                continue;
            }
            for (0..rows) |r| if (!bitwise.same(T, first[r], got[r], true)) {
                check.say("{s} m={d} n={d} trans={} conj={} xconj={}, A, x, y shifted by {d},{d},{d}: y[{d}] {any}, unshifted {any}\n", .{ @typeName(T), m, n, trans, conj, xconj, o % 4, (o >> 2) % 4, o >> 4, r, got[r], first[r] });
                return error.TestUnexpectedResult;
            };
        }
    }
}

/// Every gemv kernel, the first failure of each named.
fn checkShifted(vt: *const Vtable) !void {
    var failed: ?anyerror = null;
    shifted(f32, false, false, false, vt.sgemv_n.?) catch |e| {
        failed = e;
    };
    shifted(f32, true, false, false, vt.sgemv_t.?) catch |e| {
        failed = e;
    };
    shifted(f64, false, false, false, vt.dgemv_n.?) catch |e| {
        failed = e;
    };
    shifted(f64, true, false, false, vt.dgemv_t.?) catch |e| {
        failed = e;
    };
    inline for (.{ .{ "c", f32 }, .{ "z", f64 } }) |p| {
        inline for (variants) |v| {
            shifted(C(p[1]), v[1], v[2], v[3], portable.slotPtr(vt, .from(p[0]), v[0]).*.?) catch |e| {
                failed = e;
            };
        }
    }
    if (failed) |e| return e;
}

test "portable gemv kernels: the same bits whatever the address" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkShifted(&vt);
}

test "every linked core's gemv kernels: the same bits whatever the address" {
    try check.eachCore(checkShifted);
}

// the extended types' x87 sums
test "x86_64's dgemv_n vector block and SSE2 dgemv_t blocks give the reference's arithmetic, bit for bit" {
    // Zig's own x86_64 backend rejects these blocks' operands and VEX (the objects take LLVM)
    if (sysinfo.target.cpu.arch != .x86_64 or sysinfo.zig_backend != .stage2_llvm) return error.SkipZigTest;
    const driver = @import("driver.zig");
    const xk = @import("x86_64/kernel.zig");
    const own = @import("x86_64/uk/dgemv_t_4.zig");
    var a: [4 * 40]f64 = undefined;
    var x: [40]f64 = undefined;
    check.fill(f64, &a, 11);
    check.fill(f64, &x, 12);
    const ap = [4][*]f64{ a[0..].ptr, a[40..].ptr, a[80..].ptr, a[120..].ptr };
    const alpha = [1]f64{1.25};
    // Skylake-X's `dgemv_kernel_4x4` and `4x2` against a row-by-row model
    inline for (.{ false, true }) |wide| inline for (.{ false, true }) |fuse| {
        const B = xk.VectorBlockWith(.dgemv_n, wide, fuse);
        var n: usize = 4;
        while (n <= 40) : (n += 4) {
            inline for (.{ 4, 2 }) |cols| {
                var y: [40]f64 = undefined;
                check.fill(f64, &y, n);
                var want = y;
                for (0..n) |i| {
                    var sum = fma.madd(fuse, f64, ap[0][i], x[0], ap[1][i] * x[1]);
                    for (2..cols) |c| sum = fma.madd(fuse, f64, ap[c][i], x[c], sum);
                    want[i] = fma.madd(fuse, f64, sum, alpha[0], want[i]);
                }
                if (cols == 4) B.dgemv_kernel_4x4(@intCast(n), &ap, &x, &y, &alpha) else B.dgemv_kernel_4x2(@intCast(n), &ap, &x, &y, &alpha);
                try std.testing.expectEqualSlices(f64, &want, &y);
            }
        }
    };
    // the SSE2 forms of `dgemv_t_4.c`'s `kernel_4x2` and `4x1` against its assembly
    if (comptime !std.Target.x86.featureSetHas(sysinfo.target.cpu.features, .sse3)) return;
    const S = driver.Fallback(f64, fma.fuse);
    var n: usize = 2;
    while (n <= 40) : (n += 2) {
        var want: [2]f64 = undefined;
        var got: [2]f64 = undefined;
        own.dgemv_kernel_4x2(@intCast(n), ap[0], ap[1], &x, &want);
        S.t4x2(@intCast(n), ap[0], ap[1], &x, &got);
        try std.testing.expectEqualSlices(f64, &want, &got);
        own.dgemv_kernel_4x1(@intCast(n), ap[2], &x, &want);
        S.t4x1(@intCast(n), ap[2], &x, &got);
        try std.testing.expectEqual(want[0], got[0]);
    }
}

test "x86_64's Skylake-X sgemv (Haswell's blocks and its own) matches plain arithmetic, where this build has AVX2 and FMA" {
    if (sysinfo.target.cpu.arch != .x86_64 or sysinfo.zig_backend != .stage2_llvm) return error.SkipZigTest;
    const f = sysinfo.target.cpu.features;
    if (comptime !std.Target.x86.featureSetHasAll(f, .{ .avx2, .fma })) return error.SkipZigTest;
    const S = @import("x86_64/kernel.zig").Skylake;
    // `n`: the few-column path, the column-block path with its 2- and
    // 1-column tails (`n` past 48), and 4 columns, which takes the second
    try checkGemv(f32, false, false, false, S.entry(.sgemv_n));
    try checkGemvAt(f32, false, false, false, S.entry(.sgemv_n), &.{ .{ 37, 51 }, .{ 70, 55 }, .{ 20, 4 }, .{ 129, 49 } }, &strides, std.testing.allocator);
    // `t`: `checkGemv`'s `lda = m + 1` takes the column-block path, and
    // `lda == m` every shape of up to 8 rows across its blocks
    try checkGemv(f32, true, false, false, S.entry(.sgemv_t));
    try checkSkylakeT(S.entry(.sgemv_t));
}

/// `kern` with `lda == m` against plain arithmetic, for `m` up to 8.
fn checkSkylakeT(kern: anytype) !void {
    var a: [8 * 70]f32 = undefined;
    var x: [16]f32 = undefined;
    var y: [210]f32 = undefined;
    for (1..9) |m| for ([_]usize{ 1, 2, 3, 5, 7, 9, 15, 16, 17, 31, 33, 35, 47, 64, 69 }) |n| for (strides) |inc| {
        check.fill(f32, &a, m * 3 + n);
        check.fill(f32, &x, m);
        check.fill(f32, &y, n + 1);
        var want = y;
        for (0..n) |j| {
            var acc: f64 = 0;
            for (0..m) |i| acc += @as(f64, a[j * m + i]) * x[i * inc[0]];
            want[j * inc[1]] = @floatCast(@as(f64, want[j * inc[1]]) + 1.25 * acc);
        }
        var buf: [160]f32 = @splat(std.math.nan(f32));
        _ = kern(@intCast(m), @intCast(n), 0, 1.25, &a, @intCast(m), &x, @intCast(inc[0]), &y, @intCast(inc[1]), &buf);
        check.expectCloseSlices(f32, &want, &y, check.tolerance(f32, 4 * m + 4)) catch |e| {
            check.say("sgemv_t m={d} n={d}: ", .{ m, n });
            return e;
        };
    };
}

test {
    _ = extended;
}
