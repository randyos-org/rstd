//! Tests of `../../blas3/driver.zig`.

const std = @import("std");
const types = @import("../../types.zig");
const Vtable = @import("../../Vtable.zig");
const portable = @import("../../portable.zig");
const util = @import("../../util.zig");
const droot = @import("../../blas3/driver.zig");
const drout = @import("../../blas3/driver/routines.zig");

const Float = types.Float;
const isComplex = types.isComplex;
const Blocking = droot.Blocking;
const isZero = droot.isZero;
const scale = droot.scale;
const scaleBlock = drout.scaleBlock;

// ---- tests ---------------------------------------------------------------------

const check = @import("../check.zig");
const check3 = @import("check.zig");

const fromReal = util.fromReal;

const Routine = check3.Routine;

/// The routines `foldedCase` takes.
const routines = [_]Routine{ .gemm, .symm, .hemm, .syrk, .herk, .syr2k, .her2k, .trmm, .trsm };

/// The operands of one `foldedCase`.
fn Operands(comptime T: type) type {
    return struct { a: []T, b: []T, c: []T, ld: usize };
}

/// `r` with `s` as its `beta` (`trmm`'s and `trsm`'s `alpha`), options
/// `opt`, on `o`: `opt` is the side (`gemm`'s `transa`, the rank updates'
/// `upper`), then `upper` (`transb`, `trans`), then `trans` for `trmm` and
/// `trsm` (both conjugations for complex `gemm`; `trsm` conjugates complex).
fn callWith(comptime T: type, vt: *const Vtable, r: Routine, opt: [3]bool, s: T, m: usize, n: usize, k: usize, o: Operands(T)) !void {
    const complex = comptime isComplex(T);
    const alpha: T = if (complex) .init(-0.75, if (r == .herk) 0 else 0.25) else -0.75;
    const rank = r == .syrk or r == .herk or r == .syr2k or r == .her2k;
    const triangle = r == .trmm or r == .trsm;
    const conj = opt[2] and complex;
    const case: check3.Case = .{
        .r = r,
        .ta = if (triangle) .of(opt[2], r == .trsm and complex) else if (rank) .of(opt[1], false) else .of(opt[0], conj),
        .tb = .of(opt[1], conj),
        .left = opt[0],
        .upper = if (rank) opt[0] else opt[1],
        .m = m,
        .n = n,
        .k = k,
    };
    try check3.callAny(T, std.testing.allocator, check.io, vt, r, case, o.a.ptr, o.ld, o.b.ptr, o.ld, if (triangle) s else alpha, s, o.c.ptr, o.ld, 1);
}

/// `r` with `s` folded in against a separate pass of `s` first and then 1,
/// on a blocking shrunk so every edge comes at small sizes: the same bits.
fn foldedCase(comptime T: type, vt: *const Vtable, r: Routine, opt: [3]bool, s: T) !void {
    const al = std.testing.allocator;
    const blk = Blocking.of(T, vt);
    const m = 2 * blk.p + blk.um + 1;
    const n = 2 * blk.r + 3;
    const k = 2 * blk.q + 3;
    const ld = @max(m, n, k) + 1;
    const a = try al.alloc(T, ld * ld);
    defer al.free(a);
    const b = try al.alloc(T, ld * ld);
    defer al.free(b);
    const c = try al.alloc(T, ld * n);
    defer al.free(c);
    const c2 = try al.alloc(T, ld * n);
    defer al.free(c2);
    check.fill(T, a, 1);
    check.fill(T, b, 2);
    check.fill(T, c, 3);
    // a tame triangle for `trsm`: small off the diagonal, at least 1 on it
    const f: Float(T) = 1 / @as(Float(T), @floatFromInt(ld));
    for (0..ld) |j| for (0..ld) |i| {
        const p = &a[i + j * ld];
        const g: Float(T) = if (i == j) 1 else f;
        const h: Float(T) = if (i == j) 2 else 0;
        p.* = if (comptime isComplex(T)) .init(p.re * g + h, p.im * g) else p.* * g + h;
    };
    // a zero `s` must not read `C`
    if (isZero(T, s)) @memset(c, if (comptime isComplex(T)) .init(std.math.nan(Float(T)), 0) else std.math.nan(T));
    @memcpy(c2, c);
    const o: Operands(T) = .{ .a = a, .b = b, .c = c, .ld = ld };
    try callWith(T, vt, r, opt, s, m, n, k, o);
    const rank = r == .syrk or r == .herk or r == .syr2k or r == .her2k;
    if (rank) {
        scaleBlock(T, vt, if (opt[0]) .upper else .lower, r == .herk or r == .her2k, 0, n, 0, n, s, c2.ptr, ld);
    } else scale(T, vt, m, n, s, c2.ptr, ld);
    try callWith(T, vt, r, opt, util.one(T), m, n, k, .{ .a = a, .b = b, .c = c2, .ld = ld });
    try std.testing.expectEqualSlices(u8, std.mem.sliceAsBytes(c2), std.mem.sliceAsBytes(c));
}

test "beta folded into the first slice gives the bits of a separate pass first" {
    inline for (check.sd_types) |T| {
        var vt: Vtable = .{};
        portable.install(&vt);
        const own = Blocking.of(T, &vt);
        check3.shrinkBlocking(T, &vt, false, own.um, own.un);
        const complex = comptime isComplex(T);
        const ss = [_]T{ fromReal(T, 0), if (complex) .init(0.5, -0.25) else 0.5, fromReal(T, -2) };
        for (routines) |r| {
            if (!complex and (r == .hemm or r == .herk or r == .her2k)) continue;
            for (0..8) |bits| {
                const opt: [3]bool = .{ bits & 1 != 0, bits & 2 != 0, bits & 4 != 0 };
                for (ss) |s0| {
                    // `trmm`'s and `trsm`'s zero `alpha` returns before the drivers
                    if ((r == .trmm or r == .trsm) and isZero(T, s0)) continue;
                    // `herk`'s and `her2k`'s `beta` is real
                    const s = if (complex and (r == .herk or r == .her2k)) T.init(s0.re, 0) else s0;
                    foldedCase(T, &vt, r, opt, s) catch |e| {
                        std.debug.print("{s} {t} opt={any} s={any}\n", .{ @typeName(T), r, opt, s });
                        return e;
                    };
                }
            }
        }
    }
}
