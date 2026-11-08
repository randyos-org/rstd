//! `trmm` and `trsm`: OpenBLAS's `driver/level3/trmm_L.c`, `trmm_R.c`,
//! `trsm_L.c` and `trsm_R.c` on `driver.zig`'s blocking, packing space and
//! kernel calls; more than one worker splits `B` into serial parts.

const std = @import("std");
const types = @import("../types.zig");
const Vtable = @import("../Vtable.zig");
const portable = @import("../portable.zig");
const util = @import("../util.zig");
const threading = @import("../threading.zig");
const driver = @import("driver.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const isComplex = types.isComplex;
const Op = portable.Op;
const Blocking = driver.Blocking;
const Conj = driver.Conj;
const Error = driver.Error;
const Work = driver.Work;
const kernel = driver.kernel;
const scale = driver.scale;
const unlessOne = driver.unlessOne;
const placed = driver.placed;
const rowsCopy = driver.rowsCopy;
const colsCopy = driver.colsCopy;
const stack_bytes = driver.stack_bytes;
const raw = util.raw;
const isZero = util.isZero;

/// A copy of the triangle (`trsm_iunucopy` and the like); `offset` places
/// the block's lines against its depth.
fn triCopy(comptime T: type, vt: *const Vtable, comptime op: Op, m: usize, n: usize, a: [*]T, lda: usize, offset: isize, dst: [*]T) void {
    if (comptime isComplex(T)) {
        _ = portable.kernel(T, vt, op)(@intCast(m), @intCast(n), raw(T, a), @intCast(lda), offset, raw(T, dst));
    } else {
        _ = portable.kernel(T, vt, op)(@intCast(m), @intCast(n), a, @intCast(lda), offset, dst);
    }
}

/// A `trmm` or `trsm` kernel: `alpha` for `trmm`; `trsm`'s, a dummy
/// upstream, is `unitTag`.
fn triKernel(comptime T: type, vt: *const Vtable, comptime op: Op, m: usize, n: usize, k: usize, alpha: T, sa: [*]T, sb: [*]T, c: [*]T, ldc: usize, offset: isize) void {
    if (comptime isComplex(T)) {
        _ = portable.kernel(T, vt, op)(@intCast(m), @intCast(n), @intCast(k), alpha.re, alpha.im, raw(T, sa), raw(T, sb), raw(T, c), @intCast(ldc), offset);
    } else {
        _ = portable.kernel(T, vt, op)(@intCast(m), @intCast(n), @intCast(k), alpha, sa, sb, c, @intCast(ldc), offset);
    }
}

const fromReal = util.fromReal;

/// The `trsm` kernels' `alpha`: 1 for a unit diagonal, which they skip
/// dividing by, else -1 (upstream's `dm1`).
fn unitTag(comptime T: type, unit: bool) T {
    return fromReal(T, if (unit) 1 else -1);
}

/// The options of a `trmm` or `trsm`: values, as `Desc`'s are.
pub const Tri = struct {
    left: bool,
    upper: bool,
    trans: bool,
    conj: bool,
    unit: bool,

    /// The gemm kernel's conjugation: the triangle's side.
    fn gemmConj(t: Tri) Conj {
        return if (!t.conj) .n else if (t.left) .l else .r;
    }
};

/// Which triangle copy a driver packs with: its triangle and whether it reads
/// across rows (`t`) or down columns (`n`).
const TriCopy = struct { upper: bool, across: bool };

/// `trmm_<io><u|l><n|t><u|n>copy` (`posX`, `posY`) or `trsm_...` (`offset`),
/// as `c` and `unit` say.
fn triPack(comptime T: type, vt: *const Vtable, comptime routine: []const u8, comptime io: []const u8, c: TriCopy, unit: bool, m: usize, n: usize, a: [*]T, lda: usize, x: usize, y: usize, offset: isize, dst: [*]T) void {
    switch (c.upper) {
        inline else => |u| switch (c.across) {
            inline else => |t| switch (unit) {
                inline else => |un| {
                    const op = @field(Op, routine ++ "_" ++ io ++ (if (u) "u" else "l") ++ (if (t) "t" else "n") ++ (if (un) "u" else "n") ++ "copy");
                    if (comptime std.mem.eql(u8, routine, "trmm")) return placed(T, vt, op, m, n, a, lda, x, y, dst);
                    return triCopy(T, vt, op, m, n, a, lda, offset, dst);
                },
            },
        },
    }
}

/// `<routine>_kernel_<L|R><N|T>`, or `R` and `C` for `N` and `T` conjugated.
fn triRun(comptime T: type, vt: *const Vtable, comptime routine: []const u8, o: Tri, nk: bool, m: usize, n: usize, k: usize, alpha: T, sa: [*]T, sb: [*]T, c: [*]T, ldc: usize, offset: isize) void {
    switch (o.left) {
        inline else => |l| switch (nk) {
            inline else => |nn| switch (o.conj) {
                inline else => |cj| if (comptime cj and !isComplex(T)) unreachable else {
                    const kk: []const u8 = if (cj) (if (nn) "R" else "C") else (if (nn) "N" else "T");
                    return triKernel(T, vt, @field(Op, routine ++ "_kernel_" ++ (if (l) "L" else "R") ++ kk), m, n, k, alpha, sa, sb, c, ldc, offset);
                },
            },
        },
    }
}

/// `op(A)` from row `i`, column `l`, as the copies take it.
fn opAt(comptime T: type, trans: bool, p: [*]T, ld: usize, i: usize, l: usize) [*]T {
    return if (trans) p + l + i * ld else p + i + l * ld;
}

/// The packing space of `trmm` and `trsm`, as `Work(T).lens` gives it: a row
/// block of `P` by a slice of `Q`, and that slice by `R` columns (the
/// reference's kernels take all of it).
pub fn triLens(blk: Blocking, m: usize, n: usize, ka: usize) [3]usize {
    const depth = @min(ka, blk.q);
    return .{ @min(m, blk.p) * depth, depth * @min(n, blk.r), 0 };
}

/// `B := alpha*op(A)*B` (`left`) or `alpha*B*op(A)`, `A` triangular
/// (`upper` or lower, `unit` diagonal), `op` transposing (`trans`) and
/// conjugating (`conj`, complex) as asked; `B` is `m x n`. `B` is scaled by
/// `alpha` first, never read for a zero one (`trmm_L.c`, `trmm_R.c`). Up to
/// `want` workers run it.
pub fn trmm(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, left: bool, upper: bool, trans: bool, conj: bool, unit: bool, m: usize, n: usize, alpha: T, a: [*]const T, lda: usize, b: [*]T, ldb: usize, want: u32) Error!void {
    return triangular(allocator, io, T, vt, false, left, upper, trans, conj, unit, m, n, alpha, a, lda, b, ldb, want);
}

/// `B := alpha*inv(op(A))*B` (`left`) or `alpha*B*inv(op(A))`, as `trmm`
/// takes its options (`trsm_L.c`, `trsm_R.c`). Up to `want` workers run it.
pub fn trsm(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, left: bool, upper: bool, trans: bool, conj: bool, unit: bool, m: usize, n: usize, alpha: T, a: [*]const T, lda: usize, b: [*]T, ldb: usize, want: u32) Error!void {
    return triangular(allocator, io, T, vt, true, left, upper, trans, conj, unit, m, n, alpha, a, lda, b, ldb, want);
}

fn triangular(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, comptime solve: bool, left: bool, upper: bool, trans: bool, conj: bool, unit: bool, m: usize, n: usize, alpha: T, a: [*]const T, lda: usize, b: [*]T, ldb: usize, want: u32) Error!void {
    // the reference passes `alpha` as `beta`; the drivers fold it in
    if (isZero(T, alpha)) return scale(T, vt, m, n, alpha, b, ldb);
    const s = unlessOne(T, alpha);
    const blk = Blocking.of(T, vt);
    const ap: [*]T = @constCast(a);
    const opt: Tri = .{ .left = left, .upper = upper, .trans = trans, .conj = isComplex(T) and conj, .unit = unit };
    if (want > 1) return parallel(allocator, io, T, vt, solve, opt, blk, m, n, ap, lda, s, b, ldb, want);
    var stack: [stack_bytes]u8 align(Work(T).alignment) = undefined;
    const lens = triLens(blk, m, n, if (left) m else n);
    const sc = try util.Scratch(T, Work(T).alignment, Work(T).alignment).init(allocator, &stack, Work(T).total(lens));
    defer sc.deinit(allocator);
    triDrive(T, vt, solve, opt, blk, m, n, ap, lda, s, b, ldb, Work(T).at(sc.ptr, lens));
}

/// The `trmm` (or with `solve`, `trsm`) loops on `work`.
pub fn triDrive(comptime T: type, vt: *const Vtable, comptime solve: bool, o: Tri, blk: Blocking, m: usize, n: usize, a: [*]T, lda: usize, s: ?T, b: [*]T, ldb: usize, work: Work(T)) void {
    if (solve) {
        if (o.left) trsmL(T, vt, o, blk, m, n, a, lda, s, b, ldb, work) else trsmR(T, vt, o, blk, m, n, a, lda, s, b, ldb, work);
    } else {
        if (o.left) trmmL(T, vt, o, blk, m, n, a, lda, s, b, ldb, work) else trmmR(T, vt, o, blk, m, n, a, lda, s, b, ldb, work);
    }
}

/// `B := s*B` on `m x n` at `b`, unless `s` is null: the triangular
/// drivers' `alpha`, where each part of `B` is first read.
fn scaleIf(comptime T: type, vt: *const Vtable, s: ?T, m: usize, n: usize, b: [*]T, ldb: usize) void {
    if (s) |v| scale(T, vt, m, n, v, b, ldb);
}

/// `min_i` of the `trmm` drivers: at most `P`, whole row panels past one.
fn trmmRows(blk: Blocking, rest: usize) usize {
    const v = @min(rest, blk.p);
    return if (v > blk.um) v / blk.um * blk.um else v;
}

fn trmmL(comptime T: type, vt: *const Vtable, o: Tri, blk: Blocking, m: usize, n: usize, a: [*]T, lda: usize, s: ?T, b: [*]T, ldb: usize, w: Work(T)) void {
    const p1 = util.one(T);
    const gc = o.gemmConj();
    // `op(A)` upper: the slices go down, each row taking the depths from it (`LN`)
    const down = o.upper != o.trans;
    const tc: TriCopy = .{ .upper = down != o.trans, .across = !o.trans };
    // each slice's copy of `B` is the first read of its rows: `alpha` goes there
    var js: usize = 0;
    while (js < n) : (js += blk.r) {
        const min_j = @min(n - js, blk.r);
        if (down) {
            var min_l = @min(m, blk.q);
            var min_i = trmmRows(blk, min_l);
            triPack(T, vt, "trmm", "i", tc, o.unit, min_l, min_i, a, lda, 0, 0, 0, w.sa);
            var jjs = js;
            while (jjs < js + min_j) {
                const min_jj = blk.chunk(js + min_j - jjs);
                const pb = w.sb + min_l * (jjs - js);
                scaleIf(T, vt, s, min_l, min_jj, b + jjs * ldb, ldb);
                colsCopy(T, vt, false, min_l, min_jj, b + jjs * ldb, ldb, pb);
                triRun(T, vt, "trmm", o, down, min_i, min_jj, min_l, p1, w.sa, pb, b + jjs * ldb, ldb, 0);
                jjs += min_jj;
            }
            var is = min_i;
            while (is < min_l) : (is += min_i) {
                min_i = trmmRows(blk, min_l - is);
                triPack(T, vt, "trmm", "i", tc, o.unit, min_l, min_i, a, lda, 0, is, 0, w.sa);
                triRun(T, vt, "trmm", o, down, min_i, min_j, min_l, p1, w.sa, w.sb, b + is + js * ldb, ldb, @intCast(is));
            }
            var ls = min_l;
            while (ls < m) : (ls += blk.q) {
                min_l = @min(m - ls, blk.q);
                min_i = trmmRows(blk, ls);
                rowsCopy(T, vt, o.trans, min_l, min_i, if (o.trans) a + ls else a + ls * lda, lda, w.sa);
                jjs = js;
                while (jjs < js + min_j) {
                    const min_jj = blk.chunk(js + min_j - jjs);
                    const pb = w.sb + min_l * (jjs - js);
                    scaleIf(T, vt, s, min_l, min_jj, b + ls + jjs * ldb, ldb);
                    colsCopy(T, vt, false, min_l, min_jj, b + ls + jjs * ldb, ldb, pb);
                    kernel(T, vt, gc, min_i, min_jj, min_l, p1, w.sa, pb, b + jjs * ldb, ldb);
                    jjs += min_jj;
                }
                is = min_i;
                while (is < ls) : (is += min_i) {
                    min_i = trmmRows(blk, ls - is);
                    rowsCopy(T, vt, o.trans, min_l, min_i, if (o.trans) a + ls + is * lda else a + is + ls * lda, lda, w.sa);
                    kernel(T, vt, gc, min_i, min_j, min_l, p1, w.sa, w.sb, b + is + js * ldb, ldb);
                }
                is = ls;
                while (is < ls + min_l) : (is += min_i) {
                    min_i = trmmRows(blk, ls + min_l - is);
                    triPack(T, vt, "trmm", "i", tc, o.unit, min_l, min_i, a, lda, ls, is, 0, w.sa);
                    triRun(T, vt, "trmm", o, down, min_i, min_j, min_l, p1, w.sa, w.sb, b + is + js * ldb, ldb, @intCast(is - ls));
                }
            }
        } else {
            // `op(A)` lower: the slices go up, each row taking the depths up to it (`LT`)
            var min_l = @min(m, blk.q);
            var min_i = trmmRows(blk, min_l);
            triPack(T, vt, "trmm", "i", tc, o.unit, min_l, min_i, a, lda, m - min_l, m - min_l, 0, w.sa);
            var jjs = js;
            while (jjs < js + min_j) {
                const min_jj = blk.chunk(js + min_j - jjs);
                const pb = w.sb + min_l * (jjs - js);
                scaleIf(T, vt, s, min_l, min_jj, b + (m - min_l) + jjs * ldb, ldb);
                colsCopy(T, vt, false, min_l, min_jj, b + (m - min_l) + jjs * ldb, ldb, pb);
                triRun(T, vt, "trmm", o, down, min_i, min_jj, min_l, p1, w.sa, pb, b + (m - min_l) + jjs * ldb, ldb, 0);
                jjs += min_jj;
            }
            var is = m - min_l + min_i;
            while (is < m) : (is += min_i) {
                min_i = trmmRows(blk, m - is);
                triPack(T, vt, "trmm", "i", tc, o.unit, min_l, min_i, a, lda, m - min_l, is, 0, w.sa);
                triRun(T, vt, "trmm", o, down, min_i, min_j, min_l, p1, w.sa, w.sb, b + is + js * ldb, ldb, @intCast(is - (m - min_l)));
            }
            var ls = m - min_l;
            while (ls > 0) : (ls -= min_l) {
                min_l = @min(ls, blk.q);
                const l0 = ls - min_l;
                min_i = trmmRows(blk, min_l);
                triPack(T, vt, "trmm", "i", tc, o.unit, min_l, min_i, a, lda, l0, l0, 0, w.sa);
                jjs = js;
                while (jjs < js + min_j) {
                    const min_jj = blk.chunk(js + min_j - jjs);
                    const pb = w.sb + min_l * (jjs - js);
                    scaleIf(T, vt, s, min_l, min_jj, b + l0 + jjs * ldb, ldb);
                    colsCopy(T, vt, false, min_l, min_jj, b + l0 + jjs * ldb, ldb, pb);
                    triRun(T, vt, "trmm", o, down, min_i, min_jj, min_l, p1, w.sa, pb, b + l0 + jjs * ldb, ldb, 0);
                    jjs += min_jj;
                }
                is = l0 + min_i;
                while (is < ls) : (is += min_i) {
                    min_i = trmmRows(blk, ls - is);
                    triPack(T, vt, "trmm", "i", tc, o.unit, min_l, min_i, a, lda, l0, is, 0, w.sa);
                    triRun(T, vt, "trmm", o, down, min_i, min_j, min_l, p1, w.sa, w.sb, b + is + js * ldb, ldb, @intCast(is - l0));
                }
                is = ls;
                while (is < m) : (is += min_i) {
                    min_i = trmmRows(blk, m - is);
                    rowsCopy(T, vt, o.trans, min_l, min_i, if (o.trans) a + l0 + is * lda else a + is + l0 * lda, lda, w.sa);
                    kernel(T, vt, gc, min_i, min_j, min_l, p1, w.sa, w.sb, b + is + js * ldb, ldb);
                }
            }
        }
    }
}

fn trmmR(comptime T: type, vt: *const Vtable, o: Tri, blk: Blocking, m: usize, n: usize, a: [*]T, lda: usize, s: ?T, b: [*]T, ldb: usize, w: Work(T)) void {
    const p1 = util.one(T);
    const gc = o.gemmConj();
    // `op(A)` lower: the columns go right, each taking the depths from it (`RT`)
    const right = o.upper == o.trans;
    const tc: TriCopy = .{ .upper = right == o.trans, .across = o.trans };
    const min_i0 = @min(m, blk.p);
    if (right) {
        var js: usize = 0;
        while (js < n) : (js += blk.r) {
            const min_j = @min(n - js, blk.r);
            // the first column block reads every column of `B` first
            const s0 = if (js == 0) s else null;
            var ls = js;
            while (ls < js + min_j) : (ls += blk.q) {
                const min_l = @min(js + min_j - ls, blk.q);
                scaleIf(T, vt, s0, min_i0, min_l, b + ls * ldb, ldb);
                rowsCopy(T, vt, false, min_l, min_i0, b + ls * ldb, ldb, w.sa);
                var jjs: usize = 0;
                while (jjs < ls - js) {
                    const min_jj = blk.chunk(ls - js - jjs);
                    colsCopy(T, vt, o.trans, min_l, min_jj, opAt(T, o.trans, a, lda, ls, js + jjs), lda, w.sb + min_l * jjs);
                    kernel(T, vt, gc, min_i0, min_jj, min_l, p1, w.sa, w.sb + min_l * jjs, b + (js + jjs) * ldb, ldb);
                    jjs += min_jj;
                }
                jjs = 0;
                while (jjs < min_l) {
                    const min_jj = blk.chunk(min_l - jjs);
                    const pb = w.sb + min_l * (ls - js + jjs);
                    triPack(T, vt, "trmm", "o", tc, o.unit, min_l, min_jj, a, lda, ls, ls + jjs, 0, pb);
                    triRun(T, vt, "trmm", o, !right, min_i0, min_jj, min_l, p1, w.sa, pb, b + (ls + jjs) * ldb, ldb, -@as(isize, @intCast(jjs)));
                    jjs += min_jj;
                }
                var is = min_i0;
                while (is < m) : (is += blk.p) {
                    const min_i = @min(m - is, blk.p);
                    scaleIf(T, vt, s0, min_i, min_l, b + is + ls * ldb, ldb);
                    rowsCopy(T, vt, false, min_l, min_i, b + is + ls * ldb, ldb, w.sa);
                    kernel(T, vt, gc, min_i, ls - js, min_l, p1, w.sa, w.sb, b + is + js * ldb, ldb);
                    triRun(T, vt, "trmm", o, !right, min_i, min_l, min_l, p1, w.sa, w.sb + (ls - js) * min_l, b + is + ls * ldb, ldb, 0);
                }
            }
            ls = js + min_j;
            while (ls < n) : (ls += blk.q) {
                const min_l = @min(n - ls, blk.q);
                scaleIf(T, vt, s0, min_i0, min_l, b + ls * ldb, ldb);
                rowsCopy(T, vt, false, min_l, min_i0, b + ls * ldb, ldb, w.sa);
                var jjs = js;
                while (jjs < js + min_j) {
                    const min_jj = blk.chunk(js + min_j - jjs);
                    const pb = w.sb + min_l * (jjs - js);
                    colsCopy(T, vt, o.trans, min_l, min_jj, opAt(T, o.trans, a, lda, ls, jjs), lda, pb);
                    kernel(T, vt, gc, min_i0, min_jj, min_l, p1, w.sa, pb, b + jjs * ldb, ldb);
                    jjs += min_jj;
                }
                var is = min_i0;
                while (is < m) : (is += blk.p) {
                    const min_i = @min(m - is, blk.p);
                    scaleIf(T, vt, s0, min_i, min_l, b + is + ls * ldb, ldb);
                    rowsCopy(T, vt, false, min_l, min_i, b + is + ls * ldb, ldb, w.sa);
                    kernel(T, vt, gc, min_i, min_j, min_l, p1, w.sa, w.sb, b + is + js * ldb, ldb);
                }
            }
        }
    } else {
        // `op(A)` upper: the columns go left, each taking the depths up to it (`RN`)
        var js = n;
        while (js > 0) : (js -= @min(js, blk.r)) {
            const min_j = @min(js, blk.r);
            const j0 = js - min_j;
            const s0 = if (js == n) s else null;
            var start_ls = j0;
            while (start_ls + blk.q < js) start_ls += blk.q;
            var ls = start_ls;
            while (true) {
                const min_l = @min(js - ls, blk.q);
                scaleIf(T, vt, s0, min_i0, min_l, b + ls * ldb, ldb);
                rowsCopy(T, vt, false, min_l, min_i0, b + ls * ldb, ldb, w.sa);
                var jjs: usize = 0;
                while (jjs < min_l) {
                    const min_jj = blk.chunk(min_l - jjs);
                    triPack(T, vt, "trmm", "o", tc, o.unit, min_l, min_jj, a, lda, ls, ls + jjs, 0, w.sb + min_l * jjs);
                    triRun(T, vt, "trmm", o, !right, min_i0, min_jj, min_l, p1, w.sa, w.sb + min_l * jjs, b + (ls + jjs) * ldb, ldb, -@as(isize, @intCast(jjs)));
                    jjs += min_jj;
                }
                jjs = 0;
                while (jjs < js - ls - min_l) {
                    const min_jj = blk.chunk(js - ls - min_l - jjs);
                    const pb = w.sb + min_l * (min_l + jjs);
                    colsCopy(T, vt, o.trans, min_l, min_jj, opAt(T, o.trans, a, lda, ls, ls + min_l + jjs), lda, pb);
                    kernel(T, vt, gc, min_i0, min_jj, min_l, p1, w.sa, pb, b + (ls + min_l + jjs) * ldb, ldb);
                    jjs += min_jj;
                }
                var is = min_i0;
                while (is < m) : (is += blk.p) {
                    const min_i = @min(m - is, blk.p);
                    scaleIf(T, vt, s0, min_i, min_l, b + is + ls * ldb, ldb);
                    rowsCopy(T, vt, false, min_l, min_i, b + is + ls * ldb, ldb, w.sa);
                    triRun(T, vt, "trmm", o, !right, min_i, min_l, min_l, p1, w.sa, w.sb, b + is + ls * ldb, ldb, 0);
                    if (js - ls - min_l > 0) kernel(T, vt, gc, min_i, js - ls - min_l, min_l, p1, w.sa, w.sb + min_l * min_l, b + is + (ls + min_l) * ldb, ldb);
                }
                if (ls < j0 + blk.q) break;
                ls -= blk.q;
            }
            ls = 0;
            while (ls < j0) : (ls += blk.q) {
                const min_l = @min(j0 - ls, blk.q);
                scaleIf(T, vt, s0, min_i0, min_l, b + ls * ldb, ldb);
                rowsCopy(T, vt, false, min_l, min_i0, b + ls * ldb, ldb, w.sa);
                var jjs = j0;
                while (jjs < js) {
                    const min_jj = blk.chunk(js - jjs);
                    const pb = w.sb + min_l * (jjs - j0);
                    colsCopy(T, vt, o.trans, min_l, min_jj, opAt(T, o.trans, a, lda, ls, jjs), lda, pb);
                    kernel(T, vt, gc, min_i0, min_jj, min_l, p1, w.sa, pb, b + jjs * ldb, ldb);
                    jjs += min_jj;
                }
                var is = min_i0;
                while (is < m) : (is += blk.p) {
                    const min_i = @min(m - is, blk.p);
                    scaleIf(T, vt, s0, min_i, min_l, b + is + ls * ldb, ldb);
                    rowsCopy(T, vt, false, min_l, min_i, b + is + ls * ldb, ldb, w.sa);
                    kernel(T, vt, gc, min_i, min_j, min_l, p1, w.sa, w.sb, b + is + j0 * ldb, ldb);
                }
            }
        }
    }
}

fn trsmL(comptime T: type, vt: *const Vtable, o: Tri, blk: Blocking, m: usize, n: usize, a: [*]T, lda: usize, s: ?T, b: [*]T, ldb: usize, w: Work(T)) void {
    const m1 = fromReal(T, -1);
    const du = unitTag(T, o.unit);
    const gc = o.gemmConj();
    // `op(A)` lower: forward from the top (`LT`), else back from the bottom (`LN`)
    const forward = o.upper == o.trans;
    const tc: TriCopy = .{ .upper = forward == o.trans, .across = !o.trans };
    var js: usize = 0;
    while (js < n) : (js += blk.r) {
        const min_j = @min(n - js, blk.r);
        if (forward) {
            var ls: usize = 0;
            while (ls < m) : (ls += blk.q) {
                const min_l = @min(m - ls, blk.q);
                // the first slice first reads each part of `B`: `alpha` goes there
                const s0 = if (ls == 0) s else null;
                var min_i = @min(min_l, blk.p);
                triPack(T, vt, "trsm", "i", tc, o.unit, min_l, min_i, opAt(T, o.trans, a, lda, ls, ls), lda, 0, 0, 0, w.sa);
                var jjs = js;
                while (jjs < js + min_j) {
                    const min_jj = blk.chunk(js + min_j - jjs);
                    const pb = w.sb + min_l * (jjs - js);
                    scaleIf(T, vt, s0, min_l, min_jj, b + ls + jjs * ldb, ldb);
                    colsCopy(T, vt, false, min_l, min_jj, b + ls + jjs * ldb, ldb, pb);
                    triRun(T, vt, "trsm", o, !forward, min_i, min_jj, min_l, du, w.sa, pb, b + ls + jjs * ldb, ldb, 0);
                    jjs += min_jj;
                }
                var is = ls + min_i;
                while (is < ls + min_l) : (is += blk.p) {
                    min_i = @min(ls + min_l - is, blk.p);
                    triPack(T, vt, "trsm", "i", tc, o.unit, min_l, min_i, opAt(T, o.trans, a, lda, is, ls), lda, 0, 0, @intCast(is - ls), w.sa);
                    triRun(T, vt, "trsm", o, !forward, min_i, min_j, min_l, du, w.sa, w.sb, b + is + js * ldb, ldb, @intCast(is - ls));
                }
                is = ls + min_l;
                while (is < m) : (is += blk.p) {
                    min_i = @min(m - is, blk.p);
                    rowsCopy(T, vt, o.trans, min_l, min_i, opAt(T, o.trans, a, lda, is, ls), lda, w.sa);
                    scaleIf(T, vt, s0, min_i, min_j, b + is + js * ldb, ldb);
                    kernel(T, vt, gc, min_i, min_j, min_l, m1, w.sa, w.sb, b + is + js * ldb, ldb);
                }
            }
        } else {
            var ls = m;
            while (ls > 0) : (ls -= @min(ls, blk.q)) {
                const min_l = @min(ls, blk.q);
                const l0 = ls - min_l;
                const s0 = if (ls == m) s else null;
                var start_is = l0;
                while (start_is + blk.p < ls) start_is += blk.p;
                var min_i = @min(ls - start_is, blk.p);
                triPack(T, vt, "trsm", "i", tc, o.unit, min_l, min_i, opAt(T, o.trans, a, lda, start_is, l0), lda, 0, 0, @intCast(start_is - l0), w.sa);
                var jjs = js;
                while (jjs < js + min_j) {
                    const min_jj = blk.chunk(js + min_j - jjs);
                    const pb = w.sb + min_l * (jjs - js);
                    scaleIf(T, vt, s0, min_l, min_jj, b + l0 + jjs * ldb, ldb);
                    colsCopy(T, vt, false, min_l, min_jj, b + l0 + jjs * ldb, ldb, pb);
                    triRun(T, vt, "trsm", o, !forward, min_i, min_jj, min_l, du, w.sa, pb, b + start_is + jjs * ldb, ldb, @intCast(start_is - l0));
                    jjs += min_jj;
                }
                var is = start_is;
                while (is >= l0 + blk.p) {
                    is -= blk.p;
                    min_i = @min(ls - is, blk.p);
                    triPack(T, vt, "trsm", "i", tc, o.unit, min_l, min_i, opAt(T, o.trans, a, lda, is, l0), lda, 0, 0, @intCast(is - l0), w.sa);
                    triRun(T, vt, "trsm", o, !forward, min_i, min_j, min_l, du, w.sa, w.sb, b + is + js * ldb, ldb, @intCast(is - l0));
                }
                is = 0;
                while (is < l0) : (is += blk.p) {
                    min_i = @min(l0 - is, blk.p);
                    rowsCopy(T, vt, o.trans, min_l, min_i, opAt(T, o.trans, a, lda, is, l0), lda, w.sa);
                    scaleIf(T, vt, s0, min_i, min_j, b + is + js * ldb, ldb);
                    kernel(T, vt, gc, min_i, min_j, min_l, m1, w.sa, w.sb, b + is + js * ldb, ldb);
                }
            }
        }
    }
}

fn trsmR(comptime T: type, vt: *const Vtable, o: Tri, blk: Blocking, m: usize, n: usize, a: [*]T, lda: usize, s: ?T, b: [*]T, ldb: usize, w: Work(T)) void {
    const m1 = fromReal(T, -1);
    const du = unitTag(T, o.unit);
    const gc = o.gemmConj();
    // `op(A)` upper: forward from the left (`RN`), else back from the right (`RT`)
    const forward = o.upper != o.trans;
    const tc: TriCopy = .{ .upper = forward != o.trans, .across = o.trans };
    const min_i0 = @min(m, blk.p);
    if (forward) {
        var js: usize = 0;
        while (js < n) : (js += blk.r) {
            const min_j = @min(n - js, blk.r);
            var ls: usize = 0;
            while (ls < js) : (ls += blk.q) {
                const min_l = @min(js - ls, blk.q);
                // as `trsmL`'s, the first slice of each column block
                const s0 = if (ls == 0) s else null;
                rowsCopy(T, vt, false, min_l, min_i0, b + ls * ldb, ldb, w.sa);
                var jjs = js;
                while (jjs < js + min_j) {
                    const min_jj = blk.chunk(js + min_j - jjs);
                    const pb = w.sb + min_l * (jjs - js);
                    colsCopy(T, vt, o.trans, min_l, min_jj, opAt(T, o.trans, a, lda, ls, jjs), lda, pb);
                    scaleIf(T, vt, s0, min_i0, min_jj, b + jjs * ldb, ldb);
                    kernel(T, vt, gc, min_i0, min_jj, min_l, m1, w.sa, pb, b + jjs * ldb, ldb);
                    jjs += min_jj;
                }
                var is = min_i0;
                while (is < m) : (is += blk.p) {
                    const min_i = @min(m - is, blk.p);
                    rowsCopy(T, vt, false, min_l, min_i, b + is + ls * ldb, ldb, w.sa);
                    scaleIf(T, vt, s0, min_i, min_j, b + is + js * ldb, ldb);
                    kernel(T, vt, gc, min_i, min_j, min_l, m1, w.sa, w.sb, b + is + js * ldb, ldb);
                }
            }
            ls = js;
            while (ls < js + min_j) : (ls += blk.q) {
                const min_l = @min(js + min_j - ls, blk.q);
                const s0 = if (ls == 0) s else null;
                scaleIf(T, vt, s0, min_i0, min_l, b + ls * ldb, ldb);
                rowsCopy(T, vt, false, min_l, min_i0, b + ls * ldb, ldb, w.sa);
                triPack(T, vt, "trsm", "o", tc, o.unit, min_l, min_l, opAt(T, o.trans, a, lda, ls, ls), lda, 0, 0, 0, w.sb);
                triRun(T, vt, "trsm", o, forward, min_i0, min_l, min_l, du, w.sa, w.sb, b + ls * ldb, ldb, 0);
                var jjs: usize = 0;
                while (jjs < js + min_j - ls - min_l) {
                    const min_jj = blk.chunk(js + min_j - ls - min_l - jjs);
                    const pb = w.sb + min_l * (min_l + jjs);
                    colsCopy(T, vt, o.trans, min_l, min_jj, opAt(T, o.trans, a, lda, ls, ls + min_l + jjs), lda, pb);
                    scaleIf(T, vt, s0, min_i0, min_jj, b + (ls + min_l + jjs) * ldb, ldb);
                    kernel(T, vt, gc, min_i0, min_jj, min_l, m1, w.sa, pb, b + (ls + min_l + jjs) * ldb, ldb);
                    jjs += min_jj;
                }
                var is = min_i0;
                while (is < m) : (is += blk.p) {
                    const min_i = @min(m - is, blk.p);
                    scaleIf(T, vt, s0, min_i, min_l, b + is + ls * ldb, ldb);
                    rowsCopy(T, vt, false, min_l, min_i, b + is + ls * ldb, ldb, w.sa);
                    triRun(T, vt, "trsm", o, forward, min_i, min_l, min_l, du, w.sa, w.sb, b + is + ls * ldb, ldb, 0);
                    scaleIf(T, vt, s0, min_i, js + min_j - ls - min_l, b + is + (ls + min_l) * ldb, ldb);
                    kernel(T, vt, gc, min_i, js + min_j - ls - min_l, min_l, m1, w.sa, w.sb + min_l * min_l, b + is + (ls + min_l) * ldb, ldb);
                }
            }
        }
    } else {
        var js = n;
        while (js > 0) : (js -= @min(js, blk.r)) {
            const min_j = @min(js, blk.r);
            const j0 = js - min_j;
            var ls = js;
            while (ls < n) : (ls += blk.q) {
                const min_l = @min(n - ls, blk.q);
                const s0 = if (ls == js) s else null;
                rowsCopy(T, vt, false, min_l, min_i0, b + ls * ldb, ldb, w.sa);
                var jjs = j0;
                while (jjs < js) {
                    const min_jj = blk.chunk(js - jjs);
                    const pb = w.sb + min_l * (jjs - j0);
                    colsCopy(T, vt, o.trans, min_l, min_jj, opAt(T, o.trans, a, lda, ls, jjs), lda, pb);
                    scaleIf(T, vt, s0, min_i0, min_jj, b + jjs * ldb, ldb);
                    kernel(T, vt, gc, min_i0, min_jj, min_l, m1, w.sa, pb, b + jjs * ldb, ldb);
                    jjs += min_jj;
                }
                var is = min_i0;
                while (is < m) : (is += blk.p) {
                    const min_i = @min(m - is, blk.p);
                    rowsCopy(T, vt, false, min_l, min_i, b + is + ls * ldb, ldb, w.sa);
                    scaleIf(T, vt, s0, min_i, min_j, b + is + j0 * ldb, ldb);
                    kernel(T, vt, gc, min_i, min_j, min_l, m1, w.sa, w.sb, b + is + j0 * ldb, ldb);
                }
            }
            var start_ls = j0;
            while (start_ls + blk.q < js) start_ls += blk.q;
            ls = start_ls;
            while (true) {
                const min_l = @min(js - ls, blk.q);
                const tb = w.sb + min_l * (ls - j0);
                const s0 = if (js == n and ls == start_ls) s else null;
                scaleIf(T, vt, s0, min_i0, min_l, b + ls * ldb, ldb);
                rowsCopy(T, vt, false, min_l, min_i0, b + ls * ldb, ldb, w.sa);
                triPack(T, vt, "trsm", "o", tc, o.unit, min_l, min_l, opAt(T, o.trans, a, lda, ls, ls), lda, 0, 0, 0, tb);
                triRun(T, vt, "trsm", o, forward, min_i0, min_l, min_l, du, w.sa, tb, b + ls * ldb, ldb, 0);
                var jjs: usize = 0;
                while (jjs < ls - j0) {
                    const min_jj = blk.chunk(ls - j0 - jjs);
                    const pb = w.sb + min_l * jjs;
                    colsCopy(T, vt, o.trans, min_l, min_jj, opAt(T, o.trans, a, lda, ls, j0 + jjs), lda, pb);
                    scaleIf(T, vt, s0, min_i0, min_jj, b + (j0 + jjs) * ldb, ldb);
                    kernel(T, vt, gc, min_i0, min_jj, min_l, m1, w.sa, pb, b + (j0 + jjs) * ldb, ldb);
                    jjs += min_jj;
                }
                var is = min_i0;
                while (is < m) : (is += blk.p) {
                    const min_i = @min(m - is, blk.p);
                    scaleIf(T, vt, s0, min_i, min_l, b + is + ls * ldb, ldb);
                    rowsCopy(T, vt, false, min_l, min_i, b + is + ls * ldb, ldb, w.sa);
                    triRun(T, vt, "trsm", o, forward, min_i, min_l, min_l, du, w.sa, tb, b + is + ls * ldb, ldb, 0);
                    scaleIf(T, vt, s0, min_i, ls - j0, b + is + j0 * ldb, ldb);
                    kernel(T, vt, gc, min_i, ls - j0, min_l, m1, w.sa, w.sb, b + is + j0 * ldb, ldb);
                }
                if (ls < j0 + blk.q) break;
                ls -= blk.q;
            }
        }
    }
}

// ---- threaded --------------------------------------------------------------------

/// `trmm` or, with `solve`, `trsm` on up to `want` workers: column
/// panels of `B` on the left, row panels on the right, as OpenBLAS's
/// `gemm_thread_n` and `gemm_thread_m` divide it and BLIS splits its `jc` and
/// `ic` loops. Each part is an independent serial problem on packing space
/// of its own, so they run through `threading.forEach`; parts start on
/// `unroll_n` columns or `unroll_m` rows, so tiles are the serial ones.
fn parallel(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, comptime solve: bool, o: Tri, blk: Blocking, m: usize, n: usize, a: [*]T, lda: usize, s: ?T, b: [*]T, ldb: usize, want: u32) Error!void {
    const W = Work(T);
    const Parts = struct {
        vt: *const Vtable,
        o: Tri,
        blk: Blocking,
        m: usize,
        n: usize,
        a: [*]T,
        lda: usize,
        s: ?T,
        b: [*]T,
        ldb: usize,

        /// Part `c`'s `m`, `n` and `B`.
        fn of(p: *const @This(), c: threading.Chunk) struct { usize, usize, [*]T } {
            return if (p.o.left) .{ p.m, c.end - c.start, p.b + c.start * p.ldb } else .{ c.end - c.start, p.n, p.b + c.start };
        }

        fn lens(p: *const @This(), c: threading.Chunk) [3]usize {
            const d = p.of(c);
            return triLens(p.blk, d[0], d[1], if (p.o.left) p.m else p.n);
        }

        fn len(p: *const @This(), c: threading.Chunk) usize {
            return W.total(p.lens(c));
        }

        fn one(p: *const @This(), c: threading.Chunk, buf: []T) void {
            const d = p.of(c);
            triDrive(T, p.vt, solve, p.o, p.blk, d[0], d[1], p.a, p.lda, p.s, d[2], p.ldb, W.at(buf.ptr, p.lens(c)));
        }
    };
    const plan: threading.Plan = .{
        .chunks = .init(if (o.left) n else m, want, .even, if (o.left) blk.un else blk.um),
        .workers = want,
    };
    const parts: Parts = .{ .vt = vt, .o = o, .blk = blk, .m = m, .n = n, .a = a, .lda = lda, .s = s, .b = b, .ldb = ldb };
    try threading.forEachPart(allocator, io, &plan, T, W.alignment, &parts, Parts.len, Parts.one);
}
