//! `trsm`'s drivers: `driver/level3/trsm_L.c` and `trsm_R.c` on the shared blocking and packing.

const Vtable = @import("../../Vtable.zig");
const driver = @import("../driver.zig");
const troot = @import("../triangular.zig");

const Blocking = driver.Blocking;
const Work = driver.Work;
const kernel = driver.kernel;
const rowsCopy = driver.rowsCopy;
const colsCopy = driver.colsCopy;
const Tri = troot.Tri;
const TriCopy = troot.TriCopy;
const fromReal = troot.fromReal;
const opAt = troot.opAt;
const scaleIf = troot.scaleIf;
const triPack = troot.triPack;
const triRun = troot.triRun;
const unitTag = troot.unitTag;

pub fn trsmL(comptime T: type, vt: *const Vtable, o: Tri, blk: Blocking, m: usize, n: usize, a: [*]T, lda: usize, s: ?T, b: [*]T, ldb: usize, w: Work(T)) void {
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
                triPack(T, vt, .trsm, .i, tc, o.unit, min_l, min_i, opAt(T, o.trans, a, lda, ls, ls), lda, 0, 0, 0, w.sa);
                var jjs = js;
                while (jjs < js + min_j) {
                    const min_jj = blk.chunk(js + min_j - jjs);
                    const pb = w.sb + min_l * (jjs - js);
                    scaleIf(T, vt, s0, min_l, min_jj, b + ls + jjs * ldb, ldb);
                    colsCopy(T, vt, false, min_l, min_jj, b + ls + jjs * ldb, ldb, pb);
                    triRun(T, vt, .trsm, o, !forward, min_i, min_jj, min_l, du, w.sa, pb, b + ls + jjs * ldb, ldb, 0);
                    jjs += min_jj;
                }
                var is = ls + min_i;
                while (is < ls + min_l) : (is += blk.p) {
                    min_i = @min(ls + min_l - is, blk.p);
                    triPack(T, vt, .trsm, .i, tc, o.unit, min_l, min_i, opAt(T, o.trans, a, lda, is, ls), lda, 0, 0, @intCast(is - ls), w.sa);
                    triRun(T, vt, .trsm, o, !forward, min_i, min_j, min_l, du, w.sa, w.sb, b + is + js * ldb, ldb, @intCast(is - ls));
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
                triPack(T, vt, .trsm, .i, tc, o.unit, min_l, min_i, opAt(T, o.trans, a, lda, start_is, l0), lda, 0, 0, @intCast(start_is - l0), w.sa);
                var jjs = js;
                while (jjs < js + min_j) {
                    const min_jj = blk.chunk(js + min_j - jjs);
                    const pb = w.sb + min_l * (jjs - js);
                    scaleIf(T, vt, s0, min_l, min_jj, b + l0 + jjs * ldb, ldb);
                    colsCopy(T, vt, false, min_l, min_jj, b + l0 + jjs * ldb, ldb, pb);
                    triRun(T, vt, .trsm, o, !forward, min_i, min_jj, min_l, du, w.sa, pb, b + start_is + jjs * ldb, ldb, @intCast(start_is - l0));
                    jjs += min_jj;
                }
                var is = start_is;
                while (is >= l0 + blk.p) {
                    is -= blk.p;
                    min_i = @min(ls - is, blk.p);
                    triPack(T, vt, .trsm, .i, tc, o.unit, min_l, min_i, opAt(T, o.trans, a, lda, is, l0), lda, 0, 0, @intCast(is - l0), w.sa);
                    triRun(T, vt, .trsm, o, !forward, min_i, min_j, min_l, du, w.sa, w.sb, b + is + js * ldb, ldb, @intCast(is - l0));
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

pub fn trsmR(comptime T: type, vt: *const Vtable, o: Tri, blk: Blocking, m: usize, n: usize, a: [*]T, lda: usize, s: ?T, b: [*]T, ldb: usize, w: Work(T)) void {
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
                triPack(T, vt, .trsm, .o, tc, o.unit, min_l, min_l, opAt(T, o.trans, a, lda, ls, ls), lda, 0, 0, 0, w.sb);
                triRun(T, vt, .trsm, o, forward, min_i0, min_l, min_l, du, w.sa, w.sb, b + ls * ldb, ldb, 0);
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
                    triRun(T, vt, .trsm, o, forward, min_i, min_l, min_l, du, w.sa, w.sb, b + is + ls * ldb, ldb, 0);
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
                triPack(T, vt, .trsm, .o, tc, o.unit, min_l, min_l, opAt(T, o.trans, a, lda, ls, ls), lda, 0, 0, 0, tb);
                triRun(T, vt, .trsm, o, forward, min_i0, min_l, min_l, du, w.sa, tb, b + ls * ldb, ldb, 0);
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
                    triRun(T, vt, .trsm, o, forward, min_i, min_l, min_l, du, w.sa, tb, b + is + ls * ldb, ldb, 0);
                    scaleIf(T, vt, s0, min_i, ls - j0, b + is + j0 * ldb, ldb);
                    kernel(T, vt, gc, min_i, ls - j0, min_l, m1, w.sa, w.sb, b + is + j0 * ldb, ldb);
                }
                if (ls < j0 + blk.q) break;
                ls -= blk.q;
            }
        }
    }
}
