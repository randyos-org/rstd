//! `trmm`'s drivers: `driver/level3/trmm_L.c` and `trmm_R.c` on the shared blocking and packing.

const Vtable = @import("../../Vtable.zig");
const util = @import("../../util.zig");
const driver = @import("../driver.zig");
const troot = @import("../triangular.zig");

const Blocking = driver.Blocking;
const Work = driver.Work;
const kernel = driver.kernel;
const rowsCopy = driver.rowsCopy;
const colsCopy = driver.colsCopy;
const Tri = troot.Tri;
const TriCopy = troot.TriCopy;
const opAt = troot.opAt;
const scaleIf = troot.scaleIf;
const triPack = troot.triPack;
const triRun = troot.triRun;

/// `min_i` of the `trmm` drivers: at most `P`, whole row panels past one.
fn trmmRows(blk: Blocking, rest: usize) usize {
    const v = @min(rest, blk.p);
    return if (v > blk.um) v / blk.um * blk.um else v;
}

pub fn trmmL(comptime T: type, vt: *const Vtable, o: Tri, blk: Blocking, m: usize, n: usize, a: [*]T, lda: usize, s: ?T, b: [*]T, ldb: usize, w: Work(T)) void {
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
            triPack(T, vt, .trmm, .i, tc, o.unit, min_l, min_i, a, lda, 0, 0, 0, w.sa);
            var jjs = js;
            while (jjs < js + min_j) {
                const min_jj = blk.chunk(js + min_j - jjs);
                const pb = w.sb + min_l * (jjs - js);
                scaleIf(T, vt, s, min_l, min_jj, b + jjs * ldb, ldb);
                colsCopy(T, vt, false, min_l, min_jj, b + jjs * ldb, ldb, pb);
                triRun(T, vt, .trmm, o, down, min_i, min_jj, min_l, p1, w.sa, pb, b + jjs * ldb, ldb, 0);
                jjs += min_jj;
            }
            var is = min_i;
            while (is < min_l) : (is += min_i) {
                min_i = trmmRows(blk, min_l - is);
                triPack(T, vt, .trmm, .i, tc, o.unit, min_l, min_i, a, lda, 0, is, 0, w.sa);
                triRun(T, vt, .trmm, o, down, min_i, min_j, min_l, p1, w.sa, w.sb, b + is + js * ldb, ldb, @intCast(is));
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
                    triPack(T, vt, .trmm, .i, tc, o.unit, min_l, min_i, a, lda, ls, is, 0, w.sa);
                    triRun(T, vt, .trmm, o, down, min_i, min_j, min_l, p1, w.sa, w.sb, b + is + js * ldb, ldb, @intCast(is - ls));
                }
            }
        } else {
            // `op(A)` lower: the slices go up, each row taking the depths up to it (`LT`)
            var min_l = @min(m, blk.q);
            var min_i = trmmRows(blk, min_l);
            triPack(T, vt, .trmm, .i, tc, o.unit, min_l, min_i, a, lda, m - min_l, m - min_l, 0, w.sa);
            var jjs = js;
            while (jjs < js + min_j) {
                const min_jj = blk.chunk(js + min_j - jjs);
                const pb = w.sb + min_l * (jjs - js);
                scaleIf(T, vt, s, min_l, min_jj, b + (m - min_l) + jjs * ldb, ldb);
                colsCopy(T, vt, false, min_l, min_jj, b + (m - min_l) + jjs * ldb, ldb, pb);
                triRun(T, vt, .trmm, o, down, min_i, min_jj, min_l, p1, w.sa, pb, b + (m - min_l) + jjs * ldb, ldb, 0);
                jjs += min_jj;
            }
            var is = m - min_l + min_i;
            while (is < m) : (is += min_i) {
                min_i = trmmRows(blk, m - is);
                triPack(T, vt, .trmm, .i, tc, o.unit, min_l, min_i, a, lda, m - min_l, is, 0, w.sa);
                triRun(T, vt, .trmm, o, down, min_i, min_j, min_l, p1, w.sa, w.sb, b + is + js * ldb, ldb, @intCast(is - (m - min_l)));
            }
            var ls = m - min_l;
            while (ls > 0) : (ls -= min_l) {
                min_l = @min(ls, blk.q);
                const l0 = ls - min_l;
                min_i = trmmRows(blk, min_l);
                triPack(T, vt, .trmm, .i, tc, o.unit, min_l, min_i, a, lda, l0, l0, 0, w.sa);
                jjs = js;
                while (jjs < js + min_j) {
                    const min_jj = blk.chunk(js + min_j - jjs);
                    const pb = w.sb + min_l * (jjs - js);
                    scaleIf(T, vt, s, min_l, min_jj, b + l0 + jjs * ldb, ldb);
                    colsCopy(T, vt, false, min_l, min_jj, b + l0 + jjs * ldb, ldb, pb);
                    triRun(T, vt, .trmm, o, down, min_i, min_jj, min_l, p1, w.sa, pb, b + l0 + jjs * ldb, ldb, 0);
                    jjs += min_jj;
                }
                is = l0 + min_i;
                while (is < ls) : (is += min_i) {
                    min_i = trmmRows(blk, ls - is);
                    triPack(T, vt, .trmm, .i, tc, o.unit, min_l, min_i, a, lda, l0, is, 0, w.sa);
                    triRun(T, vt, .trmm, o, down, min_i, min_j, min_l, p1, w.sa, w.sb, b + is + js * ldb, ldb, @intCast(is - l0));
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

pub fn trmmR(comptime T: type, vt: *const Vtable, o: Tri, blk: Blocking, m: usize, n: usize, a: [*]T, lda: usize, s: ?T, b: [*]T, ldb: usize, w: Work(T)) void {
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
                    triPack(T, vt, .trmm, .o, tc, o.unit, min_l, min_jj, a, lda, ls, ls + jjs, 0, pb);
                    triRun(T, vt, .trmm, o, !right, min_i0, min_jj, min_l, p1, w.sa, pb, b + (ls + jjs) * ldb, ldb, -@as(isize, @intCast(jjs)));
                    jjs += min_jj;
                }
                var is = min_i0;
                while (is < m) : (is += blk.p) {
                    const min_i = @min(m - is, blk.p);
                    scaleIf(T, vt, s0, min_i, min_l, b + is + ls * ldb, ldb);
                    rowsCopy(T, vt, false, min_l, min_i, b + is + ls * ldb, ldb, w.sa);
                    kernel(T, vt, gc, min_i, ls - js, min_l, p1, w.sa, w.sb, b + is + js * ldb, ldb);
                    triRun(T, vt, .trmm, o, !right, min_i, min_l, min_l, p1, w.sa, w.sb + (ls - js) * min_l, b + is + ls * ldb, ldb, 0);
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
                    triPack(T, vt, .trmm, .o, tc, o.unit, min_l, min_jj, a, lda, ls, ls + jjs, 0, w.sb + min_l * jjs);
                    triRun(T, vt, .trmm, o, !right, min_i0, min_jj, min_l, p1, w.sa, w.sb + min_l * jjs, b + (ls + jjs) * ldb, ldb, -@as(isize, @intCast(jjs)));
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
                    triRun(T, vt, .trmm, o, !right, min_i, min_l, min_l, p1, w.sa, w.sb, b + is + ls * ldb, ldb, 0);
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
