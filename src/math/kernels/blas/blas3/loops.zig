//! `driver/level3/level3.c`'s js / ls / is loops, once, for the level-3
//! drivers (`driver.zig`, `driver3m.zig`, `half.zig`); each passes a
//! context with its own packing, chunks and kernel calls.

const Blocking = @import("driver.zig").Blocking;

/// The js and ls loops over columns `[n0, n1)` and depth `k`:
/// `ctx.depth(rest)` sizes each depth slice, `ctx.slice(js, min_j, ls,
/// min_l)` does it.
pub fn columns(blk: Blocking, n0: usize, n1: usize, k: usize, ctx: anytype) void {
    var js = n0;
    while (js < n1) : (js += blk.r) {
        const min_j = @min(blk.r, n1 - js);
        var ls: usize = 0;
        var min_l: usize = 0;
        while (ls < k) : (ls += min_l) {
            min_l = ctx.depth(k - ls);
            ctx.slice(js, min_j, ls, min_l);
        }
    }
}

/// The is loop of one depth slice `[ls, ls + min_l)` of columns `[js, js +
/// min_j)`, over rows `[r_lo, r_hi)`: the first row block packs `op(B)` in
/// `ctx.chunk(rest)` columns as it goes, into `sb` as one chunk when that
/// block takes every row; the others reuse it. `s` (`beta`, or null) goes to
/// each `ctx.block`. `ctx.packRows`, `ctx.packCols` and `ctx.block` take the
/// reference's arguments.
pub fn slice(ctx: anytype, blk: Blocking, r_lo: usize, r_hi: usize, js: usize, min_j: usize, ls: usize, min_l: usize, s: anytype, sa: anytype, sb: anytype) void {
    var min_i = Blocking.split(r_hi - r_lo, blk.p, blk.um);
    const once = min_i == r_hi - r_lo;
    ctx.packRows(min_l, min_i, ls, r_lo, sa);
    var jjs = js;
    while (jjs < js + min_j) {
        const min_jj = ctx.chunk(js + min_j - jjs);
        const pb = if (once) sb else sb + min_l * (jjs - js);
        ctx.packCols(min_l, min_jj, ls, jjs, pb);
        ctx.block(r_lo, min_i, jjs, min_jj, min_l, s, sa, pb);
        jjs += min_jj;
    }
    var is = r_lo + min_i;
    while (is < r_hi) : (is += min_i) {
        min_i = Blocking.split(r_hi - is, blk.p, blk.um);
        ctx.packRows(min_l, min_i, ls, is, sa);
        jjs = js;
        while (jjs < js + min_j) {
            const min_jj = ctx.chunk(js + min_j - jjs);
            ctx.block(is, min_i, jjs, min_jj, min_l, s, sa, sb + min_l * (jjs - js));
            jjs += min_jj;
        }
    }
}
