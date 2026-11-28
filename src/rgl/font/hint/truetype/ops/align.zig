//! Placing points by a reference, the CVT or each other: FreeType 2.14.3
//! `ttinterp.c` `Ins_MSIRP`, `Ins_MDAP`, `Ins_MIAP`, `Ins_MDRP`,
//! `Ins_MIRP`, `Ins_ALIGNRP`, `Ins_ALIGNPTS`, `Ins_ISECT`, `Ins_SCFS`.

const fixed = @import("../../fixed.zig");
const zone = @import("../zone.zig");
const Engine = @import("../Engine.zig");
const mv = @import("move.zig");

const Error = Engine.Error;
const Zone = zone.Zone;

fn projectDiff(e: *const Engine, a: zone.Point, b: zone.Point) i32 {
    return e.projection.project(a[0] -% b[0], a[1] -% b[1]);
}

fn dualDiff(e: *const Engine, a: zone.Point, b: zone.Point) i32 {
    return e.projection.dualProject(a[0] -% b[0], a[1] -% b[1]);
}

fn abs(v: i32) i32 {
    return if (v < 0) 0 -% v else v;
}

/// MSIRP: `zp1`'s point `args[0]` at distance `args[1]` from `rp0`.
pub fn msirp(e: *Engine) Error!void {
    const point = mv.point16(e.arg(0).*);
    const z0 = e.zpZone(0);
    const z1 = e.zpZone(1);
    const rp0 = e.gs.rp0;
    if (point >= z1.len() or rp0 >= z0.len()) return e.pedanticError(error.InvalidReference);
    // the twilight point starts at rp0's original, as Microsoft's does
    if (e.gs.gep1 == 0) {
        z1.org[point] = z0.org[rp0];
        mv.moveOriginal(e, z1, point, e.arg(1).*);
        z1.cur[point] = z1.org[point];
    }
    const distance = projectDiff(e, z1.cur[point], z0.cur[rp0]);
    mv.move(e, z1, point, e.arg(1).* -% distance);
    e.gs.rp1 = rp0;
    e.gs.rp2 = point;
    if (e.op & 1 != 0) e.gs.rp0 = point;
}

/// MDAP: touches `zp0`'s point `args[0]`, rounded for odd opcodes.
pub fn mdap(e: *Engine) Error!void {
    const point = mv.point16(e.arg(0).*);
    const z0 = e.zpZone(0);
    if (point >= z0.len()) return e.pedanticError(error.InvalidReference);
    var distance: i32 = 0;
    if (e.op & 1 != 0) {
        const cur = e.projection.project(z0.cur[point][0], z0.cur[point][1]);
        distance = e.roundDistance(cur) -% cur;
    }
    mv.move(e, z0, point, distance);
    e.gs.rp0 = point;
    e.gs.rp1 = point;
}

/// MIAP: `zp0`'s point `args[0]` to CVT `args[1]`, with the cut-in and
/// rounding for odd opcodes; a twilight point starts at the CVT distance.
pub fn miap(e: *Engine) Error!void {
    const entry: u64 = @bitCast(@as(i64, e.arg(1).*));
    const point = mv.point16(e.arg(0).*);
    defer {
        e.gs.rp0 = point;
        e.gs.rp1 = point;
    }
    const z0 = e.zpZone(0);
    if (point >= z0.len() or entry >= e.cvt.len()) return e.pedanticError(error.InvalidReference);
    var distance = e.cvt.get(@intCast(entry));
    if (e.gs.gep0 == 0) {
        z0.org[point] = .{ fixed.mul14(distance, e.gs.free[0]), fixed.mul14(distance, e.gs.free[1]) };
        z0.cur[point] = z0.org[point];
    }
    const org_dist = e.projection.project(z0.cur[point][0], z0.cur[point][1]);
    if (e.op & 1 != 0) {
        if (abs(distance -% org_dist) > e.gs.control_value_cutin) distance = org_dist;
        distance = e.roundDistance(distance);
    }
    mv.move(e, z0, point, distance -% org_dist);
}

/// The minimum distance flag (`opcode & 8`) on `distance`, by `org`'s sign.
fn minimumDistance(e: *const Engine, org: i32, distance: i32) i32 {
    if (e.op & 8 == 0) return distance;
    const md = e.gs.min_distance;
    if (org >= 0) return @max(distance, md);
    return if (distance > 0 -% md) 0 -% md else distance;
}

/// The reference points after MDRP and MIRP, also when they fail.
fn afterRelative(e: *Engine, point: u16) void {
    e.gs.rp1 = e.gs.rp0;
    e.gs.rp2 = point;
    if (e.op & 16 != 0) e.gs.rp0 = point;
}

/// MDRP: `zp1`'s point `args[0]` at its original distance from `rp0`,
/// with the single width cut-in, rounding and minimum distance by flags.
pub fn mdrp(e: *Engine) Error!void {
    const point = mv.point16(e.arg(0).*);
    defer afterRelative(e, point);
    const z0 = e.zpZone(0);
    const z1 = e.zpZone(1);
    const rp0 = e.gs.rp0;
    if (point >= z1.len() or rp0 >= z0.len()) return e.pedanticError(error.InvalidReference);
    var org_dist = if (e.gs.gep0 == 0 or e.gs.gep1 == 0)
        dualDiff(e, z1.org[point], z0.org[rp0])
    else
        fixed.mul(dualDiff(e, z1.unscaled(point), z0.unscaled(rp0)), e.unitsScale());
    const sw = e.gs.single_width_value;
    const cutin = e.gs.single_width_cutin;
    if (cutin > 0 and org_dist < sw +% cutin and org_dist > sw -% cutin) org_dist = if (org_dist >= 0) sw else 0 -% sw;
    var distance = if (e.op & 4 != 0) e.roundDistance(org_dist) else org_dist;
    distance = minimumDistance(e, org_dist, distance);
    const cur = projectDiff(e, z1.cur[point], z0.cur[rp0]);
    mv.move(e, z1, point, distance -% cur);
}

/// MIRP: `zp1`'s point `args[0]` at CVT `args[1]`'s distance from `rp0`
/// (CVT -1 reads 0), with the single width test, auto flip, cut-in,
/// rounding and minimum distance.
pub fn mirp(e: *Engine) Error!void {
    const point = mv.point16(e.arg(0).*);
    defer afterRelative(e, point);
    const entry: u64 = @bitCast(@as(i64, e.arg(1).* +% 1));
    const z0 = e.zpZone(0);
    const z1 = e.zpZone(1);
    const rp0 = e.gs.rp0;
    if (point >= z1.len() or entry >= e.cvt.len() + 1 or rp0 >= z0.len()) return e.pedanticError(error.InvalidReference);
    var cvt_dist: i32 = if (entry == 0) 0 else e.cvt.get(@intCast(entry - 1));
    const sw = e.gs.single_width_value;
    if (abs(cvt_dist -% sw) < e.gs.single_width_cutin) cvt_dist = if (cvt_dist >= 0) sw else 0 -% sw;
    if (e.gs.gep1 == 0) {
        const o = z0.org[rp0];
        z1.org[point] = .{ o[0] +% fixed.mul14(cvt_dist, e.gs.free[0]), o[1] +% fixed.mul14(cvt_dist, e.gs.free[1]) };
        z1.cur[point] = z1.org[point];
    }
    const org_dist = dualDiff(e, z1.org[point], z0.org[rp0]);
    const cur_dist = projectDiff(e, z1.cur[point], z0.cur[rp0]);
    if (e.gs.auto_flip and (org_dist ^ cvt_dist) < 0) cvt_dist = 0 -% cvt_dist;
    var distance = cvt_dist;
    if (e.op & 4 != 0) {
        // the cut-in only within one zone
        if (e.gs.gep0 == e.gs.gep1 and abs(cvt_dist -% org_dist) > e.gs.control_value_cutin) cvt_dist = org_dist;
        distance = e.roundDistance(cvt_dist);
    }
    distance = minimumDistance(e, org_dist, distance);
    mv.move(e, z1, point, distance -% cur_dist);
}

/// ALIGNRP: `loop` points of `zp1` onto `rp0` along the projection vector.
pub fn alignrp(e: *Engine) Error!void {
    defer e.gs.loop = 1;
    const loop: usize = @intCast(e.gs.loop);
    if (e.new_top < loop) return e.pedanticError(error.TooFewArguments);
    e.new_top -= loop;
    const z0 = e.zpZone(0);
    const z1 = e.zpZone(1);
    const rp0 = e.gs.rp0;
    if (rp0 >= z0.len()) return e.pedanticError(error.InvalidReference);
    var at = e.args;
    for (0..loop) |_| {
        at -= 1;
        const p = mv.point16(e.stack[at]);
        if (p >= z1.len()) {
            try e.pedanticError(error.InvalidReference);
            continue;
        }
        mv.move(e, z1, p, 0 -% projectDiff(e, z1.cur[p], z0.cur[rp0]));
    }
}

/// ALIGNPTS: `zp1`'s `args[0]` and `zp0`'s `args[1]` to their midpoint.
pub fn alignpts(e: *Engine) Error!void {
    const p1 = mv.point16(e.arg(0).*);
    const p2 = mv.point16(e.arg(1).*);
    const z0 = e.zpZone(0);
    const z1 = e.zpZone(1);
    if (p1 >= z1.len() or p2 >= z0.len()) return e.pedanticError(error.InvalidReference);
    const distance = @divTrunc(projectDiff(e, z0.cur[p2], z1.cur[p1]), 2);
    mv.move(e, z1, p1, distance);
    mv.move(e, z0, p2, 0 -% distance);
}

/// ISECT: `zp2`'s point `args[0]` to where line `a0 a1` of `zp1` meets
/// `b0 b1` of `zp0`, or the middle of their middles when they near
/// parallel (within 3 degrees, `|tan| < 1/19`).
pub fn isect(e: *Engine) Error!void {
    const point = mv.point16(e.arg(0).*);
    const a0 = mv.point16(e.arg(1).*);
    const a1 = mv.point16(e.arg(2).*);
    const b0 = mv.point16(e.arg(3).*);
    const b1 = mv.point16(e.arg(4).*);
    const z0 = e.zpZone(0);
    const z1 = e.zpZone(1);
    const z2 = e.zpZone(2);
    if (b0 >= z0.len() or b1 >= z0.len() or a0 >= z1.len() or a1 >= z1.len() or point >= z2.len()) return e.pedanticError(error.InvalidReference);
    const dbx = z0.cur[b1][0] -% z0.cur[b0][0];
    const dby = z0.cur[b1][1] -% z0.cur[b0][1];
    const dax = z1.cur[a1][0] -% z1.cur[a0][0];
    const day = z1.cur[a1][1] -% z1.cur[a0][1];
    const dx = z0.cur[b0][0] -% z1.cur[a0][0];
    const dy = z0.cur[b0][1] -% z1.cur[a0][1];
    const disc = fixed.mulDiv(dax, 0 -% dby, 0x40) +% fixed.mulDiv(day, dbx, 0x40);
    const dot = fixed.mulDiv(dax, dbx, 0x40) +% fixed.mulDiv(day, dby, 0x40);
    if (19 *% abs(disc) > abs(dot)) {
        const val = fixed.mulDiv(dx, 0 -% dby, 0x40) +% fixed.mulDiv(dy, dbx, 0x40);
        const rx = fixed.mulDiv(val, dax, disc);
        const ry = fixed.mulDiv(val, day, disc);
        z2.cur[point] = .{ z1.cur[a0][0] +% rx, z1.cur[a0][1] +% ry };
    } else {
        inline for (0..2) |c| z2.cur[point][c] = @divTrunc((z1.cur[a0][c] +% z1.cur[a1][c]) +% (z0.cur[b0][c] +% z0.cur[b1][c]), 4);
    }
    z2.tags[point] |= zone.touched_x | zone.touched_y;
}

/// SCFS: `zp2`'s point `args[0]` to coordinate `args[1]` along the
/// projection vector; a twilight point's original follows.
pub fn scfs(e: *Engine) Error!void {
    const l = mv.point16(e.arg(0).*);
    const z2 = e.zpZone(2);
    if (l >= z2.len()) return e.pedanticError(error.InvalidReference);
    const k = e.projection.project(z2.cur[l][0], z2.cur[l][1]);
    mv.move(e, z2, l, e.arg(1).* -% k);
    if (e.gs.gep2 == 0) z2.org[l] = z2.cur[l];
}
