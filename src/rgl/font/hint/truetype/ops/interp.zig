//! Interpolating and adjusting points: FreeType 2.14.3 `ttinterp.c`
//! `Ins_IP`, `Ins_IUP` with `iup_worker_shift_` and
//! `iup_worker_interpolate_`, `Ins_DELTAP`.

const fixed = @import("../../fixed.zig");
const zone = @import("../zone.zig");
const Engine = @import("../Engine.zig");
const mv = @import("move.zig");

const Error = Engine.Error;
const Zone = zone.Zone;

/// IP: `loop` points of `zp2` kept at their original fraction of the way
/// from `rp1` (`zp0`) to `rp2` (`zp1`).
pub fn ip(e: *Engine) Error!void {
    defer e.gs.loop = 1;
    const loop: usize = @intCast(e.gs.loop);
    if (e.new_top < loop) return e.pedanticError(error.TooFewArguments);
    e.new_top -= loop;
    const z0 = e.zpZone(0);
    const z1 = e.zpZone(1);
    const z2 = e.zpZone(2);
    const rp1 = e.gs.rp1;
    const rp2 = e.gs.rp2;
    if (rp1 >= z0.len()) return e.pedanticError(error.InvalidReference);
    // the twilight zone has no font units: its originals stand in
    const twilight = e.gs.gep0 == 0 or e.gs.gep1 == 0 or e.gs.gep2 == 0;
    const orus_base = if (twilight) z0.org[rp1] else z0.unscaled(rp1);
    const cur_base = z0.cur[rp1];
    const p = e.projection;
    var old_range: i32 = 0;
    var cur_range: i32 = 0;
    // fonts with bad rp2 (Monotype's [aeu]grave) get no range
    if (rp2 < z1.len()) {
        const o = if (twilight) z1.org[rp2] else z1.unscaled(rp2);
        old_range = p.dualProject(o[0] -% orus_base[0], o[1] -% orus_base[1]);
        cur_range = p.project(z1.cur[rp2][0] -% cur_base[0], z1.cur[rp2][1] -% cur_base[1]);
    }
    var at = e.args;
    for (0..loop) |_| {
        at -= 1;
        const point: u32 = @bitCast(e.stack[at]);
        if (point >= z2.len()) {
            try e.pedanticError(error.InvalidReference);
            continue;
        }
        const o = if (twilight) z2.org[point] else z2.unscaled(point);
        const org_dist = p.dualProject(o[0] -% orus_base[0], o[1] -% orus_base[1]);
        const cur_dist = p.project(z2.cur[point][0] -% cur_base[0], z2.cur[point][1] -% cur_base[1]);
        const new_dist = if (org_dist == 0) 0 else if (old_range != 0) fixed.mulDiv(org_dist, cur_range, old_range) else org_dist;
        mv.move(e, z2, point, new_dist -% cur_dist);
    }
}

/// IUP: the glyph's untouched points along x (odd) or y, shifted or
/// interpolated between their contour's touched neighbors.
pub fn iup(e: *Engine) void {
    if (e.backward_compatibility == 0x7) return;
    if (e.backward_compatibility != 0) e.backward_compatibility |= @as(u8, 1) << @intCast(e.op & 1);
    const z = &e.zones[1];
    if (z.contours.len == 0) return;
    const c: u1 = if (e.op & 1 != 0) 0 else 1;
    const mask = if (c == 0) zone.touched_x else zone.touched_y;
    const n = z.len();
    if (n == 0) return;
    var point: usize = 0;
    for (z.contours) |end| {
        var end_point: usize = end;
        const first_point = point;
        if (end_point >= n) end_point = n -% 1;
        while (point <= end_point and z.tags[point] & mask == 0) point += 1;
        if (point > end_point) continue;
        const first_touched = point;
        var cur_touched = point;
        point += 1;
        while (point <= end_point) : (point += 1) {
            if (z.tags[point] & mask != 0) {
                interpolate(z, c, cur_touched + 1, point - 1, cur_touched, point);
                cur_touched = point;
            }
        }
        if (cur_touched == first_touched) {
            shift(z, c, first_point, end_point, cur_touched);
        } else {
            interpolate(z, c, cur_touched + 1, end_point, cur_touched, first_touched);
            if (first_touched > 0) interpolate(z, c, first_point, first_touched - 1, cur_touched, first_touched);
        }
    }
}

/// Points `p1` to `p2` but `p` by `p`'s move (`iup_worker_shift_`).
fn shift(z: *Zone, c: u1, p1: usize, p2: usize, p: usize) void {
    const d = z.cur[p][c] -% z.org[p][c];
    if (d == 0) return;
    for (p1..p) |i| z.cur[i][c] +%= d;
    for (p + 1..p2 + 1) |i| z.cur[i][c] +%= d;
}

/// Points `p1` to `p2` placed between `ref1` and `ref2` as their font units
/// lie (`iup_worker_interpolate_`).
fn interpolate(z: *Zone, c: u1, p1: usize, p2: usize, r1: usize, r2: usize) void {
    if (p1 > p2) return;
    const n = z.len();
    if (r1 >= n or r2 >= n) return;
    var ref1 = r1;
    var ref2 = r2;
    var orus1 = z.unscaled(ref1)[c];
    var orus2 = z.unscaled(ref2)[c];
    if (orus1 > orus2) {
        const t = orus1;
        orus1 = orus2;
        orus2 = t;
        ref1 = r2;
        ref2 = r1;
    }
    const org1 = z.org[ref1][c];
    const org2 = z.org[ref2][c];
    const cur1 = z.cur[ref1][c];
    const cur2 = z.cur[ref2][c];
    const delta1 = cur1 -% org1;
    const delta2 = cur2 -% org2;
    const snap = cur1 == cur2 or orus1 == orus2;
    var scale: ?i32 = null;
    for (p1..p2 + 1) |i| {
        const x = z.org[i][c];
        z.cur[i][c] = if (x <= org1)
            x +% delta1
        else if (x >= org2)
            x +% delta2
        else if (snap)
            cur1
        else blk: {
            if (scale == null) scale = fixed.div(cur2 -% cur1, orus2 -% orus1);
            break :blk cur1 +% fixed.mul(z.unscaled(i)[c] -% orus1, scale.?);
        };
    }
}

/// DELTAP1 to DELTAP3: as DELTAC, moving `zp0`'s points; under backward
/// compatibility only touched-in-y or composite points, before IUP ends it.
pub fn deltaPoint(e: *Engine) Error!void {
    var nump: i64 = e.arg(0).*;
    if (nump < 0 or nump > e.new_top / 2) {
        try e.pedanticError(error.TooFewArguments);
        nump = @intCast(e.new_top / 2);
    }
    const pairs: usize = @intCast(nump);
    e.new_top -= 2 * pairs;
    var p: i32 = e.metrics.ppem - @as(i32, e.gs.delta_base);
    p -= switch (e.op) {
        0x71 => 16,
        0x72 => 32,
        else => 0,
    };
    if (p & ~@as(i32, 0xF) != 0) return;
    p <<= 4;
    const f = @as(i32, 1) << @intCast(6 - e.gs.delta_shift);
    const z0 = e.zpZone(0);
    const bc = e.backward_compatibility;
    var at = e.args;
    for (0..pairs) |_| {
        at -= 1;
        const a = mv.point16(e.stack[at]);
        at -= 1;
        var b = e.stack[at];
        if (a >= z0.len()) {
            try e.pedanticError(error.InvalidReference);
            continue;
        }
        if (b & 0xF0 != p) continue;
        b = (b & 0xF) - 8;
        if (b >= 0) b += 1;
        b *%= f;
        if (bc == 0 or (bc != 0x7 and ((e.composite and e.gs.free[1] != 0) or z0.tags[a] & zone.touched_y != 0))) mv.move(e, z0, a, b);
    }
}
