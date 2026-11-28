//! Shifting points by a reference point's move or a distance: FreeType
//! 2.14.3 `ttinterp.c` `Ins_SHP`, `Ins_SHC`, `Ins_SHZ`, `Ins_SHPIX`.

const fixed = @import("../../fixed.zig");
const zone = @import("../zone.zig");
const Engine = @import("../Engine.zig");
const mv = @import("move.zig");

const Error = Engine.Error;

/// The loop count of points below the arguments, taken off the stack; null
/// when there are too few, after resetting the loop.
fn takeLoop(e: *Engine) Error!?usize {
    const loop: usize = @intCast(e.gs.loop);
    if (e.new_top < loop) {
        e.gs.loop = 1;
        try e.pedanticError(error.TooFewArguments);
        return null;
    }
    e.new_top -= loop;
    return loop;
}

/// SHP: `loop` points of `zp2` by the reference's move. A reference out of
/// range leaves the loop count, as FreeType returns before resetting it.
pub fn shp(e: *Engine) Error!void {
    const loop = try takeLoop(e) orelse return;
    const d = try mv.displacement(e) orelse return;
    var at = e.args;
    for (0..loop) |_| {
        at -= 1;
        const p = mv.point16(e.stack[at]);
        if (p >= e.zpZone(2).len()) {
            try e.pedanticError(error.InvalidReference);
            continue;
        }
        mv.moveZp2(e, p, d.dx, d.dy, true);
    }
    e.gs.loop = 1;
}

/// SHC: contour `args[0]` of `zp2`, all of the twilight zone's points.
pub fn shc(e: *Engine) Error!void {
    const contour = mv.point16(e.arg(0).*);
    const z = e.zpZone(2);
    const twilight = e.gs.gep2 == 0;
    const bounds: usize = if (twilight) 1 else z.contours.len;
    if (contour >= bounds) return e.pedanticError(error.InvalidReference);
    const d = try mv.displacement(e) orelse return;
    // a glyph-zone number on the twilight zone has no contours to read
    if (!twilight and contour >= z.contours.len) return;
    const start: usize = if (contour == 0) 0 else @as(usize, z.contours[contour - 1]) + 1;
    const limit: usize = if (twilight) z.len() else @as(usize, z.contours[contour]) + 1;
    for (start..@min(limit, z.len())) |i| {
        if (d.zone != e.zp[2] or d.ref != i) mv.moveZp2(e, i, d.dx, d.dy, true);
    }
}

/// SHZ: every point of `zp2` but its phantoms, untouched; `args[0]` is only
/// checked.
pub fn shz(e: *Engine) Error!void {
    const n: u64 = @bitCast(@as(i64, e.arg(0).*));
    if (n >= 2) return e.pedanticError(error.InvalidReference);
    const d = try mv.displacement(e) orelse return;
    const z = e.zpZone(2);
    const limit: usize = if (e.gs.gep2 == 0)
        z.len()
    else if (e.gs.gep2 == 1 and z.contours.len > 0)
        @as(usize, z.contours[z.contours.len - 1]) + 1
    else
        0;
    for (0..@min(limit, z.len())) |i| {
        if (d.zone != e.zp[2] or d.ref != i) mv.moveZp2(e, i, d.dx, d.dy, false);
    }
}

/// SHPIX: `loop` points of `zp2` by `args[0]` along the freedom vector;
/// under backward compatibility only in y, and only for twilight, touched
/// or composite points.
pub fn shpix(e: *Engine) Error!void {
    const twilight = e.gs.gep0 == 0 or e.gs.gep1 == 0 or e.gs.gep2 == 0;
    const loop = try takeLoop(e) orelse return;
    const dx = fixed.mul14(e.arg(0).*, e.gs.free[0]);
    const dy = fixed.mul14(e.arg(0).*, e.gs.free[1]);
    const bc = e.backward_compatibility;
    var at = e.args;
    for (0..loop) |_| {
        at -= 1;
        const p = mv.point16(e.stack[at]);
        const z = e.zpZone(2);
        if (p >= z.len()) {
            try e.pedanticError(error.InvalidReference);
            continue;
        }
        if (bc == 0) {
            mv.moveZp2(e, p, dx, dy, true);
        } else if (twilight or (bc != 0x7 and ((e.composite and e.gs.free[1] != 0) or z.tags[p] & zone.touched_y != 0))) {
            mv.moveZp2(e, p, 0, dy, true);
        }
    }
    e.gs.loop = 1;
}
