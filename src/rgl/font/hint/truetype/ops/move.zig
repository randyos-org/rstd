//! Moving points: FreeType 2.14.3 `ttinterp.c` `Direct_Move`, `_X`, `_Y`,
//! their `_Orig` forms, `Move_Zp2_Point` and `Compute_Point_Displacement`,
//! with v40's backward compatibility: no x moves while it is on, no y moves
//! once IUP has run in both directions.

const fixed = @import("../../fixed.zig");
const zone = @import("../zone.zig");
const Engine = @import("../Engine.zig");

const Zone = zone.Zone;

/// The move along the freedom vector by a distance measured along the
/// projection vector (`func_move`), touching the point.
pub fn move(e: *const Engine, z: *Zone, point: usize, d: i32) void {
    const bc = e.backward_compatibility;
    const p = &z.cur[point];
    switch (e.projection.move_axis) {
        .x => {
            if (bc == 0) p[0] +%= d;
            z.tags[point] |= zone.touched_x;
        },
        .y => {
            if (bc != 0x7) p[1] +%= d;
            z.tags[point] |= zone.touched_y;
        },
        .both => {
            const v = e.projection.move;
            if (v[0] != 0) {
                if (bc == 0) p[0] +%= fixed.mul(d, v[0]);
                z.tags[point] |= zone.touched_x;
            }
            if (v[1] != 0) {
                if (bc != 0x7) p[1] +%= fixed.mul(d, v[1]);
                z.tags[point] |= zone.touched_y;
            }
        },
    }
}

/// The same on the original point, untouched (`func_move_orig`).
pub fn moveOriginal(e: *const Engine, z: *Zone, point: usize, d: i32) void {
    const p = &z.org[point];
    switch (e.projection.move_axis) {
        .x => p[0] +%= d,
        .y => p[1] +%= d,
        .both => {
            const v = e.projection.move;
            if (v[0] != 0) p[0] +%= fixed.mul(d, v[0]);
            if (v[1] != 0) p[1] +%= fixed.mul(d, v[1]);
        },
    }
}

/// `zp2`'s point by `(dx, dy)` along the freedom vector's axes
/// (`Move_Zp2_Point`).
pub fn moveZp2(e: *Engine, point: usize, dx: i32, dy: i32, touch: bool) void {
    const bc = e.backward_compatibility;
    const z = e.zpZone(2);
    if (e.gs.free[0] != 0) {
        if (bc == 0) z.cur[point][0] +%= dx;
        if (touch) z.tags[point] |= zone.touched_x;
    }
    if (e.gs.free[1] != 0) {
        if (bc != 0x7) z.cur[point][1] +%= dy;
        if (touch) z.tags[point] |= zone.touched_y;
    }
}

/// How far a reference point moved (`Compute_Point_Displacement`): `rp1`
/// in `zp0` for odd opcodes, else `rp2` in `zp1`.
pub const Displacement = struct { zone: zone.Pointer, ref: u16, dx: i32, dy: i32 };

/// The displacement, null for a reference out of range.
pub fn displacement(e: *Engine) Engine.Error!?Displacement {
    const odd = e.op & 1 != 0;
    const which = if (odd) e.zp[0] else e.zp[1];
    const p = if (odd) e.gs.rp1 else e.gs.rp2;
    const z = e.zoneAt(which);
    if (p >= z.len()) {
        try e.pedanticError(error.InvalidReference);
        return null;
    }
    const d = e.projection.project(z.cur[p][0] -% z.org[p][0], z.cur[p][1] -% z.org[p][1]);
    return .{ .zone = which, .ref = p, .dx = fixed.mul(d, e.projection.move[0]), .dy = fixed.mul(d, e.projection.move[1]) };
}

/// `v` as FreeType's `(FT_UShort)` point index.
pub fn point16(v: i32) u16 {
    return @truncate(@as(u32, @bitCast(v)));
}
