//! What the vectors give the moves and measures: FreeType 2.14.3
//! `ttinterp.c` `Compute_Funcs`, `Project`, `Dual_Project` and `Project_x`
//! and `_y`, as values instead of function pointers.

const fixed = @import("../fixed.zig");
const GraphicsState = @import("GraphicsState.zig");

const Projection = @This();

/// Which coordinate a projection or move reads, else both.
pub const Axis = enum { x, y, both };

proj: [2]i32,
dual: [2]i32,
proj_axis: Axis,
dual_axis: Axis,
/// `Direct_Move_X`, `_Y`, or the general move.
move_axis: Axis,
/// The freedom vector over its projection, in 16.16 (`moveVector`).
move: [2]i32,
/// The freedom and projection vectors' dot product, in 2.14.
f_dot_p: i32,

pub fn compute(gs: GraphicsState) Projection {
    const f_dot_p = (gs.proj[0] * gs.free[0] + gs.proj[1] * gs.free[1] + 0x2000) >> 14;
    // commonly collinear, else prohibitively orthogonal, else scaled
    const move: [2]i32 = if (f_dot_p >= 0x3FFE)
        .{ gs.free[0] * 4, gs.free[1] * 4 }
    else if (-0x400 < f_dot_p and f_dot_p < 0x400)
        .{ 0, 0 }
    else
        .{ @divTrunc(gs.free[0] * 0x10000, f_dot_p), @divTrunc(gs.free[1] * 0x10000, f_dot_p) };
    return .{
        .proj = gs.proj,
        .dual = gs.dual,
        .proj_axis = axisOf(gs.proj),
        .dual_axis = axisOf(gs.dual),
        .move_axis = if (f_dot_p >= 0x3FFE and gs.free[0] == 0x4000) .x else if (f_dot_p >= 0x3FFE and gs.free[1] == 0x4000) .y else .both,
        .move = move,
        .f_dot_p = f_dot_p,
    };
}

fn axisOf(v: [2]i32) Axis {
    return if (v[0] == 0x4000) .x else if (v[1] == 0x4000) .y else .both;
}

/// The projection of `(dx, dy)` on the projection vector (`func_project`).
pub fn project(p: Projection, dx: i32, dy: i32) i32 {
    return along(p.proj_axis, p.proj, dx, dy);
}

/// The same on the dual projection vector (`func_dualproj`).
pub fn dualProject(p: Projection, dx: i32, dy: i32) i32 {
    return along(p.dual_axis, p.dual, dx, dy);
}

fn along(axis: Axis, v: [2]i32, dx: i32, dy: i32) i32 {
    return switch (axis) {
        .x => dx,
        .y => dy,
        .both => fixed.dot14(dx, dy, v[0], v[1]),
    };
}
