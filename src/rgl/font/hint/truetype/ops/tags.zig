//! Instructions that change point tags without moving points: FreeType
//! 2.14.3 `ttinterp.c` `Ins_FLIPPT`, `Ins_FLIPRGON`, `Ins_FLIPRGOFF`,
//! `Ins_UTP`.

const zone = @import("../zone.zig");
const Engine = @import("../Engine.zig");

const Error = Engine.Error;

/// v40's backward compatibility once IUP has run in both directions
/// freezes the outline's tags.
fn frozen(e: *const Engine) bool {
    return e.backward_compatibility == 0x7;
}

fn point(v: i32) u16 {
    return @truncate(@as(u32, @bitCast(v)));
}

/// FLIPPT: flips the on-curve tag of `loop` glyph points below it.
pub fn flipPoints(e: *Engine) Error!void {
    defer e.gs.loop = 1;
    const loop: usize = @intCast(e.gs.loop);
    if (e.new_top < loop) return e.pedanticError(error.TooFewArguments);
    e.new_top -= loop;
    if (frozen(e)) return;
    const pts = &e.zones[1];
    var at = e.args;
    for (0..loop) |_| {
        at -= 1;
        const p = point(e.stack[at]);
        if (p >= pts.len()) {
            try e.pedanticError(error.InvalidReference);
            continue;
        }
        pts.tags[p] ^= zone.on_curve;
    }
}

/// FLIPRGON and FLIPRGOFF: glyph points `args[0]` to `args[1]` on or off
/// the curve.
pub fn flipRange(e: *Engine, comptime on: bool) Error!void {
    if (frozen(e)) return;
    const k = point(e.arg(1).*);
    const l = point(e.arg(0).*);
    const pts = &e.zones[1];
    if (k >= pts.len() or l >= pts.len()) return e.pedanticError(error.InvalidReference);
    for (l..@as(usize, k) + 1) |i| {
        if (on) pts.tags[i] |= zone.on_curve else pts.tags[i] &= ~zone.on_curve;
    }
}

/// UTP: clears `zp0`'s point's touched tags along the freedom vector.
pub fn untouch(e: *Engine) Error!void {
    const p = point(e.arg(0).*);
    const z = e.zpZone(0);
    if (p >= z.len()) return e.pedanticError(error.InvalidReference);
    var mask: u8 = 0xFF;
    if (e.gs.free[0] != 0) mask &= ~zone.touched_x;
    if (e.gs.free[1] != 0) mask &= ~zone.touched_y;
    z.tags[p] &= mask;
}
