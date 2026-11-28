//! Managing the graphics state: FreeType 2.14.3 `ttinterp.c` `Ins_SxyTCA`
//! to `Ins_SDPVTL`, `Ins_SRP0` to `Ins_SZPS`, `Ins_SLOOP`, the rounding
//! modes and `SetSuperRound`, `Ins_SMD` to `Ins_SDS`, `Ins_INSTCTRL`,
//! `Ins_SCANCTRL`, `Ins_SCANTYPE`.

const fixed = @import("../../fixed.zig");
const zone = @import("../zone.zig");
const round = @import("../round.zig");
const Engine = @import("../Engine.zig");

const Error = Engine.Error;

fn update(e: *Engine) void {
    e.projection = .compute(e.gs);
}

/// SVTCA, SPVTCA, SFVTCA: an axis, x for odd opcodes.
pub fn toAxis(e: *Engine) void {
    const a: i32 = @as(i32, e.op & 1) << 14;
    const b = a ^ 0x4000;
    if (e.op < 4) {
        e.gs.proj = .{ a, b };
        e.gs.dual = .{ a, b };
    }
    if (e.op & 2 == 0) e.gs.free = .{ a, b };
    update(e);
}

/// FreeType's `Normalize`: a zero vector leaves `v` as it was.
fn normalize(x: i32, y: i32, v: *[2]i32) void {
    if (x == 0 and y == 0) return;
    v.* = fixed.normalize14(x, y);
}

/// The vector along (odd opcode: across) the line from `p2` to `p1`; along
/// x when they coincide.
fn lineVector(op: u8, p1: zone.Point, p2: zone.Point, v: *[2]i32) u8 {
    var a = p1[0] -% p2[0];
    var b = p1[1] -% p2[1];
    var o = op;
    if (a == 0 and b == 0) {
        a = 0x4000;
        o = 0;
    }
    if (o & 1 != 0) {
        const c = b;
        b = a;
        a = 0 -% c;
    }
    normalize(a, b, v);
    return o;
}

/// SPVTL and SFVTL (`Ins_SxVTL`): from `zp2`'s point `args[1]` to `zp1`'s
/// `args[0]`.
pub fn toLine(e: *Engine) Error!void {
    const a1: u16 = @truncate(@as(u32, @bitCast(e.arg(1).*)));
    const a2: u16 = @truncate(@as(u32, @bitCast(e.arg(0).*)));
    if (a1 >= e.zpZone(2).len() or a2 >= e.zpZone(1).len()) return e.pedanticError(error.InvalidReference);
    const proj = e.op < 8;
    _ = lineVector(e.op, e.zpZone(1).cur[a2], e.zpZone(2).cur[a1], if (proj) &e.gs.proj else &e.gs.free);
    if (proj) e.gs.dual = e.gs.proj;
    update(e);
}

/// SDPVTL: the dual vector from the original points, the projection
/// vector from the current ones.
pub fn dualToLine(e: *Engine) Error!void {
    const p1: u16 = @truncate(@as(u32, @bitCast(e.arg(1).*)));
    const p2: u16 = @truncate(@as(u32, @bitCast(e.arg(0).*)));
    if (p2 >= e.zpZone(1).len() or p1 >= e.zpZone(2).len()) return e.pedanticError(error.InvalidReference);
    const op = lineVector(e.op, e.zpZone(1).org[p2], e.zpZone(2).org[p1], &e.gs.dual);
    _ = lineVector(op, e.zpZone(1).cur[p2], e.zpZone(2).cur[p1], &e.gs.proj);
    update(e);
}

/// SPVFS and SFVFS: from the low 16 bits of `args[0]`, `args[1]`.
pub fn fromStack(e: *Engine) void {
    const x: i16 = @truncate(e.arg(0).*);
    const y: i16 = @truncate(e.arg(1).*);
    if (e.op == 0x0A) {
        normalize(x, y, &e.gs.proj);
        e.gs.dual = e.gs.proj;
    } else normalize(x, y, &e.gs.free);
    update(e);
}

/// GPV and GFV.
pub fn get(e: *Engine) void {
    const v = if (e.op == 0x0C) e.gs.proj else e.gs.free;
    e.arg(0).* = v[0];
    e.arg(1).* = v[1];
}

pub fn freedomToProjection(e: *Engine) void {
    e.gs.free = e.gs.proj;
    update(e);
}

/// SRP0, SRP1, SRP2.
pub fn referencePoint(e: *Engine, comptime i: u2) void {
    const p: u16 = @truncate(@as(u32, @bitCast(e.arg(0).*)));
    switch (i) {
        0 => e.gs.rp0 = p,
        1 => e.gs.rp1 = p,
        else => e.gs.rp2 = p,
    }
}

/// SZP0, SZP1, SZP2, and SZPS (`i` 3): 0 the twilight zone, 1 the glyph's.
pub fn zonePointer(e: *Engine, comptime i: u2) Error!void {
    const n = e.arg(0).*;
    const z: zone.Pointer = switch (n) {
        0 => .twilight,
        1 => .glyph,
        else => return e.pedanticError(error.InvalidReference),
    };
    const g: u16 = @intCast(n);
    if (i == 3) {
        e.zp = @splat(z);
        e.gs.gep0 = g;
        e.gs.gep1 = g;
        e.gs.gep2 = g;
        return;
    }
    e.zp[i] = z;
    switch (i) {
        0 => e.gs.gep0 = g,
        1 => e.gs.gep1 = g,
        else => e.gs.gep2 = g,
    }
}

pub fn setLoop(e: *Engine) Error!void {
    const n = e.arg(0).*;
    if (n < 0) return error.BadArgument;
    // FreeType's limit of 16 bits
    e.gs.loop = @min(n, 0xFFFF);
}

/// RTG, RTHG, RTDG, RDTG, RUTG, ROFF.
pub fn setRound(e: *Engine, comptime mode: round.Mode) void {
    e.gs.round_state = mode;
    e.round_mode = mode;
}

/// SROUND (`grid` 0x4000) and S45ROUND (0x2D41): `SetSuperRound`.
pub fn superRound(e: *Engine, comptime mode: round.Mode) void {
    const grid: i32 = if (mode == .super) 0x4000 else 0x2D41;
    const sel = e.arg(0).*;
    const period: i32 = switch (sel & 0xC0) {
        0 => @divTrunc(grid, 2),
        0x80 => grid * 2,
        else => grid,
    };
    const phase: i32 = switch (sel & 0x30) {
        0 => 0,
        0x10 => @divTrunc(period, 4),
        0x20 => @divTrunc(period, 2),
        else => @divTrunc(period * 3, 4),
    };
    const threshold: i32 = if (sel & 0x0F == 0) period - 1 else @divTrunc(((sel & 0x0F) - 4) * period, 8);
    e.super_round = .{ .period = period >> 8, .phase = phase >> 8, .threshold = threshold >> 8 };
    setRound(e, mode);
}

/// SMD, SCVTCI, SSWCI: a 26.6 setting.
pub fn setDistance(e: *Engine) void {
    const v = e.arg(0).*;
    switch (e.op) {
        0x1A => e.gs.min_distance = v,
        0x1D => e.gs.control_value_cutin = v,
        else => e.gs.single_width_cutin = v,
    }
}

/// SSW: font units, scaled.
pub fn singleWidth(e: *Engine) void {
    e.gs.single_width_value = fixed.mul(e.arg(0).*, e.metrics.scale);
}

pub fn deltaBase(e: *Engine) void {
    e.gs.delta_base = @truncate(@as(u32, @bitCast(e.arg(0).*)));
}

pub fn deltaShift(e: *Engine) Error!void {
    const n: u32 = @bitCast(e.arg(0).*);
    if (n > 6) return error.BadArgument;
    e.gs.delta_shift = @intCast(n);
}

/// INSTCTRL: selector `args[1]` (1 to 3), value `args[0]` (0 or its flag);
/// in `prep` only, but for a glyph program's selector 3 under v40.
pub fn instructControl(e: *Engine) Error!void {
    const k: u32 = @bitCast(e.arg(1).*);
    const l: u32 = @bitCast(e.arg(0).*);
    if (k < 1 or k > 3) return e.pedanticError(error.InvalidReference);
    const kf = @as(u32, 1) << @intCast(k - 1);
    if (l != 0 and l != kf) return e.pedanticError(error.InvalidReference);
    if (e.ini == .control_value) {
        e.gs.instruct_control &= ~@as(u8, @intCast(kf));
        e.gs.instruct_control |= @intCast(l);
    } else if (e.ini == .glyph and k == 3) {
        if (e.version == .v40) e.backward_compatibility = @as(u8, @intCast(l & 4)) ^ 4;
    } else return e.pedanticError(error.InvalidReference);
}

/// SCANCTRL: dropout control by the ppem threshold in the low byte and the
/// condition bits above it.
pub fn scanControl(e: *Engine) void {
    const n = e.arg(0).*;
    const a = n & 0xFF;
    if (a == 0xFF) {
        e.gs.scan_control = true;
        return;
    }
    if (a == 0) {
        e.gs.scan_control = false;
        return;
    }
    const m = e.metrics;
    if (n & 0x100 != 0 and m.ppem <= a) e.gs.scan_control = true;
    if (n & 0x200 != 0 and m.rotated) e.gs.scan_control = true;
    if (n & 0x400 != 0 and m.stretched) e.gs.scan_control = true;
    if (n & 0x800 != 0 and m.ppem > a) e.gs.scan_control = false;
    if (n & 0x1000 != 0 and m.rotated) e.gs.scan_control = false;
    if (n & 0x2000 != 0 and m.stretched) e.gs.scan_control = false;
}

pub fn scanType(e: *Engine) void {
    const n = e.arg(0).*;
    if (n >= 0) e.gs.scan_type = n & 0xFFFF;
}
