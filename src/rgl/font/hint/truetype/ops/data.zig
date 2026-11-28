//! The storage area, the CVT, measures and queries: FreeType 2.14.3
//! `ttinterp.c` `Ins_RS`, `Ins_WS`, `Ins_WCVTP`, `Ins_WCVTF`, `Ins_RCVT`,
//! `Ins_DELTAC`, `Ins_GC`, `Ins_MD`, `Ins_MPPEM`, `Ins_MPS`, `Ins_GETINFO`,
//! `Ins_GETVARIATION`, `Ins_GETDATA`.

const fixed = @import("../../fixed.zig");
const Engine = @import("../Engine.zig");
const calls = @import("calls.zig");

const Error = Engine.Error;

/// `v` as FreeType's unsigned index, out of range for negatives.
fn index(v: i32) u64 {
    return @bitCast(@as(i64, v));
}

pub fn readStorage(e: *Engine) Error!void {
    const i = index(e.arg(0).*);
    if (i >= e.storage.len()) {
        try e.pedanticError(error.InvalidReference);
        e.arg(0).* = 0;
    } else e.arg(0).* = e.storage.get(@intCast(i));
}

pub fn writeStorage(e: *Engine) Error!void {
    const i = index(e.arg(0).*);
    if (i >= e.storage.len()) return e.pedanticError(error.InvalidReference);
    e.storage.set(e.ini == .glyph, @intCast(i), e.arg(1).*);
}

/// WCVTP: a 26.6 value.
pub fn writeCvtPixels(e: *Engine) Error!void {
    const i = index(e.arg(0).*);
    if (i >= e.cvt.len()) return e.pedanticError(error.InvalidReference);
    e.cvt.set(e.ini == .glyph, @intCast(i), e.arg(1).*);
}

/// WCVTF: font units, scaled; written in place, as FreeType does.
pub fn writeCvtFunits(e: *Engine) Error!void {
    const i = index(e.arg(0).*);
    if (i >= e.cvt.len()) return e.pedanticError(error.InvalidReference);
    e.cvt.setInPlace(@intCast(i), fixed.mul(e.arg(1).*, e.metrics.scale));
}

pub fn readCvt(e: *Engine) Error!void {
    const i = index(e.arg(0).*);
    if (i >= e.cvt.len()) {
        try e.pedanticError(error.InvalidReference);
        e.arg(0).* = 0;
    } else e.arg(0).* = e.cvt.get(@intCast(i));
}

/// DELTAC1 to DELTAC3: `args[0]` pairs below the count, each a CVT index
/// over an argument whose high nibble picks the ppem and low one the step.
pub fn deltaCvt(e: *Engine) Error!void {
    var nump: i64 = e.arg(0).*;
    if (nump < 0 or nump > e.new_top / 2) {
        try e.pedanticError(error.TooFewArguments);
        nump = @intCast(e.new_top / 2);
    }
    const pairs: usize = @intCast(nump);
    e.new_top -= 2 * pairs;
    var p: i32 = e.metrics.ppem - @as(i32, e.gs.delta_base);
    p -= switch (e.op) {
        0x74 => 16,
        0x75 => 32,
        else => 0,
    };
    if (p & ~@as(i32, 0xF) != 0) return;
    p <<= 4;
    const f = @as(i32, 1) << @intCast(6 - e.gs.delta_shift);
    var at = e.args;
    for (0..pairs) |_| {
        at -= 1;
        const a = index(e.stack[at]);
        at -= 1;
        var b = e.stack[at];
        if (a >= e.cvt.len()) {
            try e.pedanticError(error.InvalidReference);
            continue;
        }
        if (b & 0xF0 != p) continue;
        b = (b & 0xF) - 8;
        if (b >= 0) b += 1;
        const i: usize = @intCast(a);
        e.cvt.set(e.ini == .glyph, i, e.cvt.get(i) +% b * f);
    }
}

/// GC: `zp2`'s point `args[0]` projected, current (even) or original (odd).
pub fn getCoordinate(e: *Engine) Error!void {
    const l = index(e.arg(0).*);
    const z = e.zpZone(2);
    if (l >= z.len()) {
        try e.pedanticError(error.InvalidReference);
        e.arg(0).* = 0;
        return;
    }
    const i: usize = @intCast(l);
    e.arg(0).* = if (e.op & 1 != 0) e.projection.dualProject(z.org[i][0], z.org[i][1]) else e.projection.project(z.cur[i][0], z.cur[i][1]);
}

/// MD: the distance from `zp1`'s point `args[1]` to `zp0`'s `args[0]`:
/// current (odd), else original, from font units outside the twilight zone.
pub fn measureDistance(e: *Engine) Error!void {
    const k: u16 = @truncate(@as(u32, @bitCast(e.arg(1).*)));
    const l: u16 = @truncate(@as(u32, @bitCast(e.arg(0).*)));
    const z0 = e.zpZone(0);
    const z1 = e.zpZone(1);
    if (l >= z0.len() or k >= z1.len()) {
        try e.pedanticError(error.InvalidReference);
        e.arg(0).* = 0;
        return;
    }
    const p = e.projection;
    e.arg(0).* = if (e.op & 1 != 0)
        p.project(z0.cur[l][0] -% z1.cur[k][0], z0.cur[l][1] -% z1.cur[k][1])
    else if (e.gs.gep0 == 0 or e.gs.gep1 == 0)
        p.dualProject(z0.org[l][0] -% z1.org[k][0], z0.org[l][1] -% z1.org[k][1])
    else blk: {
        const a = z0.unscaled(l);
        const b = z1.unscaled(k);
        const d = p.dualProject(a[0] -% b[0], a[1] -% b[1]);
        break :blk fixed.mul(d, e.unitsScale());
    };
}

pub fn ppem(e: *Engine) void {
    e.arg(0).* = e.metrics.ppem;
}

/// MPS: the ppem under v35, as GDI's 12 would mislead; the point size
/// under v40.
pub fn pointSize(e: *Engine) void {
    e.arg(0).* = if (e.version == .v35) e.metrics.ppem else e.metrics.point_size;
}

/// GETINFO: the interpreter version and the facts `args[0]` selects.
pub fn getInfo(e: *Engine) void {
    const sel = e.arg(0).*;
    var k: i32 = 0;
    if (sel & 1 != 0) k = if (e.version == .v35) 35 else 40;
    if (sel & 2 != 0 and e.metrics.rotated) k |= 1 << 8;
    if (sel & 4 != 0 and e.metrics.stretched) k |= 1 << 9;
    if (sel & 8 != 0 and e.axis_count != 0) k |= 1 << 10;
    // v40 runs with `grayscale` off (`ttgload.c` `tt_loader_init`)
    if (sel & 32 != 0 and e.version == .v35 and e.target != .mono) k |= 1 << 12;
    if (e.version == .v40 and e.target != .mono) {
        if (sel & 64 != 0) k |= 1 << 13;
        if (sel & 256 != 0 and e.target == .lcd_v) k |= 1 << 15;
        if (sel & 1024 != 0) k |= 1 << 17;
        if (sel & 2048 != 0) k |= 1 << 18;
        if (sel & 4096 != 0 and e.target != .lcd and e.target != .lcd_v) k |= 1 << 19;
    }
    e.arg(0).* = k;
}

/// GETVARIATION: each axis's normalized coordinate in 2.14; for a font
/// that is not variable, its IDEF.
pub fn getVariation(e: *Engine) Error!void {
    if (e.axis_count == 0) return calls.unknown(e);
    const n: usize = e.axis_count;
    if (n > e.stack.len - e.top) return error.StackOverflow;
    for (0..n) |i| e.arg(i).* = if (i < e.coords.len) e.coords[i] >> 2 else 0;
    e.new_top += n;
}

/// GETDATA: 17 for a variable font, else its IDEF.
pub fn getData(e: *Engine) Error!void {
    // FreeType keeps the table's push here too, a stale slot
    if (e.axis_count == 0) return calls.unknown(e);
    e.arg(0).* = 17;
}
