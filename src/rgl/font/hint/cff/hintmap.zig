//! Stem hints, hint masks and the hint map that moves y coordinates to the
//! pixel grid: FreeType 2.14.3 `src/psaux/pshints.c` (`cf2_hint_init`,
//! `cf2_hintmap_map`, `cf2_hintmap_insertHint`, `cf2_hintmap_adjustHints`,
//! `cf2_hintmap_build`), with stem darkening off; after fontations
//! `skrifa/src/outline/cff/hint.rs` (`Hint`, `HintMap`, commit 4f115482,
//! MIT). Values are 16.16.

const std = @import("std");
const fixed = @import("../fixed.zig");
const Blues = @import("blues.zig").Blues;

/// Stem hints a glyph may have, h and v together (`CF2_MAX_HINTS`).
pub const max_hints = 96;
/// Edges a hint map holds (`CF2_MAX_HINT_EDGES`).
const max_edges = 2 * max_hints;
const min_counter = @import("blues.zig").min_counter;

/// A stem from HSTEM or HSTEMHM, and where the hint map placed it.
pub const StemHint = struct {
    min: i32,
    max: i32,
    used: bool = false,
    ds_min: i32 = 0,
    ds_max: i32 = 0,
};

/// One edge of a stem hint (`CF2_HintRec`).
pub const Hint = struct {
    pub const ghost_bottom: u8 = 0x1;
    pub const ghost_top: u8 = 0x2;
    pub const pair_bottom: u8 = 0x4;
    pub const pair_top: u8 = 0x8;
    pub const locked: u8 = 0x10;
    pub const synthetic: u8 = 0x20;

    flags: u8 = 0,
    index: u8 = 0,
    cs: i32 = 0,
    ds: i32 = 0,
    scale: i32 = 0,

    pub fn isBottom(h: Hint) bool {
        return h.flags & (ghost_bottom | pair_bottom) != 0;
    }

    pub fn isTop(h: Hint) bool {
        return h.flags & (ghost_top | pair_top) != 0;
    }

    fn isPair(h: Hint) bool {
        return h.flags & (pair_bottom | pair_top) != 0;
    }

    fn isLocked(h: Hint) bool {
        return h.flags & locked != 0;
    }

    /// The bottom or top edge of `stem` (`cf2_hint_init`); widths -21 and
    /// -20 are ghost hints.
    fn init(stem: StemHint, index: usize, origin: i32, scale: i32, bottom: bool) Hint {
        var h: Hint = .{};
        const width = stem.max -% stem.min;
        if (width == -21 << 16) {
            if (bottom) h = .{ .cs = stem.max, .flags = ghost_bottom };
        } else if (width == -20 << 16) {
            if (!bottom) h = .{ .cs = stem.min, .flags = ghost_top };
        } else if (width < 0) {
            // an inverted pair: its edges swapped
            h = if (bottom) .{ .cs = stem.max, .flags = pair_bottom } else .{ .cs = stem.min, .flags = pair_top };
        } else {
            h = if (bottom) .{ .cs = stem.min, .flags = pair_bottom } else .{ .cs = stem.max, .flags = pair_top };
        }
        h.cs +%= origin;
        h.scale = scale;
        h.index = @intCast(index);
        if (h.flags != 0 and stem.used) {
            h.ds = if (h.isTop()) stem.ds_max else stem.ds_min;
            h.flags |= locked;
        } else {
            h.ds = fixed.mul(h.cs, scale);
        }
        return h;
    }
};

/// Which stem hints are active, most significant bit first.
pub const HintMask = struct {
    bytes: [max_hints / 8]u8 = @splat(0xFF),

    pub fn fromBytes(b: []const u8) HintMask {
        var m: HintMask = .{ .bytes = @splat(0) };
        @memcpy(m.bytes[0..b.len], b);
        return m;
    }

    fn get(m: HintMask, i: usize) bool {
        return m.bytes[i >> 3] & bit(i) != 0;
    }

    fn clear(m: *HintMask, i: usize) void {
        m.bytes[i >> 3] &= ~bit(i);
    }

    fn bit(i: usize) u8 {
        return @as(u8, 0x80) >> @intCast(i & 7);
    }
};

/// The sorted hint edges, each mapping character space to device space
/// from its own coordinate up.
pub const HintMap = struct {
    edges: [max_edges]Hint = undefined,
    len: usize = 0,
    valid: bool = false,
    scale: i32,
    /// Uniform scaling for an unhinted load.
    hinted: bool,

    /// `cf2_hintmap_map`.
    pub fn map(m: *const HintMap, cs: i32) i32 {
        if (m.len == 0 or !m.hinted) return fixed.mul(cs, m.scale);
        var i: usize = 0;
        while (i < m.len - 1 and cs >= m.edges[i + 1].cs) i += 1;
        while (i > 0 and cs < m.edges[i].cs) i -= 1;
        // below the first edge: the uniform scale
        if (i == 0 and cs < m.edges[0].cs) return fixed.mul(cs -% m.edges[0].cs, m.scale) +% m.edges[0].ds;
        return fixed.mul(cs -% m.edges[i].cs, m.edges[i].scale) +% m.edges[i].ds;
    }

    /// `cf2_hintmap_insertHint`; overlapping hints are dropped.
    fn insert(m: *HintMap, bottom: Hint, top: Hint, initial: ?*const HintMap) void {
        var pair = true;
        var first = bottom;
        var second = top;
        if (bottom.flags == 0) {
            first = top;
            pair = false;
        } else if (top.flags == 0) {
            pair = false;
        }
        if (pair and top.cs < bottom.cs) return;
        var at: usize = 0;
        while (at < m.len and m.edges[at].cs < first.cs) at += 1;
        if (at < m.len) {
            const cur = m.edges[at];
            if (cur.cs == first.cs or (pair and cur.cs <= second.cs) or cur.flags & Hint.pair_top != 0) return;
        }
        if (initial) |im| if (!first.isLocked()) {
            if (pair) {
                // the stem's middle from the initial map, its width at the nominal scale
                const mid = im.map(first.cs +% @divTrunc(second.cs -% first.cs, 2));
                const half = fixed.mul(@divTrunc(second.cs -% first.cs, 2), m.scale);
                first.ds = mid -% half;
                second.ds = mid +% half;
            } else {
                first.ds = im.map(first.cs);
            }
        };
        if (at > 0 and first.ds < m.edges[at - 1].ds) return;
        if (at < m.len) {
            if (pair) {
                if (second.ds > m.edges[at].ds) return;
            } else if (first.ds > m.edges[at].ds) return;
        }
        const n: usize = if (pair) 2 else 1;
        if (m.len + n > max_edges) return;
        @memmove(m.edges[at + n .. m.len + n], m.edges[at..m.len]);
        m.edges[at] = first;
        if (pair) m.edges[at + 1] = second;
        m.len += n;
    }

    /// `cf2_hintmap_adjustHints`: one edge of each unlocked hint onto the
    /// pixel grid, keeping `min_counter` from its neighbors.
    fn adjust(m: *HintMap) void {
        var saved: [max_edges]struct { j: usize, move: i32 } = undefined;
        var n_saved: usize = 0;
        var i: usize = 0;
        while (i < m.len) {
            const pair = m.edges[i].isPair();
            const j = if (pair) i + 1 else i;
            if (!m.edges[i].isLocked()) {
                const frac_down = m.edges[i].ds & 0xFFFF;
                const frac_up = m.edges[j].ds & 0xFFFF;
                const down_move_down = -frac_down;
                const up_move_down = -frac_up;
                const down_move_up: i32 = if (frac_down == 0) 0 else (1 << 16) - frac_down;
                const up_move_up: i32 = if (frac_up == 0) 0 else (1 << 16) - frac_up;
                const move_up = @min(down_move_up, up_move_up);
                const move_down = @max(down_move_down, up_move_down);
                var save = false;
                var adjustment: i32 = 0;
                const room_down = i == 0 or m.edges[i - 1].ds <= m.edges[i].ds +% move_down -% min_counter;
                if (j >= m.len - 1 or m.edges[j + 1].ds >= m.edges[j].ds +% move_up +% min_counter) {
                    adjustment = if (room_down and -move_down < move_up) move_down else move_up;
                } else if (room_down) {
                    adjustment = move_down;
                    save = move_up < -move_down;
                } else {
                    save = true;
                }
                // a move not optimal may be retried once the edges above have moved
                if (save and j < m.len - 1 and !m.edges[j + 1].isLocked()) {
                    saved[n_saved] = .{ .j = j, .move = move_up - adjustment };
                    n_saved += 1;
                }
                m.edges[i].ds +%= adjustment;
                if (pair) m.edges[j].ds +%= adjustment;
            }
            if (i > 0 and m.edges[i].cs != m.edges[i - 1].cs) m.setScale(i - 1);
            if (pair) {
                if (m.edges[j].cs != m.edges[j - 1].cs) m.setScale(j - 1);
                i += 1;
            }
            i += 1;
        }
        var k = n_saved;
        while (k > 0) {
            k -= 1;
            const s = saved[k];
            if (m.edges[s.j + 1].ds >= m.edges[s.j].ds +% s.move +% min_counter) {
                m.edges[s.j].ds +%= s.move;
                if (m.edges[s.j].isPair()) m.edges[s.j - 1].ds +%= s.move;
            }
        }
    }

    /// The scale from edge `i` up to edge `i + 1`.
    fn setScale(m: *HintMap, i: usize) void {
        const a = m.edges[i + 1];
        const b = m.edges[i];
        m.edges[i].scale = fixed.div(a.ds -% b.ds, a.cs -% b.cs);
    }

    /// `cf2_hintmap_build`: the map for the stems `mask` names. Building
    /// the initial map (`initial` null) takes the captured stems only and
    /// locks the baseline; any other build records where it put each stem.
    pub fn build(m: *HintMap, blues: *const Blues, stems: []StemHint, mask_in: HintMask, initial: ?*HintMap, origin: i32) void {
        if (initial) |im| if (!im.valid) im.build(blues, stems, .{}, null, origin);
        const im: ?*const HintMap = initial;
        m.len = 0;
        var mask = mask_in;
        if (blues.em_box) {
            const e = blues.emBoxEdges();
            m.insert(e[0], .{}, im);
            m.insert(.{}, e[1], im);
        }
        // captured or already placed stems first
        for (stems, 0..) |stem, i| {
            if (!mask.get(i)) continue;
            var bottom = Hint.init(stem, i, origin, m.scale, true);
            var top = Hint.init(stem, i, origin, m.scale, false);
            if (bottom.isLocked() or top.isLocked() or blues.capture(&bottom, &top)) {
                m.insert(bottom, top, im);
                mask.clear(i);
            }
        }
        if (initial == null) {
            // an edge at 0 unless one spans it
            if (m.len == 0 or m.edges[0].cs > 0 or m.edges[m.len - 1].cs < 0) {
                m.insert(.{ .flags = Hint.ghost_bottom | Hint.locked | Hint.synthetic, .scale = m.scale }, .{}, null);
            }
        } else for (stems, 0..) |stem, i| {
            if (!mask.get(i)) continue;
            m.insert(Hint.init(stem, i, origin, m.scale, true), Hint.init(stem, i, origin, m.scale, false), im);
        }
        m.adjust();
        if (initial != null) for (m.edges[0..m.len]) |e| {
            if (e.flags & Hint.synthetic != 0) continue;
            const s = &stems[e.index];
            if (e.isTop()) s.ds_max = e.ds else s.ds_min = e.ds;
            s.used = true;
        };
        m.valid = true;
    }
};
