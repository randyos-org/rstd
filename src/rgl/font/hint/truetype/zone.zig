//! The interpreter's zones: the twilight zone's scratch points and the
//! glyph's (FreeType 2.14.3 `ttobjs.h` `TT_GlyphZoneRec`; fontations
//! `skrifa/src/outline/glyf/hint/zone.rs`, commit 4f115482, MIT).

/// A point in 26.6 pixels, or in font units for `Zone.orus`.
pub const Point = [2]i32;

/// Which zone a zone pointer selects (`SZP0` and the rest).
pub const Pointer = enum(u1) { twilight = 0, glyph = 1 };

/// FreeType's point tags: on the curve, and touched in x or y.
pub const on_curve: u8 = 0x01;
pub const touched_x: u8 = 0x08;
pub const touched_y: u8 = 0x10;

pub const Zone = struct {
    /// The points in font units; empty for the twilight zone, whose are 0.
    orus: []const Point = &.{},
    /// The points as scaled, before any instruction.
    org: []Point = &.{},
    /// The points as instructions move them.
    cur: []Point = &.{},
    tags: []u8 = &.{},
    /// Each contour's last point.
    contours: []const u16 = &.{},

    pub fn len(z: Zone) usize {
        return z.cur.len;
    }

    /// Point `i` in font units; (0, 0) past `orus`.
    pub fn unscaled(z: Zone, i: usize) Point {
        return if (i < z.orus.len) z.orus[i] else .{ 0, 0 };
    }
};
