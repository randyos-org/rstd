//! The autohinter's character clusters without HarfBuzz: FreeType 2.14.3
//! `src/autofit/afshaper.c` (`af_shaper_get_cluster_nohb`) and `afblue.h`
//! (`GET_UTF8_CHAR`). FreeType loads HarfBuzz at run time when it can; the
//! references here were made without it.

const Sfnt = @import("../../sfnt/Sfnt.zig");

/// `GET_UTF8_CHAR`: no validation, as FreeType's.
fn utf8(s: []const u8, at: *usize) u32 {
    var ch: u32 = s[at.*];
    at.* += 1;
    if (ch >= 0x80) {
        var len: u32 = undefined;
        if (ch < 0xE0) {
            len = 1;
            ch &= 0x1F;
        } else if (ch < 0xF0) {
            len = 2;
            ch &= 0x0F;
        } else {
            len = 3;
            ch &= 0x07;
        }
        while (len > 0 and at.* < s.len) : (len -= 1) {
            ch = (ch << 6) | (s[at.*] & 0x3F);
            at.* += 1;
        }
    }
    return ch;
}

pub const Cluster = struct {
    /// Glyphs: 0 for a cluster of several characters, else 1.
    count: u32,
    /// The glyph, 0 if unmapped.
    glyph: u16,
};

/// The cluster at `at` in the space-separated `s`, `at` moved past it.
pub fn cluster(sfnt: *const Sfnt, s: []const u8, at: *usize) Cluster {
    while (at.* < s.len and s[at.*] == ' ') at.* += 1;
    if (at.* >= s.len) return .{ .count = 0, .glyph = 0 };
    const ch = utf8(s, at);
    var more = false;
    while (at.* < s.len and s[at.*] != ' ') {
        _ = utf8(s, at);
        more = true;
    }
    if (more) return .{ .count = 0, .glyph = 0 };
    const g = sfnt.getGlyphIndex(ch) orelse return .{ .count = 1, .glyph = 0 };
    return .{ .count = 1, .glyph = @backingInt(g) };
}
