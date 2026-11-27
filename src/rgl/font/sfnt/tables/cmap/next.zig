//! The next mapped code after a given one, as FreeType's `FT_Get_Next_Char`
//! walks a charmap: codes in increasing order, skipping those mapped to
//! glyph 0.

const types = @import("../../types.zig");
const CmapSubtable = @import("__root__.zig").CmapSubtable;

const GlyphIndex = types.GlyphIndex;

pub const Mapping = struct { codepoint: u32, glyph: GlyphIndex };

/// The first code above `after` mapped to a glyph other than 0, if any.
pub fn nextChar(sub: CmapSubtable, after: u32) ?Mapping {
    if (after >= 0x10FFFF) return null;
    const from = after + 1;
    switch (sub) {
        .fmt4 => |t| {
            // the first segment ending at or past `from`; segments are sorted
            var lo: usize = 0;
            var hi: usize = t.seg_count;
            while (lo < hi) {
                const mid = lo + (hi - lo) / 2;
                if (t.getEndCode(mid) < from) lo = mid + 1 else hi = mid;
            }
            for (lo..t.seg_count) |i| {
                var cp: u32 = @max(from, t.getStartCode(i));
                while (cp <= t.getEndCode(i)) : (cp += 1) {
                    if (scan(sub, cp)) |m| return m;
                }
            }
            return null;
        },
        inline .fmt8, .fmt12, .fmt13 => |t| {
            var lo: usize = 0;
            var hi: usize = t.numGroups();
            while (lo < hi) {
                const mid = lo + (hi - lo) / 2;
                if (t.getGroup(mid).end_char_code < from) lo = mid + 1 else hi = mid;
            }
            for (lo..t.numGroups()) |i| {
                const g = t.getGroup(i);
                var cp: u32 = @max(from, g.start_char_code);
                while (cp <= g.end_char_code) : (cp += 1) {
                    if (scan(sub, cp)) |m| return m;
                }
            }
            return null;
        },
        inline .fmt0, .fmt6, .fmt10 => |t| {
            var cp: u32 = @max(from, t.first());
            while (cp < t.first() + t.count()) : (cp += 1) {
                if (scan(sub, cp)) |m| return m;
            }
            return null;
        },
        .fmt2 => {
            var cp = from;
            while (cp <= 0xFFFF) : (cp += 1) {
                if (scan(sub, cp)) |m| return m;
            }
            return null;
        },
        .fmt14, .unsupported => return null,
    }
}

fn scan(sub: CmapSubtable, cp: u32) ?Mapping {
    const g = sub.getGlyphIndex(cp);
    return if (g == .notdef) null else .{ .codepoint = cp, .glyph = g };
}
