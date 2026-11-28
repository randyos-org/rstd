//! A CFF font's charset, each glyph's SID, as FreeType 2.14.3 reads it
//! (`src/cff/cffload.c` `cff_charset_load`), and the glyph a Standard
//! Encoding code names (`src/psaux/cffdecode.c`
//! `cff_lookup_glyph_by_stdcharcode`), which SEAC accents use.

const std = @import("std");
const CffTable = @import("CffTable.zig");
const data = @import("charset_data.zig");

/// Each glyph's SID in order, read without allocating.
const Sids = struct {
    cff: CffTable,
    format: u8,
    /// The charset offset: 0, 1, 2 for the predefined ones.
    offset: usize,
    count: usize,
    j: usize = 0,
    at: usize = 0,
    sid: u32 = 0,
    left: u32 = 0,

    /// Null if FreeType fails to load the charset (and so the font).
    fn init(cff: CffTable) ?Sids {
        const count = cff.charstrings_index.count();
        const off = cff.charsetOffset() orelse return null;
        var s: Sids = .{ .cff = cff, .format = 0, .offset = off, .count = count };
        switch (off) {
            0 => if (count > 229) return null,
            1 => if (count > data.expert_charset.len) return null,
            2 => if (count > data.expert_subset_charset.len) return null,
            else => {
                const t = cff.table_data;
                if (off >= t.len) return null;
                s.format = t[off];
                s.at = off + 1;
                switch (s.format) {
                    0 => if (count > 0 and s.at + 2 * (count - 1) > t.len) return null,
                    1, 2 => {},
                    else => return null,
                }
            },
        }
        return s;
    }

    /// The next glyph's SID; null after the last, or for a range past the data.
    fn next(s: *Sids) ?u16 {
        if (s.j >= s.count) return null;
        defer s.j += 1;
        switch (s.offset) {
            0 => return @intCast(s.j),
            1 => return data.expert_charset[s.j],
            2 => return data.expert_subset_charset[s.j],
            else => {},
        }
        if (s.j == 0) return 0;
        const t = s.cff.table_data;
        if (s.format == 0) {
            defer s.at += 2;
            return std.mem.readInt(u16, t[s.at..][0..2], .big);
        }
        if (s.left == 0) {
            const wide = s.format == 2;
            if (s.at + @as(usize, if (wide) 4 else 3) > t.len) return null;
            s.sid = std.mem.readInt(u16, t[s.at..][0..2], .big);
            var nleft: u32 = if (wide) std.mem.readInt(u16, t[s.at + 2 ..][0..2], .big) else t[s.at + 2];
            s.at += if (wide) 4 else 3;
            // a range past SID 0xFFFF is trimmed
            if (s.sid > 0xFFFF - nleft) nleft = 0xFFFF - s.sid;
            s.left = nleft + 1;
        }
        defer {
            s.sid += 1;
            s.left -= 1;
        }
        return @intCast(s.sid);
    }
};

/// The first glyph whose SID is Standard Encoding `code`'s; null for a
/// code outside 0..255, or no such glyph.
pub fn glyphForStdCode(cff: CffTable, code: i32) ?u16 {
    if (code < 0 or code > 255) return null;
    const sid = data.standard_encoding[@intCast(code)];
    var s = Sids.init(cff) orelse return null;
    while (s.next()) |v| {
        if (v == sid) return @intCast(s.j - 1);
    }
    return null;
}
