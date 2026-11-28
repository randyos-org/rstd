//! FreeType's tricky fonts: old CJK faces whose glyphs are built by their
//! bytecode, which must never be autohinted. FreeType 2.14.3
//! `src/truetype/ttobjs.c` (`tt_skip_pdffont_random_tag`,
//! `tt_check_trickyness_family`, `tt_synth_sfnt_checksum`,
//! `tt_check_trickyness_sfnt_ids`, `tt_check_trickyness`), the family name
//! as `src/sfnt/sfobjs.c` (`sfnt_load_face`, `tt_face_get_name`) picks it.

const std = @import("std");
const Sfnt = @import("../sfnt/Sfnt.zig");

/// `trick_names`; a `?` stands for a non-ASCII character.
const names = [_][]const u8{
    "cpop",               "DFGirl-W6-WIN-BF", "DFGothic-EB",   "DFGyoSho-Lt",   "DFHei",
    "DFHSGothic-W5",      "DFHSMincho-W3",    "DFHSMincho-W7", "DFKaiSho-SB",   "DFKaiShu",
    "DFKai-SB",           "DFMing",           "DLC",           "HuaTianKaiTi?", "HuaTianSongTi?",
    "Ming(for ISO10646)", "MingLiU",          "MingMedium",    "PMingLiU",      "MingLi43",
};

/// Per face, the checksum and length of its `cvt `, `fpgm` and `prep`
/// tables (`sfnt_id`); data from that file, regenerated rather than edited.
const ids = [_][3][2]u32{
    // MingLiU 1995
    .{ .{ 0x05BCF058, 0x000002E4 }, .{ 0x28233BF1, 0x000087C4 }, .{ 0xA344A1EA, 0x000001E1 } },
    // MingLiU 1996-
    .{ .{ 0x05BCF058, 0x000002E4 }, .{ 0x28233BF1, 0x000087C4 }, .{ 0xA344A1EB, 0x000001E1 } },
    // DFGothic-EB
    .{ .{ 0x12C3EBB2, 0x00000350 }, .{ 0xB680EE64, 0x000087A7 }, .{ 0xCE939563, 0x00000758 } },
    // DFGyoSho-Lt
    .{ .{ 0x11E5EAD4, 0x00000350 }, .{ 0xCE5956E9, 0x0000BC85 }, .{ 0x8272F416, 0x00000045 } },
    // DFHei-Md-HK-BF
    .{ .{ 0x1257EB46, 0x00000350 }, .{ 0xF699D160, 0x0000715F }, .{ 0xD222F568, 0x000003BC } },
    // DFHSGothic-W5
    .{ .{ 0x1262EB4E, 0x00000350 }, .{ 0xE86A5D64, 0x00007940 }, .{ 0x7850F729, 0x000005FF } },
    // DFHSMincho-W3
    .{ .{ 0x122DEB0A, 0x00000350 }, .{ 0x3D16328A, 0x0000859B }, .{ 0xA93FC33B, 0x000002CB } },
    // DFHSMincho-W7
    .{ .{ 0x125FEB26, 0x00000350 }, .{ 0xA5ACC982, 0x00007EE1 }, .{ 0x90999196, 0x0000041F } },
    // DFKaiShu
    .{ .{ 0x11E5EAD4, 0x00000350 }, .{ 0x5A30CA3B, 0x00009063 }, .{ 0x13A42602, 0x0000007E } },
    // DFKaiShu, variant
    .{ .{ 0x11E5EAD4, 0x00000350 }, .{ 0xA6E78C01, 0x00008998 }, .{ 0x13A42602, 0x0000007E } },
    // DFKaiShu-Md-HK-BF
    .{ .{ 0x11E5EAD4, 0x00000360 }, .{ 0x9DB282B2, 0x0000C06E }, .{ 0x53E6D7CA, 0x00000082 } },
    // DFMing-Bd-HK-BF
    .{ .{ 0x1243EB18, 0x00000350 }, .{ 0xBA0A8C30, 0x000074AD }, .{ 0xF3D83409, 0x0000037B } },
    // DLCLiShu
    .{ .{ 0x07DCF546, 0x00000308 }, .{ 0x40FE7C90, 0x00008E2A }, .{ 0x608174B5, 0x0000007A } },
    // DLCHayBold
    .{ .{ 0xEB891238, 0x00000308 }, .{ 0xD2E4DCD4, 0x0000676F }, .{ 0x8EA5F293, 0x000003B8 } },
    // HuaTianKaiTi
    .{ .{ 0xFFFBFFFC, 0x00000008 }, .{ 0x9C9E48B8, 0x0000BEA2 }, .{ 0x70020112, 0x00000008 } },
    // HuaTianSongTi
    .{ .{ 0xFFFBFFFC, 0x00000008 }, .{ 0x0A5A0483, 0x00017C39 }, .{ 0x70020112, 0x00000008 } },
    // NEC fadpop7.ttf
    .{ .{ 0x00000000, 0x00000000 }, .{ 0x40C92555, 0x000000E5 }, .{ 0xA39B58E3, 0x0000117C } },
    // NEC fadrei5.ttf
    .{ .{ 0x00000000, 0x00000000 }, .{ 0x33C41652, 0x000000E5 }, .{ 0x26D6C52A, 0x00000F6A } },
    // NEC fangot7.ttf
    .{ .{ 0x00000000, 0x00000000 }, .{ 0x6DB1651D, 0x0000019D }, .{ 0x6C6E4B03, 0x00002492 } },
    // NEC fangyo5.ttf
    .{ .{ 0x00000000, 0x00000000 }, .{ 0x40C92555, 0x000000E5 }, .{ 0xDE51FAD0, 0x0000117C } },
    // NEC fankyo5.ttf
    .{ .{ 0x00000000, 0x00000000 }, .{ 0x85E47664, 0x000000E5 }, .{ 0xA6C62831, 0x00001CAA } },
    // NEC fanrgo5.ttf
    .{ .{ 0x00000000, 0x00000000 }, .{ 0x2D891CFD, 0x0000019D }, .{ 0xA0604633, 0x00001DE8 } },
    // NEC fangot5.ttc
    .{ .{ 0x00000000, 0x00000000 }, .{ 0x40AA774C, 0x000001CB }, .{ 0x9B5CAA96, 0x00001F9A } },
    // NEC fanmin3.ttc
    .{ .{ 0x00000000, 0x00000000 }, .{ 0x0D3DE9CB, 0x00000141 }, .{ 0xD4127766, 0x00002280 } },
    // NEC FA-Gothic, 1996
    .{ .{ 0x00000000, 0x00000000 }, .{ 0x4A692698, 0x000001F0 }, .{ 0x340D4346, 0x00001FCA } },
    // NEC FA-Minchou, 1996
    .{ .{ 0x00000000, 0x00000000 }, .{ 0xCD34C604, 0x00000166 }, .{ 0x6CF31046, 0x000022B0 } },
    // NEC FA-RoundGothicB, 1996
    .{ .{ 0x00000000, 0x00000000 }, .{ 0x5DA75315, 0x0000019D }, .{ 0x40745A5F, 0x000022E0 } },
    // NEC FA-RoundGothicM, 1996
    .{ .{ 0x00000000, 0x00000000 }, .{ 0xF055FC48, 0x000001C2 }, .{ 0x3900DED3, 0x00001E18 } },
    // MINGLI.TTF, 1992
    .{ .{ 0x00170003, 0x00000060 }, .{ 0xDBB4306E, 0x000058AA }, .{ 0xD643482A, 0x00000035 } },
    // DFHei-Bd-WIN-HK-BF, issue #1087
    .{ .{ 0x1269EB58, 0x00000350 }, .{ 0x5CD5957A, 0x00006A4E }, .{ 0xF758323A, 0x00000380 } },
    // DFMing-Md-WIN-HK-BF, issue #1087
    .{ .{ 0x122FEB0B, 0x00000350 }, .{ 0x7F10919A, 0x000070A9 }, .{ 0x7CD7E7B7, 0x0000025C } },
};

/// `tt_check_trickyness`.
pub fn isTricky(sfnt: *const Sfnt) bool {
    // only the TrueType driver checks
    if (sfnt.glyf == null) return false;
    var buf: [256]u8 = undefined;
    if (familyName(sfnt, &buf)) |name| {
        // `tt_skip_pdffont_random_tag`: a PDF subset's `ABCDEF+`
        var n = name;
        if (n.len > 7 and n[6] == '+') {
            for (n[0..6]) |c| {
                if (c < 'A' or c > 'Z') break;
            } else n = n[7..];
        }
        for (names) |t| if (std.mem.indexOf(u8, n, t) != null) return true;
    }
    return checksumsMatch(sfnt);
}

/// `tt_synth_sfnt_checksum`.
pub fn checksum(data: []const u8) u32 {
    var sum: u32 = 0;
    var i: usize = 0;
    while (i + 4 <= data.len) : (i += 4) sum +%= std.mem.readInt(u32, data[i..][0..4], .big);
    var shift: u5 = 24;
    while (i < data.len) : (i += 1) {
        sum +%= @as(u32, data[i]) << shift;
        shift -%= 8;
    }
    return sum;
}

/// `tt_check_trickyness_sfnt_ids`.
fn checksumsMatch(sfnt: *const Sfnt) bool {
    const tables = [3]?[]const u8{ sfnt.cvt, sfnt.fpgm, sfnt.prep };
    var matched: [ids.len]u8 = @splat(0);
    for (tables, 0..) |t, k| {
        const data = t orelse continue;
        const sum = checksum(data);
        for (ids, 0..) |id, j| {
            if (data.len != id[k][1]) continue;
            if (id[k][0] == sum) matched[j] += 1;
            if (matched[j] == 3) return true;
        }
    }
    // a table a face lacks matches an entry of length 0
    for (ids, 0..) |id, j| {
        for (tables, 0..) |t, k| {
            if (t == null and id[k][1] == 0) matched[j] += 1;
        }
        if (matched[j] == 3) return true;
    }
    return false;
}

/// The family name, ASCII with `?` for the rest, as FreeType's
/// `family_name`; null without one.
fn familyName(sfnt: *const Sfnt, buf: []u8) ?[]const u8 {
    const wws = if (sfnt.os2) |os2| @as(u16, @bitCast(os2.fsSelection())) & 256 != 0 else false;
    const order: []const u16 = if (wws) &.{ 16, 1 } else &.{ 21, 16, 1 };
    for (order) |id| if (nameString(sfnt, id, buf)) |s| return s;
    return null;
}

/// `tt_face_get_name`: Windows Unicode (English preferred) unless only a
/// non-English one and a Macintosh one exist, then Macintosh (English
/// preferred, else Roman), then Unicode or ISO.
fn nameString(sfnt: *const Sfnt, name_id: u16, buf: []u8) ?[]const u8 {
    const table = sfnt.name orelse return null;
    var apple_roman: ?usize = null;
    var apple_english: ?usize = null;
    var win: ?usize = null;
    var unicode: ?usize = null;
    var is_english = false;
    for (0..table.count()) |n| {
        const rec = table.getNameRecord(n);
        if (rec.name_id != name_id or rec.length == 0) continue;
        switch (rec.platform_id) {
            0, 2 => unicode = n,
            1 => {
                if (rec.language_id == 0) {
                    apple_english = n;
                } else if (rec.encoding_id == 0) apple_roman = n;
            },
            3 => if (win == null or rec.language_id & 0x3FF == 0x009) switch (rec.encoding_id) {
                0, 1, 10 => {
                    is_english = rec.language_id & 0x3FF == 0x009;
                    win = n;
                },
                else => {},
            },
            else => {},
        }
    }
    const apple = apple_english orelse apple_roman;
    const pick: struct { n: usize, utf16: bool } = if (win != null and !(apple != null and !is_english))
        .{ .n = win.?, .utf16 = true }
    else if (apple) |a|
        .{ .n = a, .utf16 = false }
    else if (unicode) |u|
        .{ .n = u, .utf16 = true }
    else
        return null;
    const raw = table.getRawString(table.getNameRecord(pick.n)) orelse return null;
    var len: usize = 0;
    var i: usize = 0;
    while (len < buf.len) {
        var code: u16 = undefined;
        if (pick.utf16) {
            if (i + 2 > raw.len) break;
            code = std.mem.readInt(u16, raw[i..][0..2], .big);
            i += 2;
        } else {
            if (i >= raw.len) break;
            code = raw[i];
            i += 1;
        }
        if (code == 0) break;
        buf[len] = if (code < 32 or code > 127) '?' else @intCast(code);
        len += 1;
    }
    return buf[0..len];
}
