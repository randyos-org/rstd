//! The autohinter's scripts: each one's Unicode ranges, the ranges of its non-base characters, its hinting order and the standard characters its stem widths come from. Data from FreeType 2.14.3 (tag `VER-2-14-3`) `src/autofit/afscript.h` and `afranges.c`;
//! regenerate it from those files with every macro expanded (`AF_CONFIG_OPTION_CJK` and `AF_CONFIG_OPTION_INDIC` on, as FreeType builds by default) rather than editing it.

pub const Script = struct {
    tag: []const u8,
    /// Vertical hinting goes from the top down (`AF_HINTING_TOP_TO_BOTTOM`).
    top_to_bottom: bool,
    standard_chars: []const u8,
    ranges: []const [2]u32,
    nonbase: []const [2]u32,
};

/// In `af_script_classes` order.
pub const scripts = [_]Script{
    .{ .tag = "adlm", .top_to_bottom = false, .standard_chars = "\xF0\x9E\xA4\x8C \xF0\x9E\xA4\xAE", .ranges = &.{
        .{ 0x1E900, 0x1E95F },
    }, .nonbase = &.{
        .{ 0x1D944, 0x1E94A },
    } },
    .{ .tag = "arab", .top_to_bottom = false, .standard_chars = "\xD9\x84 \xD8\xAD \xD9\x80", .ranges = &.{
        .{ 0x600, 0x6FF },     .{ 0x750, 0x7FF },     .{ 0x870, 0x89F }, .{ 0x8A0, 0x8FF }, .{ 0xFB50, 0xFDFF }, .{ 0xFE70, 0xFEFF },
        .{ 0x10EC0, 0x10EFF }, .{ 0x1EE00, 0x1EEFF },
    }, .nonbase = &.{
        .{ 0x600, 0x605 },   .{ 0x610, 0x61A },   .{ 0x64B, 0x65F },     .{ 0x670, 0x670 },   .{ 0x6D6, 0x6DC },   .{ 0x6DF, 0x6E4 },
        .{ 0x6E7, 0x6E8 },   .{ 0x6EA, 0x6ED },   .{ 0x897, 0x89F },     .{ 0x8CA, 0x8E1 },   .{ 0x8E3, 0x8FF },   .{ 0xFBB2, 0xFBC1 },
        .{ 0xFE70, 0xFE70 }, .{ 0xFE72, 0xFE72 }, .{ 0xFE74, 0xFE74 },   .{ 0xFE76, 0xFE76 }, .{ 0xFE78, 0xFE78 }, .{ 0xFE7A, 0xFE7A },
        .{ 0xFE7C, 0xFE7C }, .{ 0xFE7E, 0xFE7E }, .{ 0x10EFD, 0x10EFF },
    } },
    .{ .tag = "armn", .top_to_bottom = false, .standard_chars = "\xD5\xBD \xD5\x8D", .ranges = &.{
        .{ 0x530, 0x58F }, .{ 0xFB13, 0xFB17 },
    }, .nonbase = &.{
        .{ 0x559, 0x55F },
    } },
    .{ .tag = "avst", .top_to_bottom = false, .standard_chars = "\xF0\x90\xAC\x9A", .ranges = &.{
        .{ 0x10B00, 0x10B3F },
    }, .nonbase = &.{
        .{ 0x10B39, 0x10B3F },
    } },
    .{ .tag = "bamu", .top_to_bottom = false, .standard_chars = "\xEA\x9B\x81 \xEA\x9B\xAF", .ranges = &.{
        .{ 0xA6A0, 0xA6FF },
    }, .nonbase = &.{
        .{ 0xA6F0, 0xA6F1 },
    } },
    .{ .tag = "beng", .top_to_bottom = true, .standard_chars = "\xE0\xA7\xA6 \xE0\xA7\xAA", .ranges = &.{
        .{ 0x980, 0x9FF },
    }, .nonbase = &.{
        .{ 0x981, 0x981 }, .{ 0x9BC, 0x9BC }, .{ 0x9C1, 0x9C4 }, .{ 0x9CD, 0x9CD }, .{ 0x9E2, 0x9E3 }, .{ 0x9FE, 0x9FE },
    } },
    .{ .tag = "buhd", .top_to_bottom = false, .standard_chars = "\xE1\x9D\x8B \xE1\x9D\x8F", .ranges = &.{
        .{ 0x1740, 0x175F },
    }, .nonbase = &.{
        .{ 0x1752, 0x1753 },
    } },
    .{ .tag = "cakm", .top_to_bottom = false, .standard_chars = "\xF0\x91\x84\xA4 \xF0\x91\x84\x89 \xF0\x91\x84\x9B", .ranges = &.{
        .{ 0x11100, 0x1114F },
    }, .nonbase = &.{
        .{ 0x11100, 0x11102 }, .{ 0x11127, 0x11134 }, .{ 0x11146, 0x11146 },
    } },
    .{ .tag = "cans", .top_to_bottom = false, .standard_chars = "\xE1\x91\x8C \xE1\x93\x9A", .ranges = &.{
        .{ 0x1400, 0x167F }, .{ 0x18B0, 0x18FF }, .{ 0x11AB0, 0x11ABF },
    }, .nonbase = &.{} },
    .{ .tag = "cari", .top_to_bottom = false, .standard_chars = "\xF0\x90\x8A\xAB \xF0\x90\x8B\x89", .ranges = &.{
        .{ 0x102A0, 0x102DF },
    }, .nonbase = &.{} },
    .{ .tag = "cher", .top_to_bottom = false, .standard_chars = "\xE1\x8E\xA4 \xE1\x8F\x85 \xEA\xAE\x95", .ranges = &.{
        .{ 0x13A0, 0x13FF }, .{ 0xAB70, 0xABBF },
    }, .nonbase = &.{} },
    .{ .tag = "copt", .top_to_bottom = false, .standard_chars = "\xE2\xB2\x9E \xE2\xB2\x9F", .ranges = &.{
        .{ 0x2C80, 0x2CFF },
    }, .nonbase = &.{
        .{ 0x2CEF, 0x2CF1 },
    } },
    .{ .tag = "cprt", .top_to_bottom = false, .standard_chars = "\xF0\x90\xA0\x85 \xF0\x90\xA0\xA3", .ranges = &.{
        .{ 0x10800, 0x1083F },
    }, .nonbase = &.{} },
    .{ .tag = "cyrl", .top_to_bottom = false, .standard_chars = "\xD0\xBE \xD0\x9E", .ranges = &.{
        .{ 0x400, 0x4FF }, .{ 0x500, 0x52F }, .{ 0x2DE0, 0x2DFF }, .{ 0xA640, 0xA69F }, .{ 0x1C80, 0x1C8F }, .{ 0x1E030, 0x1E08F },
    }, .nonbase = &.{
        .{ 0x483, 0x489 }, .{ 0x2DE0, 0x2DFF }, .{ 0xA66F, 0xA67F }, .{ 0xA69E, 0xA69F },
    } },
    .{ .tag = "deva", .top_to_bottom = true, .standard_chars = "\xE0\xA4\xA0 \xE0\xA4\xB5 \xE0\xA4\x9F", .ranges = &.{
        .{ 0x900, 0x93B },     .{ 0x93D, 0x950 }, .{ 0x953, 0x963 }, .{ 0x966, 0x97F }, .{ 0x20B9, 0x20B9 }, .{ 0xA8E0, 0xA8FF },
        .{ 0x11B00, 0x11B5F },
    }, .nonbase = &.{
        .{ 0x900, 0x902 },   .{ 0x93A, 0x93A },   .{ 0x941, 0x948 }, .{ 0x94D, 0x94D }, .{ 0x953, 0x957 }, .{ 0x962, 0x963 },
        .{ 0xA8E0, 0xA8F1 }, .{ 0xA8FF, 0xA8FF },
    } },
    .{ .tag = "dsrt", .top_to_bottom = false, .standard_chars = "\xF0\x90\x90\x84 \xF0\x90\x90\xAC", .ranges = &.{
        .{ 0x10400, 0x1044F },
    }, .nonbase = &.{} },
    .{ .tag = "ethi", .top_to_bottom = false, .standard_chars = "\xE1\x8B\x90", .ranges = &.{
        .{ 0x1200, 0x137F }, .{ 0x1380, 0x139F }, .{ 0x2D80, 0x2DDF }, .{ 0xAB00, 0xAB2F }, .{ 0x1E7E0, 0x1E7FF },
    }, .nonbase = &.{
        .{ 0x135D, 0x135F },
    } },
    .{ .tag = "geor", .top_to_bottom = false, .standard_chars = "\xE1\x83\x98 \xE1\x83\x94 \xE1\x83\x90 \xE1\xB2\xBF", .ranges = &.{
        .{ 0x10D0, 0x10FF }, .{ 0x1C90, 0x1CBF },
    }, .nonbase = &.{} },
    .{ .tag = "geok", .top_to_bottom = false, .standard_chars = "\xE1\x82\xB6 \xE1\x82\xB1 \xE2\xB4\x99", .ranges = &.{
        .{ 0x10A0, 0x10CD }, .{ 0x2D00, 0x2D2D },
    }, .nonbase = &.{} },
    .{ .tag = "glag", .top_to_bottom = false, .standard_chars = "\xE2\xB0\x95 \xE2\xB1\x85", .ranges = &.{
        .{ 0x2C00, 0x2C5F }, .{ 0x1E000, 0x1E02F },
    }, .nonbase = &.{
        .{ 0x1E000, 0x1E02F },
    } },
    .{ .tag = "goth", .top_to_bottom = true, .standard_chars = "\xF0\x90\x8C\xB4 \xF0\x90\x8C\xBE \xF0\x90\x8D\x83", .ranges = &.{
        .{ 0x10330, 0x1034F },
    }, .nonbase = &.{} },
    .{ .tag = "grek", .top_to_bottom = false, .standard_chars = "\xCE\xBF \xCE\x9F", .ranges = &.{
        .{ 0x370, 0x3FF }, .{ 0x1F00, 0x1FFF },
    }, .nonbase = &.{
        .{ 0x37A, 0x37A },   .{ 0x384, 0x385 }, .{ 0x1FBD, 0x1FC1 }, .{ 0x1FCD, 0x1FCF }, .{ 0x1FDD, 0x1FDF }, .{ 0x1FED, 0x1FEF },
        .{ 0x1FFD, 0x1FFE },
    } },
    .{ .tag = "gujr", .top_to_bottom = false, .standard_chars = "\xE0\xAA\x9F \xE0\xAB\xA6", .ranges = &.{
        .{ 0xA80, 0xAFF },
    }, .nonbase = &.{
        .{ 0xA81, 0xA82 }, .{ 0xABC, 0xABC }, .{ 0xAC1, 0xAC8 }, .{ 0xACD, 0xACD }, .{ 0xAE2, 0xAE3 }, .{ 0xAFA, 0xAFF },
    } },
    .{ .tag = "guru", .top_to_bottom = true, .standard_chars = "\xE0\xA8\xA0 \xE0\xA8\xB0 \xE0\xA9\xA6", .ranges = &.{
        .{ 0xA00, 0xA7F },
    }, .nonbase = &.{
        .{ 0xA01, 0xA02 }, .{ 0xA3C, 0xA3C }, .{ 0xA41, 0xA51 }, .{ 0xA70, 0xA71 }, .{ 0xA75, 0xA75 },
    } },
    .{ .tag = "hebr", .top_to_bottom = false, .standard_chars = "\xD7\x9D", .ranges = &.{
        .{ 0x590, 0x5FF }, .{ 0xFB1D, 0xFB4F },
    }, .nonbase = &.{
        .{ 0x591, 0x5BF }, .{ 0x5C1, 0x5C2 }, .{ 0x5C4, 0x5C5 }, .{ 0x5C7, 0x5C7 }, .{ 0xFB1E, 0xFB1E },
    } },
    .{ .tag = "kali", .top_to_bottom = false, .standard_chars = "\xEA\xA4\x8D \xEA\xA4\x80", .ranges = &.{
        .{ 0xA900, 0xA92F },
    }, .nonbase = &.{
        .{ 0xA926, 0xA92D },
    } },
    .{ .tag = "khmr", .top_to_bottom = false, .standard_chars = "\xE1\x9F\xA0", .ranges = &.{
        .{ 0x1780, 0x17FF },
    }, .nonbase = &.{
        .{ 0x17B7, 0x17BD }, .{ 0x17C6, 0x17C6 }, .{ 0x17C9, 0x17D3 }, .{ 0x17DD, 0x17DD },
    } },
    .{ .tag = "khms", .top_to_bottom = false, .standard_chars = "\xE1\xA7\xA1 \xE1\xA7\xAA", .ranges = &.{
        .{ 0x19E0, 0x19FF },
    }, .nonbase = &.{} },
    .{ .tag = "knda", .top_to_bottom = false, .standard_chars = "\xE0\xB3\xA6 \xE0\xB2\xAC", .ranges = &.{
        .{ 0xC80, 0xCFF },
    }, .nonbase = &.{
        .{ 0xC81, 0xC81 }, .{ 0xCBC, 0xCBC }, .{ 0xCBF, 0xCBF }, .{ 0xCC6, 0xCC6 }, .{ 0xCCC, 0xCCD }, .{ 0xCE2, 0xCE3 },
    } },
    .{ .tag = "lao", .top_to_bottom = false, .standard_chars = "\xE0\xBB\x90", .ranges = &.{
        .{ 0xE80, 0xEFF },
    }, .nonbase = &.{
        .{ 0xEB1, 0xEB1 }, .{ 0xEB4, 0xEBC }, .{ 0xEC8, 0xECE },
    } },
    .{ .tag = "latn", .top_to_bottom = false, .standard_chars = "o O 0", .ranges = &.{
        .{ 0x20, 0x7F },     .{ 0xA0, 0xA9 },     .{ 0xAB, 0xB1 },     .{ 0xB4, 0xB8 },       .{ 0xBB, 0xFF },       .{ 0x100, 0x17F },
        .{ 0x180, 0x24F },   .{ 0x250, 0x2AF },   .{ 0x2B9, 0x2DF },   .{ 0x2E5, 0x2FF },     .{ 0x300, 0x36F },     .{ 0x1AB0, 0x1ABE },
        .{ 0x1D00, 0x1D2B }, .{ 0x1D6B, 0x1D77 }, .{ 0x1D79, 0x1D7F }, .{ 0x1D80, 0x1D9A },   .{ 0x1DC0, 0x1DFF },   .{ 0x1E00, 0x1EFF },
        .{ 0x2000, 0x206F }, .{ 0x20A0, 0x20B8 }, .{ 0x20BA, 0x20CF }, .{ 0x2150, 0x218F },   .{ 0x2C60, 0x2C7B },   .{ 0x2C7E, 0x2C7F },
        .{ 0x2E00, 0x2E7F }, .{ 0xA720, 0xA76F }, .{ 0xA771, 0xA7F0 }, .{ 0xA7F2, 0xA7F7 },   .{ 0xA7FA, 0xA7FF },   .{ 0xAB30, 0xAB5B },
        .{ 0xAB60, 0xAB68 }, .{ 0xAB6A, 0xAB6F }, .{ 0xFB00, 0xFB06 }, .{ 0x1D400, 0x1D7FF }, .{ 0x1DF00, 0x1DFFF },
    }, .nonbase = &.{
        .{ 0x5E, 0x60 },     .{ 0x7E, 0x7E },     .{ 0xA8, 0xA9 },     .{ 0xAE, 0xB0 },     .{ 0xB4, 0xB4 },     .{ 0xB8, 0xB8 },
        .{ 0xBC, 0xBE },     .{ 0x2B9, 0x2DF },   .{ 0x2E5, 0x2FF },   .{ 0x300, 0x36F },   .{ 0x1AB0, 0x1AEB }, .{ 0x1DC0, 0x1DFF },
        .{ 0x2017, 0x2017 }, .{ 0x203E, 0x203E }, .{ 0xA788, 0xA788 }, .{ 0xA7F8, 0xA7FA },
    } },
    .{ .tag = "latb", .top_to_bottom = false, .standard_chars = "\xE2\x82\x92 \xE2\x82\x80", .ranges = &.{
        .{ 0x1D62, 0x1D6A }, .{ 0x2080, 0x209C }, .{ 0x2C7C, 0x2C7C },
    }, .nonbase = &.{} },
    .{ .tag = "latp", .top_to_bottom = false, .standard_chars = "\xE1\xB5\x92 \xE1\xB4\xBC \xE2\x81\xB0", .ranges = &.{
        .{ 0xAA, 0xAA },     .{ 0xB2, 0xB3 },     .{ 0xB9, 0xBA },     .{ 0x2B0, 0x2B8 },     .{ 0x2E0, 0x2E4 },   .{ 0x1D2C, 0x1D61 },
        .{ 0x1D78, 0x1D78 }, .{ 0x1D9B, 0x1DBF }, .{ 0x2070, 0x207F }, .{ 0x2C7D, 0x2C7D },   .{ 0xA770, 0xA770 }, .{ 0xA7F1, 0xA7F1 },
        .{ 0xA7F8, 0xA7F9 }, .{ 0xAB5C, 0xAB5F }, .{ 0xAB69, 0xAB69 }, .{ 0x10780, 0x107FB },
    }, .nonbase = &.{} },
    .{ .tag = "lisu", .top_to_bottom = false, .standard_chars = "\xEA\x93\xB3", .ranges = &.{
        .{ 0xA4D0, 0xA4FF }, .{ 0x11FB0, 0x11FBF },
    }, .nonbase = &.{} },
    .{ .tag = "mlym", .top_to_bottom = false, .standard_chars = "\xE0\xB4\xA0 \xE0\xB4\xB1", .ranges = &.{
        .{ 0xD00, 0xD7F },
    }, .nonbase = &.{
        .{ 0xD00, 0xD01 }, .{ 0xD3B, 0xD3C }, .{ 0xD4D, 0xD4E }, .{ 0xD62, 0xD63 },
    } },
    .{ .tag = "medf", .top_to_bottom = false, .standard_chars = "\xF0\x96\xB9\xA1 \xF0\x96\xB9\x9B \xF0\x96\xB9\xAF", .ranges = &.{
        .{ 0x16E40, 0x16E9F },
    }, .nonbase = &.{} },
    .{ .tag = "mong", .top_to_bottom = true, .standard_chars = "\xE1\xA1\x82 \xE1\xA0\xAA", .ranges = &.{
        .{ 0x1800, 0x18AF }, .{ 0x11660, 0x1167F },
    }, .nonbase = &.{
        .{ 0x1885, 0x1886 }, .{ 0x18A9, 0x18A9 },
    } },
    .{ .tag = "mymr", .top_to_bottom = false, .standard_chars = "\xE1\x80\x9D \xE1\x80\x84 \xE1\x80\x82", .ranges = &.{
        .{ 0x1000, 0x109F }, .{ 0xA9E0, 0xA9FF }, .{ 0xAA60, 0xAA7F }, .{ 0x116D0, 0x116FF },
    }, .nonbase = &.{
        .{ 0x102D, 0x1030 }, .{ 0x1032, 0x1037 }, .{ 0x103A, 0x103A }, .{ 0x103D, 0x103E }, .{ 0x1058, 0x1059 }, .{ 0x105E, 0x1060 },
        .{ 0x1071, 0x1074 }, .{ 0x1082, 0x1082 }, .{ 0x1085, 0x1086 }, .{ 0x108D, 0x108D }, .{ 0xA9E5, 0xA9E5 }, .{ 0xAA7C, 0xAA7C },
    } },
    .{ .tag = "nkoo", .top_to_bottom = false, .standard_chars = "\xDF\x8B \xDF\x80", .ranges = &.{
        .{ 0x7C0, 0x7FF },
    }, .nonbase = &.{
        .{ 0x7EB, 0x7F5 }, .{ 0x7FD, 0x7FD },
    } },
    .{ .tag = "none", .top_to_bottom = false, .standard_chars = "", .ranges = &.{}, .nonbase = &.{} },
    .{ .tag = "olck", .top_to_bottom = false, .standard_chars = "\xE1\xB1\x9B", .ranges = &.{
        .{ 0x1C50, 0x1C7F },
    }, .nonbase = &.{} },
    .{ .tag = "orkh", .top_to_bottom = false, .standard_chars = "\xF0\x90\xB0\x97", .ranges = &.{
        .{ 0x10C00, 0x10C4F },
    }, .nonbase = &.{} },
    .{ .tag = "osge", .top_to_bottom = false, .standard_chars = "\xF0\x90\x93\x82 \xF0\x90\x93\xAA", .ranges = &.{
        .{ 0x104B0, 0x104FF },
    }, .nonbase = &.{} },
    .{ .tag = "osma", .top_to_bottom = false, .standard_chars = "\xF0\x90\x92\x86 \xF0\x90\x92\xA0", .ranges = &.{
        .{ 0x10480, 0x104AF },
    }, .nonbase = &.{} },
    .{ .tag = "rohg", .top_to_bottom = false, .standard_chars = "\xF0\x90\xB4\xB0", .ranges = &.{
        .{ 0x10D00, 0x10D3F },
    }, .nonbase = &.{} },
    .{ .tag = "saur", .top_to_bottom = false, .standard_chars = "\xEA\xA2\x9D \xEA\xA3\x90", .ranges = &.{
        .{ 0xA880, 0xA8DF },
    }, .nonbase = &.{
        .{ 0xA880, 0xA881 }, .{ 0xA8B4, 0xA8C5 },
    } },
    .{ .tag = "shaw", .top_to_bottom = false, .standard_chars = "\xF0\x90\x91\xB4", .ranges = &.{
        .{ 0x10450, 0x1047F },
    }, .nonbase = &.{} },
    .{ .tag = "sinh", .top_to_bottom = false, .standard_chars = "\xE0\xB6\xA7", .ranges = &.{
        .{ 0xD80, 0xDFF },
    }, .nonbase = &.{
        .{ 0xD81, 0xD81 }, .{ 0xDCA, 0xDCA }, .{ 0xDD2, 0xDD6 },
    } },
    .{ .tag = "sund", .top_to_bottom = false, .standard_chars = "\xE1\xAE\xB0", .ranges = &.{
        .{ 0x1B80, 0x1BBF }, .{ 0x1CC0, 0x1CCF },
    }, .nonbase = &.{
        .{ 0x1B80, 0x1B82 }, .{ 0x1BA1, 0x1BAD },
    } },
    .{ .tag = "taml", .top_to_bottom = false, .standard_chars = "\xE0\xAF\xA6", .ranges = &.{
        .{ 0xB80, 0xBFF }, .{ 0x11FC0, 0x11FFF },
    }, .nonbase = &.{
        .{ 0xB82, 0xB82 }, .{ 0xBC0, 0xBC2 }, .{ 0xBCD, 0xBCD },
    } },
    .{ .tag = "tavt", .top_to_bottom = false, .standard_chars = "\xEA\xAA\x92 \xEA\xAA\xAB", .ranges = &.{
        .{ 0xAA80, 0xAADF },
    }, .nonbase = &.{
        .{ 0xAAB0, 0xAAB0 }, .{ 0xAAB2, 0xAAB4 }, .{ 0xAAB7, 0xAAB8 }, .{ 0xAABE, 0xAABF }, .{ 0xAAC1, 0xAAC1 },
    } },
    .{ .tag = "telu", .top_to_bottom = false, .standard_chars = "\xE0\xB1\xA6 \xE0\xB1\xA7", .ranges = &.{
        .{ 0xC00, 0xC7F },
    }, .nonbase = &.{
        .{ 0xC00, 0xC00 }, .{ 0xC04, 0xC04 }, .{ 0xC3C, 0xC3C }, .{ 0xC3E, 0xC40 }, .{ 0xC46, 0xC56 }, .{ 0xC62, 0xC63 },
    } },
    .{ .tag = "tfng", .top_to_bottom = false, .standard_chars = "\xE2\xB5\x94", .ranges = &.{
        .{ 0x2D30, 0x2D7F },
    }, .nonbase = &.{} },
    .{ .tag = "thai", .top_to_bottom = false, .standard_chars = "\xE0\xB8\xB2 \xE0\xB9\x85 \xE0\xB9\x90", .ranges = &.{
        .{ 0xE00, 0xE7F },
    }, .nonbase = &.{
        .{ 0xE31, 0xE31 }, .{ 0xE34, 0xE3A }, .{ 0xE47, 0xE4E },
    } },
    .{ .tag = "vaii", .top_to_bottom = false, .standard_chars = "\xEA\x98\x93 \xEA\x96\x9C \xEA\x96\xB4", .ranges = &.{
        .{ 0xA500, 0xA63F },
    }, .nonbase = &.{} },
    .{ .tag = "limb", .top_to_bottom = false, .standard_chars = "o", .ranges = &.{
        .{ 0x1900, 0x194F },
    }, .nonbase = &.{
        .{ 0x1920, 0x1922 }, .{ 0x1927, 0x1934 }, .{ 0x1937, 0x193B },
    } },
    .{ .tag = "orya", .top_to_bottom = false, .standard_chars = "o", .ranges = &.{
        .{ 0xB00, 0xB7F },
    }, .nonbase = &.{
        .{ 0xB01, 0xB02 }, .{ 0xB3C, 0xB3C }, .{ 0xB3F, 0xB3F }, .{ 0xB41, 0xB44 }, .{ 0xB4D, 0xB56 }, .{ 0xB62, 0xB63 },
    } },
    .{ .tag = "sylo", .top_to_bottom = false, .standard_chars = "o", .ranges = &.{
        .{ 0xA800, 0xA82F },
    }, .nonbase = &.{
        .{ 0xA802, 0xA802 }, .{ 0xA806, 0xA806 }, .{ 0xA80B, 0xA80B }, .{ 0xA825, 0xA826 }, .{ 0xA82C, 0xA82C },
    } },
    .{ .tag = "tibt", .top_to_bottom = false, .standard_chars = "o", .ranges = &.{
        .{ 0xF00, 0xFFF },
    }, .nonbase = &.{
        .{ 0xF18, 0xF19 }, .{ 0xF35, 0xF35 }, .{ 0xF37, 0xF37 }, .{ 0xF39, 0xF39 }, .{ 0xF3E, 0xF3F }, .{ 0xF71, 0xF7E },
        .{ 0xF80, 0xF84 }, .{ 0xF86, 0xF87 }, .{ 0xF8D, 0xFBC },
    } },
    .{ .tag = "hani", .top_to_bottom = false, .standard_chars = "\xE7\x94\xB0 \xE5\x9B\x97", .ranges = &.{
        .{ 0x1100, 0x11FF },   .{ 0x2E80, 0x2EFF },   .{ 0x2F00, 0x2FDF },   .{ 0x2FF0, 0x2FFF },   .{ 0x3000, 0x303F },   .{ 0x3040, 0x309F },
        .{ 0x30A0, 0x30FF },   .{ 0x3100, 0x312F },   .{ 0x3130, 0x318F },   .{ 0x3190, 0x319F },   .{ 0x31A0, 0x31BF },   .{ 0x31C0, 0x31EF },
        .{ 0x31F0, 0x31FF },   .{ 0x3300, 0x33FF },   .{ 0x3400, 0x4DBF },   .{ 0x4DC0, 0x4DFF },   .{ 0x4E00, 0x9FFF },   .{ 0xA960, 0xA97F },
        .{ 0xAC00, 0xD7AF },   .{ 0xD7B0, 0xD7FF },   .{ 0xF900, 0xFAFF },   .{ 0xFE10, 0xFE1F },   .{ 0xFE30, 0xFE4F },   .{ 0xFF00, 0xFFEF },
        .{ 0x1AFF0, 0x1AFFF }, .{ 0x1B000, 0x1B0FF }, .{ 0x1B100, 0x1B12F }, .{ 0x1B130, 0x1B16F }, .{ 0x1D300, 0x1D35F }, .{ 0x20000, 0x2A6DF },
        .{ 0x2A700, 0x2B73F }, .{ 0x2B740, 0x2B81F }, .{ 0x2B820, 0x2CEAF }, .{ 0x2CEB0, 0x2EBEF }, .{ 0x2EBF0, 0x2EE5D }, .{ 0x2F800, 0x2FA1F },
        .{ 0x30000, 0x3134A }, .{ 0x31350, 0x323AF }, .{ 0x323B0, 0x33479 },
    }, .nonbase = &.{
        .{ 0x302A, 0x302F }, .{ 0x3190, 0x319F },
    } },
};
