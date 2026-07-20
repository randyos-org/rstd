//! ANSI escape code utilities (ECMA-48)
//!
//! see: https://stackoverflow.com/a/33206814/13230486
//! see: https://en.wikipedia.org/wiki/ANSI_escape_code#CSIsection
//! see: https://en.wikipedia.org/wiki/ANSI.SYS

const std = @import("std");
const log = std.log.scoped(.ansi);

// ASCII control characters
// aka `C0` control set

/// Null (\0)
pub const NUL = "\x00";
/// Start of header
pub const SOH = "\x01";
/// Start of text
pub const STX = "\x02";
/// End of text
pub const ETX = "\x03";
/// End of transmission
pub const EOT = "\x04";
/// Enquiry
pub const ENQ = "\x05";
/// Acknowledge
pub const ACK = "\x06";
/// Bell (\a)
pub const BEL = "\x07";
/// Backspace (\b)
pub const BS = "\x08";
/// Horizontal tab (\t)
pub const HT = "\x09";
/// Line feed (\n)
pub const LF = "\x0a";
/// Vertical tab (\v)
pub const VT = "\x0b";
/// Form feed (\f)
pub const FF = "\x0c";
/// Carriage return (\r)
pub const CR = "\x0d";
/// Shift out
pub const SO = "\x0e";
/// Shift in
pub const SI = "\x0f";
/// Data link escape
pub const DLE = "\x10";
/// Device control 1
pub const DC1 = "\x11";
/// Device control 2
pub const DC2 = "\x12";
/// Device control 3
pub const DC3 = "\x13";
/// Device control 4
pub const DC4 = "\x14";
/// Negative acknowledge
pub const NAK = "\x15";
/// Synchronous idle
pub const SYN = "\x16";
/// End of transmission block
pub const ETB = "\x17";
/// Cancel
pub const CAN = "\x18";
/// End of medium
pub const EM = "\x19";
/// Substitute
pub const SUB = "\x1a";
/// Escape (\e)
pub const ESC = "\x1b";
/// File separator
pub const FS = "\x1c";
/// Group separator
pub const GS = "\x1d";
/// Record separator
pub const RS = "\x1e";
/// Unit separator
pub const US = "\x1f";
/// Space
pub const SP = "\x20";
/// Delete
pub const DEL = "\x7f";

/// Carriage return + line feed combination
pub const CRLF = "\r\n";

// `C1` aka `Fe` control set
/// Control Sequence Introducer (CSI) character used in ANSI escape codes
pub const CSI = ESC ++ "[";
/// Operating System Command (OSC) character used in ANSI escape codes
pub const OSC = ESC ++ "]";

/// CSI suffixes
pub const CsiSuffix = enum {
    /// Cursor Up (A)
    CUU,
    /// Cursor Down (B)
    CUD,
    /// Cursor Forward (C)
    CUF,
    /// Cursor Back (D)
    CUB,
    /// Cursor Next Line (E)
    CNL,
    /// Cursor Previous Line (F)
    CPL,
    /// Cursor Horizontal Absolute (G)
    CHA,
    /// Cursor Position (H)
    CUP,
    /// Erase in Display (J)
    ED,
    /// Erase in Line (K)
    EL,
    /// Scroll Up (S)
    SU,
    /// Scroll Down (T)
    SD,
    /// Select Graphic Rendition (m)
    SGR,

    /// Horizontal and Vertical Position (f)
    /// This one is less reliable than CUP
    HVP,

    // no arguments
    /// Save Cursor Position (s)
    SCP,
    /// Restore Cursor Position (u)
    RCP,
    /// Device Status Report (6n)
    DSR,

    // response codes
    /// Cursor Position Report (response to DSR) (R)
    CPR,

    /// The literal byte(s) terminating a CSI sequence for this command
    pub fn suffix(self: CsiSuffix) []const u8 {
        return switch (self) {
            .CUU => "A",
            .CUD => "B",
            .CUF => "C",
            .CUB => "D",
            .CNL => "E",
            .CPL => "F",
            .CHA => "G",
            .CUP => "H",
            .ED => "J",
            .EL => "K",
            .SU => "S",
            .SD => "T",
            .SGR => "m",
            .HVP => "f",
            .SCP => "s",
            .RCP => "u",
            .DSR => "6n",
            .CPR => "R",
        };
    }
};

/// Select Graphic Rendition (SGR) CSI suffix
pub const SGR = CsiSuffix.suffix(.SGR);
/// Clear screen CSI suffix
pub const CLS = "2" ++ CsiSuffix.suffix(.ED);
/// Cursor to home (top left) CSI suffix
pub const HOME = CsiSuffix.suffix(.CUP);

pub const CsiCodes = struct {
    /// Cursor Up (A)
    pub const CUU = "A";
    /// Cursor Down (B)
    pub const CUD = "B";
    /// Cursor Forward (C)
    pub const CUF = "C";
    /// Cursor Back (D)
    pub const CUB = "D";
    /// Cursor Next Line (E)
    pub const CNL = "E";
    /// Cursor Previous Line (F)
    pub const CPL = "F";
    /// Cursor Horizontal Absolute (G)
    pub const CHA = "G";
    /// Cursor Position (H)
    pub const CUP = "H";
    /// Erase in Display (J)
    pub const ED = "J";
    /// Erase in Line (K)
    pub const EL = "K";
    /// Scroll Up (S)
    pub const SU = "S";
    /// Scroll Down (T)
    pub const SD = "T";
    /// Select Graphic Rendition (m)
    pub const SGR = "m";
    /// Horizontal and Vertical Position (f)
    pub const HVP = "f";
    /// Save Cursor Position (s)
    pub const SCP = "s";
    /// Restore Cursor Position (u)
    pub const RCP = "u";
    /// Device Status Report (6n)
    pub const DSR = "6n";
    /// Cursor Position Report (response to DSR) (R)
    pub const CPR = "R";
};

const CursorPos = struct {
    n: u32 = 1,
};

const CursorPos2D = struct {
    row: u32 = 1,
    col: u32 = 1,
};

const AbsCursorPos2D = struct {
    row: u32 = 0,
    col: u32 = 0,
};

const CursorEraseMode = enum(u32) {
    ToEnd = 0,
    ToStart = 1,
    Whole = 2,
    WholeScrollback = 3,
};

/// Writes a CSI escape sequence (e.g. "\x1b[3;4H") into `buf`, joining
/// `args` with ';' and terminating with `cmd`'s suffix. Returns the
/// written slice of `buf`.
fn csiSeq(buf: []u8, cmd: CsiSuffix, args: []const u32) ![]const u8 {
    @memcpy(buf[0..CSI.len], CSI);
    var len: usize = CSI.len;
    for (args, 0..) |arg, i| {
        if (i != 0) {
            buf[len] = ';';
            len += 1;
        }
        len += (try std.fmt.bufPrint(buf[len..], "{d}", .{arg})).len;
    }
    const suf = cmd.suffix();
    @memcpy(buf[len..][0..suf.len], suf);
    len += suf.len;
    return buf[0..len];
}

/// Cursor Up `n` lines
pub fn cuu(buf: []u8, n: CursorPos) ![]const u8 {
    return csiSeq(buf, .CUU, &.{n.n});
}
pub const cursorUp = cuu;

/// Cursor Down `n` lines
pub fn cud(buf: []u8, n: CursorPos) ![]const u8 {
    return csiSeq(buf, .CUD, &.{n.n});
}
pub const cursorDown = cud;

/// Cursor Forward `n` columns
pub fn cuf(buf: []u8, n: CursorPos) ![]const u8 {
    return csiSeq(buf, .CUF, &.{n.n});
}
pub const cursorForward = cuf;

/// Cursor Back `n` columns
pub fn cub(buf: []u8, n: CursorPos) ![]const u8 {
    return csiSeq(buf, .CUB, &.{n.n});
}
pub const cursorBack = cub;

/// Cursor Next Line, `n` lines down, column 1
pub fn cnl(buf: []u8, n: CursorPos) ![]const u8 {
    return csiSeq(buf, .CNL, &.{n.n});
}
pub const cursorNextLine = cnl;

/// Cursor Previous Line, `n` lines up, column 1
pub fn cpl(buf: []u8, n: CursorPos) ![]const u8 {
    return csiSeq(buf, .CPL, &.{n.n});
}
pub const cursorPreviousLine = cpl;

/// Cursor Horizontal Absolute: move to column `col` on the current line
pub fn cha(buf: []u8, col: CursorPos) ![]const u8 {
    return csiSeq(buf, .CHA, &.{col.n});
}
pub const cursorHorizAbsolute = cha;

/// Cursor Position: move to (`row`, `col`), both 1-based
pub fn cup(buf: []u8, rc: AbsCursorPos2D) ![]const u8 {
    return csiSeq(buf, .CUP, &.{ rc.row, rc.col });
}
pub const cursorPosition = cup;

/// Erase in Display (0 = cursor to end, 1 = start to cursor, 2 = whole screen)
pub fn ed(buf: []u8, mode: CursorEraseMode) ![]const u8 {
    return csiSeq(buf, .ED, &.{@intFromEnum(mode)});
}
pub const eraseDisplay = ed;
pub fn clearScreen(buf: []u8) ![]const u8 {
    return ed(buf, .Whole);
}

/// Erase in Line (0 = cursor to end, 1 = start to cursor, 2 = whole line)
pub fn el(buf: []u8, mode: CursorEraseMode) ![]const u8 {
    return csiSeq(buf, .EL, &.{@intFromEnum(mode)});
}
pub const eraseLine = el;

/// Scroll Up `n` lines
pub fn su(buf: []u8, n: CursorPos) ![]const u8 {
    return csiSeq(buf, .SU, &.{n.n});
}
pub const scrollUp = su;

/// Scroll Down `n` lines
pub fn sd(buf: []u8, n: CursorPos) ![]const u8 {
    return csiSeq(buf, .SD, &.{n.n});
}
pub const scrollDown = sd;

/// Select Graphic Rendition: apply the given (possibly empty) list of raw
/// numeric codes, e.g. `sgr(buf, &.{ 1, 31 })` for bold red. Prefer
/// `sgrStyled` for a typed, chainable version that also covers the
/// variable-argument true-color/indexed-color codes.
pub fn sgr(buf: []u8, codes: []const u32) ![]const u8 {
    return csiSeq(buf, .SGR, codes);
}
pub const textFormatCode = sgr;

/// A single Select Graphic Rendition style. Most styles are fixed codes
/// with no arguments; the true-color and 256-color-palette variants carry
/// the (runtime) color value they encode.
pub const SgrStyle = union(enum) {
    reset,
    bold,
    dim,
    italic,
    underline,
    blink,
    invert,
    hidden,
    strike,

    /// Some terminals treat this as "not bold", but ECMA-48 actually defines
    /// 21 as double underline -- not reliable for turning off bold; prefer `no_dim`.
    no_bold,
    /// The spec's real intensity reset -- turns off both bold and dim, since
    /// they're the same "intensity" attribute (21/`no_bold` is not reliable for this).
    no_dim,
    no_italic,
    no_underline,
    no_blink,
    no_invert,
    no_hidden,
    no_strike,

    fg_black,
    fg_dark_red,
    fg_dark_green,
    fg_dark_yellow,
    fg_dark_blue,
    fg_dark_magenta,
    fg_dark_cyan,
    fg_light_gray,
    /// Reset the foreground color to the terminal's default
    fg_default,
    fg_dark_gray,
    fg_red,
    fg_green,
    fg_yellow,
    fg_blue,
    fg_magenta,
    fg_cyan,
    fg_white,
    /// 256-color palette index
    fg_indexed: u8,
    /// 24-bit true color
    fg_rgb: RGB,

    bg_black,
    bg_dark_red,
    bg_dark_green,
    bg_dark_yellow,
    bg_dark_blue,
    bg_dark_magenta,
    bg_dark_cyan,
    bg_light_gray,
    /// Reset the background color to the terminal's default
    bg_default,
    bg_dark_gray,
    bg_red,
    bg_green,
    bg_yellow,
    bg_blue,
    bg_magenta,
    bg_cyan,
    bg_white,
    /// 256-color palette index
    bg_indexed: u8,
    /// 24-bit true color
    bg_rgb: RGB,

    pub const RGB = struct { r: u8, g: u8, b: u8 };

    pub fn fgRgb(r: u8, g: u8, b: u8) SgrStyle {
        return .{ .fg_rgb = .{ .r = r, .g = g, .b = b } };
    }

    pub fn bgRgb(r: u8, g: u8, b: u8) SgrStyle {
        return .{ .bg_rgb = .{ .r = r, .g = g, .b = b } };
    }

    /// This style's numeric SGR parameter(s) (e.g. "1" or "38;2;10;20;30"),
    /// with no leading/trailing ';' -- `sgrStyled` joins multiple styles.
    /// `buf` only needs to be large enough for the true-color/indexed-color
    /// cases; the fixed styles don't touch it.
    fn params(self: SgrStyle, buf: []u8) ![]const u8 {
        return switch (self) {
            .reset => "0",
            .bold => "1",
            .dim => "2",
            .italic => "3",
            .underline => "4",
            .blink => "5",
            .invert => "7",
            .hidden => "8",
            .strike => "9",

            .no_bold => "21",
            .no_dim => "22",
            .no_italic => "23",
            .no_underline => "24",
            .no_blink => "25",
            .no_invert => "27",
            .no_hidden => "28",
            .no_strike => "29",

            .fg_black => "30",
            .fg_dark_red => "31",
            .fg_dark_green => "32",
            .fg_dark_yellow => "33",
            .fg_dark_blue => "34",
            .fg_dark_magenta => "35",
            .fg_dark_cyan => "36",
            .fg_light_gray => "37",

            .fg_default => "39",

            .fg_dark_gray => "90",
            .fg_red => "91",
            .fg_green => "92",
            .fg_yellow => "93",
            .fg_blue => "94",
            .fg_magenta => "95",
            .fg_cyan => "96",
            .fg_white => "97",

            .fg_indexed => |n| try std.fmt.bufPrint(buf, "38;5;{d}", .{n}),
            .fg_rgb => |rgb| try std.fmt.bufPrint(buf, "38;2;{d};{d};{d}", .{ rgb.r, rgb.g, rgb.b }),

            .bg_black => "40",
            .bg_dark_red => "41",
            .bg_dark_green => "42",
            .bg_dark_yellow => "43",
            .bg_dark_blue => "44",
            .bg_dark_magenta => "45",
            .bg_dark_cyan => "46",
            .bg_light_gray => "47",

            .bg_default => "49",

            .bg_dark_gray => "100",
            .bg_red => "101",
            .bg_green => "102",
            .bg_yellow => "103",
            .bg_blue => "104",
            .bg_magenta => "105",
            .bg_cyan => "106",
            .bg_white => "107",

            .bg_indexed => |n| try std.fmt.bufPrint(buf, "48;5;{d}", .{n}),
            .bg_rgb => |rgb| try std.fmt.bufPrint(buf, "48;2;{d};{d};{d}", .{ rgb.r, rgb.g, rgb.b }),
        };
    }
};

/// Select Graphic Rendition, from a chain of typed styles, e.g.
/// `sgrStyled(buf, &.{ .bold, .bg_blue, ansi.SgrStyle.fgRgb(100, 150, 200) })`.
pub fn sgrStyled(buf: []u8, styles: []const SgrStyle) ![]const u8 {
    @memcpy(buf[0..CSI.len], CSI);
    var len: usize = CSI.len;
    // Long enough for the widest case, "48;2;255;255;255".
    var param_buf: [16]u8 = undefined;
    for (styles, 0..) |style, i| {
        if (i != 0) {
            buf[len] = ';';
            len += 1;
        }
        const p = try style.params(&param_buf);
        @memcpy(buf[len..][0..p.len], p);
        len += p.len;
    }
    buf[len] = CsiSuffix.suffix(.SGR);
    len += 1;
    return buf[0..len];
}
pub const styleText = sgrStyled;

/// Horizontal and Vertical Position: move to (`row`, `col`), both 1-based.
/// Equivalent to `cup`, but less widely supported.
pub fn hvp(buf: []u8, rc: AbsCursorPos2D) ![]const u8 {
    return csiSeq(buf, .HVP, &.{ rc.row, rc.col });
}
pub const moveToHorizVertPosition = hvp;

/// Save Cursor Position (no arguments, so no buffer is needed)
pub const scp = CSI ++ "s";
/// Restore Cursor Position (no arguments, so no buffer is needed)
pub const rcp = CSI ++ "u";
/// Device Status Report (no arguments, so no buffer is needed). The
/// terminal replies with a CPR sequence, which this module doesn't
/// generate since it's a response, not something we send.
pub const dsr = CSI ++ "6n";

test "csiSeq" {
    var buf: [32]u8 = undefined;
    try std.testing.expectEqualStrings("\x1b[3;4H", try cup(&buf, .{ .row = 3, .col = 4 }));
    try std.testing.expectEqualStrings("\x1b[2J", try ed(&buf, .Whole));
    try std.testing.expectEqualStrings("\x1b[1;31m", try sgr(&buf, &.{ 1, 31 }));
    try std.testing.expectEqualStrings("\x1b[m", try sgr(&buf, &.{}));
    try std.testing.expectEqualStrings("\x1b[6n", dsr);
}

test "sgrStyled" {
    var buf: [32]u8 = undefined;
    try std.testing.expectEqualStrings(
        "\x1b[1;44;38;2;100;150;200m",
        try sgrStyled(&buf, &.{ .bold, .bg_blue, SgrStyle.fgRgb(100, 150, 200) }),
    );
    try std.testing.expectEqualStrings("\x1b[38;5;201m", try sgrStyled(&buf, &.{.{ .fg_indexed = 201 }}));
    try std.testing.expectEqualStrings("\x1b[m", try sgrStyled(&buf, &.{}));
}

pub const SgrCode = struct {
    pub const reset = "0";
    pub const bold = "1";
    pub const dim = "2";
    pub const italic = "3";
    pub const underline = "4";
    pub const blink = "5";
    pub const invert = "7";
    pub const hidden = "8";
    pub const strike = "9";

    pub const no_bold = "21";
    pub const no_dim = "22";
    pub const no_italic = "23";
    pub const no_underline = "24";
    pub const no_blink = "25";
    pub const no_invert = "27";
    pub const no_hidden = "28";
    pub const no_strike = "29";

    pub const fg_black = "30";
    pub const fg_dark_red = "31";
    pub const fg_dark_green = "32";
    pub const fg_dark_yellow = "33";
    pub const fg_dark_blue = "34";
    pub const fg_dark_magenta = "35";
    pub const fg_dark_cyan = "36";
    pub const fg_light_gray = "37";

    pub const fg_default = "39";

    pub const fg_dark_gray = "90";
    pub const fg_red = "91";
    pub const fg_green = "92";
    pub const fg_yellow = "93";
    pub const fg_blue = "94";
    pub const fg_magenta = "95";
    pub const fg_cyan = "96";
    pub const fg_white = "97";

    pub const fg_indexed = "38;5";
    pub const fg_rgb = "38;2";

    pub const bg_black = "40";
    pub const bg_dark_red = "41";
    pub const bg_dark_green = "42";
    pub const bg_dark_yellow = "43";
    pub const bg_dark_blue = "44";
    pub const bg_dark_magenta = "45";
    pub const bg_dark_cyan = "46";
    pub const bg_light_gray = "47";

    pub const bg_default = "49";

    pub const bg_dark_gray = "100";
    pub const bg_red = "101";
    pub const bg_green = "102";
    pub const bg_yellow = "103";
    pub const bg_blue = "104";
    pub const bg_magenta = "105";
    pub const bg_cyan = "106";
    pub const bg_white = "107";

    pub const bg_indexed = "48;5";
    pub const bg_rgb = "48;2";
};
