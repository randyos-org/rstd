//! Generic SpreadsheetML helpers: cell references, XML escaping, the shared
//! string table, and a worksheet `sheetData` writer.
const std = @import("std");
const zigstd = @import("../zigstd.zig");

const Allocator = zigstd.Allocator;
const Io = zigstd.StdIo;

pub const ns_main = "http://schemas.openxmlformats.org/spreadsheetml/2006/main";
pub const ns_rel = "http://schemas.openxmlformats.org/officeDocument/2006/relationships";
pub const ns_pkg_rel = "http://schemas.openxmlformats.org/package/2006/relationships";
pub const xml_decl = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n";

/// Column letters for a 0-based column index (0 -> "A", 26 -> "AA").
pub fn colName(buf: *[3]u8, col: u32) []const u8 {
    zigstd.assert(col < 16384);
    var n = col + 1;
    var i: usize = buf.len;
    while (n > 0) {
        n -= 1;
        i -= 1;
        buf[i] = 'A' + @as(u8, @intCast(n % 26));
        n /= 26;
    }
    return buf[i..];
}

/// A1-style reference for a 0-based column and 1-based row.
pub const Ref = struct {
    col: u32,
    row: u32,
    abs_col: bool = false,
    abs_row: bool = false,

    /// The same cell without `$` markers.
    pub fn relative(self: Ref) Ref {
        return .{ .col = self.col, .row = self.row };
    }

    pub fn format(self: Ref, w: *Io.Writer) Io.Writer.Error!void {
        var buf: [3]u8 = undefined;
        if (self.abs_col) try w.writeByte('$');
        try w.writeAll(colName(&buf, self.col));
        if (self.abs_row) try w.writeByte('$');
        try w.print("{d}", .{self.row});
    }
};

/// Range between two references, e.g. `A1:C10`.
pub const Range = struct {
    first: Ref,
    last: Ref,

    pub fn format(self: Range, w: *Io.Writer) Io.Writer.Error!void {
        try w.print("{f}:{f}", .{ self.first, self.last });
    }
};

/// Writes `s` with the five XML special characters escaped.
pub fn writeEscaped(w: *Io.Writer, s: []const u8) Io.Writer.Error!void {
    var start: usize = 0;
    for (s, 0..) |c, i| {
        const rep: []const u8 = switch (c) {
            '&' => "&amp;",
            '<' => "&lt;",
            '>' => "&gt;",
            '"' => "&quot;",
            '\'' => "&apos;",
            else => continue,
        };
        try w.writeAll(s[start..i]);
        try w.writeAll(rep);
        start = i + 1;
    }
    try w.writeAll(s[start..]);
}

/// `{f}` formatter that XML-escapes a string.
pub fn esc(s: []const u8) std.fmt.Alt([]const u8, formatEscaped) {
    return .{ .data = s };
}

fn formatEscaped(s: []const u8, w: *Io.Writer) Io.Writer.Error!void {
    return writeEscaped(w, s);
}

/// Workbook-wide shared string table. Keys are borrowed; callers keep them
/// alive until the table is written.
pub const SharedStrings = struct {
    map: std.array_hash_map.String(void) = .empty,
    /// Total references (the `count` attribute), as opposed to unique strings.
    refs: u32 = 0,

    pub fn index(self: *SharedStrings, gpa: Allocator, s: []const u8) !u32 {
        const gop = try self.map.getOrPut(gpa, s);
        self.refs += 1;
        return @intCast(gop.index);
    }

    pub fn write(self: *const SharedStrings, w: *Io.Writer) Io.Writer.Error!void {
        try w.print(xml_decl ++ "<sst xmlns=\"" ++ ns_main ++ "\" count=\"{d}\" uniqueCount=\"{d}\">", .{
            self.refs, self.map.count(),
        });
        for (self.map.keys()) |s| {
            // Preserve leading/trailing spaces exactly.
            try w.print("<si><t xml:space=\"preserve\">{f}</t></si>", .{esc(s)});
        }
        try w.writeAll("</sst>");
    }
};

/// Streams `<sheetData>` rows and cells. Rows and cells must be emitted in
/// ascending order, as the file format requires.
pub const SheetData = struct {
    w: *Io.Writer,
    gpa: Allocator,
    sst: *SharedStrings,
    row: u32 = 0,

    pub fn begin(self: *SheetData) !void {
        try self.w.writeAll("<sheetData>");
    }

    pub fn end(self: *SheetData) !void {
        if (self.row != 0) try self.w.writeAll("</row>");
        try self.w.writeAll("</sheetData>");
    }

    pub fn startRow(self: *SheetData, row: u32) !void {
        zigstd.assert(row > self.row);
        if (self.row != 0) try self.w.writeAll("</row>");
        try self.w.print("<row r=\"{d}\">", .{row});
        self.row = row;
    }

    fn open(self: *SheetData, col: u32, style: u32, t: ?[]const u8) !void {
        try self.w.print("<c r=\"{f}\"", .{Ref{ .col = col, .row = self.row }});
        if (style != 0) try self.w.print(" s=\"{d}\"", .{style});
        if (t) |ty| try self.w.print(" t=\"{s}\"", .{ty});
        try self.w.writeByte('>');
    }

    pub fn string(self: *SheetData, col: u32, style: u32, s: []const u8) !void {
        const idx = try self.sst.index(self.gpa, s);
        try self.open(col, style, "s");
        try self.w.print("<v>{d}</v></c>", .{idx});
    }

    pub fn number(self: *SheetData, col: u32, style: u32, v: f64) !void {
        try self.open(col, style, null);
        try self.w.print("<v>{d}</v></c>", .{v});
    }

    pub fn boolean(self: *SheetData, col: u32, style: u32, v: bool) !void {
        try self.open(col, style, "b");
        try self.w.print("<v>{d}</v></c>", .{@intFromBool(v)});
    }

    /// Formula cell with no cached value; the workbook recalculates on load.
    /// `f` is the formula without its leading `=`.
    pub fn formula(self: *SheetData, col: u32, style: u32, f: []const u8) !void {
        try self.open(col, style, null);
        try self.w.print("<f>{f}</f></c>", .{esc(f)});
    }

    /// Styled empty cell.
    pub fn blank(self: *SheetData, col: u32, style: u32) !void {
        try self.w.print("<c r=\"{f}\" s=\"{d}\"/>", .{ Ref{ .col = col, .row = self.row }, style });
    }
};

/// Excel serial date (1900 date system) for a proleptic Gregorian date.
pub fn serialDate(year: i32, month: u32, day: u32) i32 {
    // days_from_civil (H. Hinnant), relative to 1970-01-01.
    const y = if (month <= 2) year - 1 else year;
    const era = @divFloor(y, 400);
    const yoe: i32 = y - era * 400;
    const m: i32 = @intCast(month);
    const mp = @mod(m + 9, 12);
    const doy = @divFloor(153 * mp + 2, 5) + @as(i32, @intCast(day)) - 1;
    const doe = yoe * 365 + @divFloor(yoe, 4) - @divFloor(yoe, 100) + doy;
    const unix_days = era * 146097 + doe - 719468;
    return unix_days + 25569;
}

test colName {
    var buf: [3]u8 = undefined;
    try zigstd.testExpectEqualStrings("A", colName(&buf, 0));
    try zigstd.testExpectEqualStrings("Z", colName(&buf, 25));
    try zigstd.testExpectEqualStrings("AA", colName(&buf, 26));
    try zigstd.testExpectEqualStrings("EF", colName(&buf, 135));
    try zigstd.testExpectEqualStrings("XFD", colName(&buf, 16383));
}

test serialDate {
    try zigstd.testExpectEqual(@as(i32, 25569), serialDate(1970, 1, 1));
    try zigstd.testExpectEqual(@as(i32, 45658), serialDate(2025, 1, 1));
    try zigstd.testExpectEqual(@as(i32, 46307), serialDate(2026, 10, 12));
}

test writeEscaped {
    var aw: Io.Writer.Allocating = .init(zigstd.testAllocator);
    defer aw.deinit();
    try writeEscaped(&aw.writer, "a<b & \"c\">'d'");
    try zigstd.testExpectEqualStrings("a&lt;b &amp; &quot;c&quot;&gt;&apos;d&apos;", aw.written());
}
