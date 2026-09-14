const bytes = @import("../../../../bytes.zig");
const zigstd = @import("../../../../zigstd.zig");

/// A DICT: operand(s)-then-operator tokens (Top/Font/Private DICTs all use
/// this encoding). table_data is exactly this DICT's own byte range.
const Dict = @This();

table_data: []const u8,

/// DICT operators this module reads. Escaped ("12 N") operators are encoded
/// as 0x100 | N. Non-exhaustive - the spec defines many more we don't need.
pub const Operator = enum(u16) {
    charset = 15,
    encoding = 16,
    char_strings = 17,
    private = 18,
    /// Only meaningful within a Private DICT, not a Top/Font DICT.
    subrs = 19,
    charstring_type = 0x100 | 6,
    /// Presence (regardless of operand values) marks a CID-keyed font.
    ros = 0x100 | 30,
    fd_array = 0x100 | 36,
    fd_select = 0x100 | 37,
    _,
};

/// The raw, un-decoded operand bytes immediately preceding `op`'s
/// occurrence, or null if this DICT has no such operator. Per spec an
/// operator should appear at most once; this returns the first match.
pub fn getOperandBytes(self: Dict, op: Operator) ?[]const u8 {
    var pos: usize = 0;
    while (pos < self.table_data.len) {
        const start = pos;
        while (pos < self.table_data.len and self.table_data[pos] >= 28) {
            skipDictOperand(self.table_data, &pos);
        }
        const operand_end = pos;
        if (pos >= self.table_data.len) break;

        var code: u16 = self.table_data[pos];
        pos += 1;
        if (code == 12) {
            code = 0x100 | @as(u16, self.table_data[pos]);
            pos += 1;
        }

        if (code == @intFromEnum(op)) return self.table_data[start..operand_end];
    }
    return null;
}

pub fn has(self: Dict, op: Operator) bool {
    return self.getOperandBytes(op) != null;
}

/// Decodes up to `out.len` leading integer operands for `op`, returning
/// how many were found. Real-number operands aren't supported - no
/// operator this module reads (offsets, sizes, counts, the CID marker)
/// ever uses one.
pub fn getInts(self: Dict, op: Operator, out: []i32) usize {
    const operands = self.getOperandBytes(op) orelse return 0;
    var pos: usize = 0;
    var n: usize = 0;
    while (n < out.len and pos < operands.len) : (n += 1) {
        out[n] = decodeDictInt(operands, &pos);
    }
    return n;
}

pub fn getInt(self: Dict, op: Operator) ?i32 {
    var out: [1]i32 = undefined;
    return if (self.getInts(op, &out) == 1) out[0] else null;
}

/// Decodes one DICT integer operand starting at data[pos.*], advancing
/// pos.* past it. Caller must have already checked data[pos.*] is an
/// integer-operand lead byte (28, 29, or 32-254) - not a real-number lead
/// (30) or an operator byte (0-21, 12).
fn decodeDictInt(data: []const u8, pos: *usize) i32 {
    const b0 = data[pos.*];
    pos.* += 1;
    return switch (b0) {
        32...246 => @as(i32, b0) - 139,
        247...250 => blk: {
            const b1 = data[pos.*];
            pos.* += 1;
            break :blk (@as(i32, b0) - 247) * 256 + @as(i32, b1) + 108;
        },
        251...254 => blk: {
            const b1 = data[pos.*];
            pos.* += 1;
            break :blk -(@as(i32, b0) - 251) * 256 - @as(i32, b1) - 108;
        },
        28 => blk: {
            const v = bytes.getI16BE(data[pos.*..][0..2]);
            pos.* += 2;
            break :blk @as(i32, v);
        },
        29 => blk: {
            const v = bytes.getI32BE(data[pos.*..][0..4]);
            pos.* += 4;
            break :blk v;
        },
        else => unreachable,
    };
}

/// Skips one DICT operand (integer or real) starting at data[pos.*],
/// without decoding it - used while scanning past operators we don't care
/// about. Real numbers (lead byte 30) are nibble-encoded and terminated by
/// a 0xF nibble, so they can't be skipped by a fixed byte count like the
/// integer forms.
fn skipDictOperand(data: []const u8, pos: *usize) void {
    if (data[pos.*] == 30) {
        pos.* += 1;
        while (pos.* < data.len) {
            const v = data[pos.*];
            pos.* += 1;
            if ((v & 0x0F) == 0x0F or (v >> 4) == 0x0F) break;
        }
    } else {
        _ = decodeDictInt(data, pos);
    }
}

/// Test-support: appends a DICT integer operand always using the fixed
/// 5-byte form (opcode 29 + 4-byte big-endian value), so a DICT's total
/// byte length doesn't depend on the (often forward-referencing) operand
/// values it encodes. Reused by __root__.zig's tests.
pub fn appendDictInt(buf: []u8, at: usize, value: i32) usize {
    buf[at] = 29;
    bytes.setI32BE(buf[at + 1 ..], value);
    return at + 5;
}

test "Dict.getInt finds a one-byte and an escaped operator, skipping a real-number operand" {
    // Encodes: [real (skipped, irrelevant)] charset=15
    // then the operators we actually want: [int 42] CharStrings=17   [int -3] 12 30 (ROS)
    var buf: [64]u8 = undefined;
    var pos: usize = 0;

    // A minimal real-number operand: lead byte 30, then a single terminator nibble pair (0xFF).
    buf[pos] = 30;
    pos += 1;
    buf[pos] = 0xFF; // both nibbles 0xF: terminates immediately
    pos += 1;
    buf[pos] = 15; // charset
    pos += 1;

    pos = appendDictInt(&buf, pos, 42);
    buf[pos] = 17; // CharStrings
    pos += 1;

    pos = appendDictInt(&buf, pos, -3);
    buf[pos] = 12;
    buf[pos + 1] = 30; // ROS (escaped)
    pos += 2;

    const dict: Dict = .{ .table_data = buf[0..pos] };
    try zigstd.testExpectEqual(true, dict.has(.charset));
    try zigstd.testExpectEqual(@as(?i32, 42), dict.getInt(.char_strings));
    try zigstd.testExpectEqual(@as(?i32, -3), dict.getInt(.ros));
    try zigstd.testExpectEqual(@as(?i32, null), dict.getInt(.fd_array));
}
