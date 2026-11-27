//! Reads instructions from TrueType bytecode.
//!
//! Ported from fontations `read-fonts/src/tables/glyf/bytecode/decode.rs`
//! and `instruction.rs` (commit 4f115482, MIT).

const opcode = @import("opcode.zig");

const Decoder = @This();

bytecode: []const u8,
/// The offset of the next instruction.
pc: usize = 0,

/// Operands past the end of the bytecode (`Code_Overflow`).
pub const Error = error{CodeOverflow};

/// One decoded instruction.
pub const Instruction = struct {
    op: u8,
    /// A push's operands.
    operands: Operands = .{},
    /// Its offset in the bytecode.
    pc: usize,
};

/// A push instruction's inline operands: bytes, or big-endian words.
pub const Operands = struct {
    bytes: []const u8 = &.{},
    words: bool = false,

    pub fn len(o: Operands) usize {
        return if (o.words) o.bytes.len / 2 else o.bytes.len;
    }

    /// Operand `i`: a byte zero-extended, a word sign-extended.
    pub fn get(o: Operands, i: usize) i32 {
        if (!o.words) return o.bytes[i];
        const w: u16 = (@as(u16, o.bytes[2 * i]) << 8) | o.bytes[2 * i + 1];
        return @as(i16, @bitCast(w));
    }
};

/// The next instruction, null at the end of the bytecode.
pub fn next(d: *Decoder) Error!?Instruction {
    if (d.pc >= d.bytecode.len) return null;
    const pc = d.pc;
    const op = d.bytecode[pc];
    var len: usize = undefined;
    var count_len: usize = 0;
    const l = opcode.lengths[op];
    if (l < 0) {
        if (pc + 1 >= d.bytecode.len) return error.CodeOverflow;
        len = @as(usize, @abs(l)) * d.bytecode[pc + 1] + 2;
        count_len = 1;
    } else len = @intCast(l);
    if (len > d.bytecode.len - pc) return error.CodeOverflow;
    d.pc = pc + len;
    var ins: Instruction = .{ .op = op, .pc = pc };
    if (opcode.isPush(op)) ins.operands = .{ .bytes = d.bytecode[pc + 1 + count_len .. pc + len], .words = opcode.isPushWords(op) };
    return ins;
}
