//! Managing the stack and pushing onto it: FreeType 2.14.3 `ttinterp.c`
//! `Ins_DUP` to `Ins_ROLL` and `Ins_NPUSHB` to `Ins_PUSHW`.

const opcode = @import("../opcode.zig");
const Engine = @import("../Engine.zig");
const Decoder = @import("../Decoder.zig");

const Error = Engine.Error;

pub fn dup(e: *Engine) void {
    e.arg(1).* = e.arg(0).*;
}

pub fn clear(e: *Engine) void {
    e.new_top = 0;
}

pub fn swap(e: *Engine) void {
    const l = e.arg(0).*;
    e.arg(0).* = e.arg(1).*;
    e.arg(1).* = l;
}

pub fn depth(e: *Engine) void {
    e.arg(0).* = @intCast(e.top);
}

/// CINDEX: copies the element `L` below the index; 0 for a bad `L`.
pub fn cindex(e: *Engine) Error!void {
    const l = e.arg(0).*;
    if (l <= 0 or l > e.args) {
        try e.pedanticError(error.InvalidReference);
        e.arg(0).* = 0;
    } else e.arg(0).* = e.stack[e.args - @as(usize, @intCast(l))];
}

/// MINDEX: moves the element `L` below the index to the top.
pub fn mindex(e: *Engine) Error!void {
    const l = e.arg(0).*;
    if (l <= 0 or l > e.args) return e.pedanticError(error.InvalidReference);
    const ul: usize = @intCast(l);
    const from = e.args - ul;
    const k = e.stack[from];
    @memmove(e.stack[from .. e.args - 1], e.stack[from + 1 .. e.args]);
    e.stack[e.args - 1] = k;
}

pub fn roll(e: *Engine) void {
    const a = e.arg(2).*;
    const b = e.arg(1).*;
    const c = e.arg(0).*;
    e.arg(2).* = c;
    e.arg(1).* = a;
    e.arg(0).* = b;
}

/// The PUSH family: the operands after the opcode, bytes or signed words.
pub fn push(e: *Engine) Error!void {
    var d: Decoder = .{ .bytecode = e.code[@backingInt(e.cur)], .pc = e.ip };
    const ins = (try d.next()).?;
    const n = ins.operands.len();
    // NPUSHB and NPUSHW push past the table's count
    if (e.op == 0x40 or e.op == 0x41) {
        if (n > e.stack.len - e.top) return error.StackOverflow;
        e.new_top += n;
    }
    for (0..n) |i| e.arg(i).* = ins.operands.get(i);
    e.length = d.pc - e.ip;
}

test "push lengths and stack effects agree" {
    for (0..256) |op| {
        const b: u8 = @intCast(op);
        if (!opcode.isPush(b) or b == 0x40 or b == 0x41) continue;
        const n: usize = if (opcode.isPushWords(b)) b - opcode.pushw_first + 1 else b - opcode.pushb_first + 1;
        try @import("std").testing.expectEqual(n, opcode.stack_effects[b].push);
    }
}
