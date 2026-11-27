//! Control flow: FreeType 2.14.3 `ttinterp.c` `Ins_IF`, `Ins_ELSE`,
//! `Ins_JMPR`, `Ins_JROT`, `Ins_JROF`.

const Engine = @import("../Engine.zig");

const Error = Engine.Error;

/// IF: on false, skips to the matching ELSE or EIF.
pub fn @"if"(e: *Engine) Error!void {
    if (e.arg(0).* != 0) return;
    var depth: u32 = 1;
    while (true) {
        try e.skipCode();
        switch (e.op) {
            0x58 => depth += 1,
            0x1B => if (depth == 1) return,
            0x59 => {
                depth -= 1;
                if (depth == 0) return;
            },
            else => {},
        }
    }
}

/// ELSE, reached by running the true branch: skips to the matching EIF.
pub fn @"else"(e: *Engine) Error!void {
    var depth: u32 = 1;
    while (depth != 0) {
        try e.skipCode();
        switch (e.op) {
            0x58 => depth += 1,
            0x59 => depth -= 1,
            else => {},
        }
    }
}

/// JMPR by `args[0]` from this instruction; a zero jump on an empty stack,
/// before the program's start or past its function's ENDF is an error.
pub fn jump(e: *Engine) Error!void {
    const off = e.arg(0).*;
    if (off == 0 and e.args == 0) return error.BadArgument;
    const ip = @as(i64, @intCast(e.ip)) + off;
    if (ip < 0 or (e.call_top > 0 and ip > e.calls[e.call_top - 1].def.end)) return error.BadArgument;
    e.ip = @intCast(ip);
    e.length = 0;
    if (off < 0) {
        e.neg_jump_counter += 1;
        if (e.neg_jump_counter > e.loop_max) return error.ExecutionTooLong;
    }
}

/// JROT and JROF: JMPR when `args[1]` is (not) 0.
pub fn jumpIf(e: *Engine, comptime when: bool) Error!void {
    if ((e.arg(1).* != 0) == when) try jump(e);
}
