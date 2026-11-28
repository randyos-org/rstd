//! Defining and calling functions and instructions: FreeType 2.14.3
//! `ttinterp.c` `Ins_FDEF`, `Ins_ENDF`, `Ins_CALL`, `Ins_LOOPCALL`,
//! `Ins_IDEF`, `Ins_UNKNOWN`.

const definitions = @import("../definitions.zig");
const Engine = @import("../Engine.zig");

const Error = Engine.Error;
const Definition = definitions.Definition;

/// FDEF (`instruction` false) or IDEF: records the code up to ENDF.
pub fn define(e: *Engine, comptime instruction: bool) Error!void {
    if (e.ini == .glyph) return error.DefinitionInGlyphProgram;
    const t = if (instruction) &e.instructions else &e.functions;
    const n: i64 = e.arg(0).*;
    // a negative number matches none, as FreeType's unsigned compare
    const key: u64 = @bitCast(n);
    const d = t.slot(key) orelse return error.TooManyDefinitions;
    if (n < 0 or n > (if (instruction) 0xFF else 0xFFFF)) return error.TooManyDefinitions;
    const number: u32 = @intCast(n);
    d.* = .{ .program = e.cur, .start = e.ip + 1, .number = number, .active = true };
    t.max_number = @max(t.max_number, number);
    while (true) {
        try e.skipCode();
        switch (e.op) {
            0x89, 0x2C => return error.NestedDefinition,
            0x2D => {
                d.end = e.ip;
                return;
            },
            else => {},
        }
    }
}

/// ENDF: the next LOOPCALL iteration, else back to the caller.
pub fn endf(e: *Engine) Error!void {
    if (e.call_top == 0) return error.EndfWithoutCall;
    e.call_top -= 1;
    const c = &e.calls[e.call_top];
    c.count -= 1;
    if (c.count > 0) {
        e.call_top += 1;
        e.ip = c.def.start;
        e.length = 0;
    } else try e.goTo(c.caller, c.return_ip);
}

/// The active function `f`, as `Ins_CALL` finds it.
fn function(e: *Engine, f: i32) Error!*const Definition {
    const key: u64 = @bitCast(@as(i64, f));
    if (key > e.functions.max_number) return error.InvalidDefinition;
    const d = e.functions.find(key) orelse return error.InvalidDefinition;
    if (!d.active) return error.InvalidDefinition;
    return d;
}

fn enter(e: *Engine, d: *const Definition, count: i32) Error!void {
    e.calls[e.call_top] = .{ .caller = e.cur, .return_ip = e.ip + 1, .count = count, .def = d };
    e.call_top += 1;
    try e.goTo(d.program, d.start);
}

pub fn call(e: *Engine) Error!void {
    const d = try function(e, e.arg(0).*);
    if (e.call_top >= Engine.call_depth) return error.StackOverflow;
    try enter(e, d, 1);
}

/// LOOPCALL: `args[1]` called `args[0]` times; none for a count of 0 or less.
pub fn loopcall(e: *Engine) Error!void {
    const d = try function(e, e.arg(1).*);
    if (e.call_top >= Engine.call_depth) return error.StackOverflow;
    const count = e.arg(0).*;
    if (count <= 0) return;
    try enter(e, d, count);
    e.loopcall_counter += @intCast(count);
    if (e.loopcall_counter > e.loop_max) return error.ExecutionTooLong;
}

/// An opcode with no instruction of its own: its IDEF, if any.
pub fn unknown(e: *Engine) Error!void {
    for (e.instructions.defs[0..e.instructions.count]) |*d| {
        if (d.active and @as(u8, @truncate(d.number)) == e.op) {
            if (e.call_top >= Engine.call_depth) return error.StackOverflow;
            return enter(e, d, 1);
        }
    }
    return error.UnhandledOpcode;
}
