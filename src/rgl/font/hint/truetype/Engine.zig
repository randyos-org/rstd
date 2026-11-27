//! The TrueType bytecode interpreter: FreeType 2.14.3 `src/truetype/ttinterp.c`
//! (`TT_RunIns`, `TT_Run_Context`), its execution model kept: the stack
//! checked against each instruction's counts before it runs (short stacks
//! zero-filled unless pedantic), any error ending the program, glyph
//! programs starting from the state `prep` left. Laid out after fontations
//! `skrifa/src/outline/glyf/hint/engine/` (commit 4f115482, MIT).
//!
//! The caller owns every buffer, sized from `maxp` (the per-size setup,
//! #808). Values are `i32`: FreeType's `FT_Long` is
//! as wide only on 32-bit hosts, so the two can differ once a value passes
//! 32 bits, which no real font's program does.

const std = @import("std");
const opcode = @import("opcode.zig");
const definitions = @import("definitions.zig");
const program_mod = @import("program.zig");
const zone = @import("zone.zig");
const round = @import("round.zig");
const errors = @import("errors.zig");
const dispatch = @import("dispatch.zig");
const GraphicsState = @import("GraphicsState.zig");
const Projection = @import("Projection.zig");
const CowSlice = @import("CowSlice.zig");

const Engine = @This();

pub const Error = errors.Error;
pub const Program = program_mod.Program;
pub const Zone = zone.Zone;

/// FreeType's `interpreter-version` property.
pub const Version = enum { v35, v40 };

/// The render mode the load targets (`FT_LOAD_TARGET_*`).
pub const Target = enum { normal, light, mono, lcd, lcd_v };

/// The size's metrics; square pixels only.
pub const Metrics = struct {
    ppem: i32,
    /// Font units to 26.6 pixels, 16.16.
    scale: i32,
    /// The size in 26.6 points, `MPS`'s answer under v40.
    point_size: i32,
    rotated: bool = false,
    stretched: bool = false,
};

/// One active call (`TT_CallRec`).
pub const Call = struct {
    caller: Program,
    /// Where the caller resumes.
    return_ip: usize,
    /// LOOPCALL's iterations left, this one included.
    count: i32,
    def: *const definitions.Definition,
};

/// FreeType's call stack depth (`ttinterp.c` `Init_Context`).
pub const call_depth = 32;

/// `TT_CONFIG_OPTION_MAX_RUNNABLE_OPCODES`, FreeType's default.
pub const max_instructions = 1_000_000;

version: Version = .v40,
target: Target = .normal,
pedantic: bool = false,
metrics: Metrics,
/// The face's glyph count, which caps the loop limits.
glyph_count: u32,
/// Variation axes; 0 unless the font is variable (`face->blend`).
axis_count: u16 = 0,
/// Normalized coordinates in 16.16; missing ones read as 0.
coords: []const i32 = &.{},
/// A composite glyph, whose `orus` are already scaled: MD measures them at
/// 1.0 (fontations `graphics.rs` `unscaled_to_pixels`).
composite: bool = false,

stack: []i32,
functions: definitions.Table,
instructions: definitions.Table,
cvt: CowSlice,
storage: CowSlice,
zones: [2]Zone = .{ .{}, .{} },
/// Each program's bytecode, by `Program`.
code: [3][]const u8 = .{ &.{}, &.{}, &.{} },

/// The graphics state, which persists between runs.
gs: GraphicsState = .{},
/// v40's backward compatibility: 4 on, plus 1 and 2 once IUP has run in x
/// and y (`ttinterp.h`); the loader sets bit 4.
backward_compatibility: u8 = 0,
super_round: round.Super = .{},

// one run's state
ini: Program = .font,
cur: Program = .font,
ip: usize = 0,
/// The current instruction's length; 0 once it jumped.
length: usize = 1,
op: u8 = 0,
top: usize = 0,
/// Where the current instruction's arguments start.
args: usize = 0,
/// The stack's depth once it is done.
new_top: usize = 0,
calls: [call_depth]Call = undefined,
call_top: usize = 0,
zp: [3]zone.Pointer = @splat(.glyph),
projection: Projection = .compute(.{}),
round_mode: round.Mode = .grid,
loopcall_counter: u64 = 0,
neg_jump_counter: u64 = 0,
loop_max: u64 = 0,

/// Argument `i` of the current instruction.
pub fn arg(e: *Engine, i: usize) *i32 {
    return &e.stack[e.args + i];
}

pub fn zoneAt(e: *Engine, which: zone.Pointer) *Zone {
    return &e.zones[@backingInt(which)];
}

/// Zone pointer `i`'s zone (`zp0` to `zp2`).
pub fn zpZone(e: *Engine, i: usize) *Zone {
    return e.zoneAt(e.zp[i]);
}

/// Runs `program` from its start (`TT_Run_Context`, then `TT_RunIns`):
/// zone pointers on the glyph zone, rounding to the grid, the loop limits
/// from the glyph's points and the CVT. The font program clears the
/// definitions first, as `tt_size_init_bytecode` does.
pub fn run(e: *Engine, program: Program) Error!void {
    if (program == .font) {
        e.functions.reset();
        e.instructions.reset();
    }
    e.ini = program;
    e.cur = program;
    e.ip = 0;
    e.cvt.copied = false;
    e.storage.copied = false;
    e.zp = @splat(.glyph);
    // twilight points past 2 (points + CVT), at least 30, go unused for the
    // run: slow malformed bytecode otherwise (ttfautohint stays inside)
    const twilight = e.zones[0];
    defer e.zones[0] = twilight;
    const tw_max: usize = @min(0xFFFF, @max(30, 2 * (e.zones[1].len() + e.cvt.len())));
    if (twilight.len() > tw_max) {
        e.zones[0].org = twilight.org[0..tw_max];
        e.zones[0].cur = twilight.cur[0..tw_max];
        e.zones[0].tags = twilight.tags[0..tw_max];
    }
    e.loopcall_counter = 0;
    e.neg_jump_counter = 0;
    const points: u64 = e.zones[1].len();
    const cvt_len: u64 = e.cvt.len();
    e.loop_max = if (points > 0) @max(50, 10 * points) + @max(50, cvt_len / 10) else 300 + 22 * cvt_len;
    e.loop_max = @min(e.loop_max, 100 * @as(u64, e.glyph_count));
    e.round_mode = .grid;
    e.projection = .compute(e.gs);
    e.backward_compatibility &= ~@as(u8, 3);
    e.top = 0;
    e.call_top = 0;
    if (e.code[@backingInt(program)].len == 0) return;
    var count: u32 = 0;
    while (true) {
        count += 1;
        if (count > max_instructions) return error.ExecutionTooLong;
        e.op = e.code[@backingInt(e.cur)][e.ip];
        e.length = 1;
        const fx = opcode.stack_effects[e.op];
        if (e.top < fx.pop) {
            if (e.pedantic) return error.TooFewArguments;
            if (e.stack.len < fx.pop) return error.StackOverflow;
            @memset(e.stack[0..fx.pop], 0);
            e.args = 0;
        } else e.args = e.top - fx.pop;
        e.new_top = e.args + fx.push;
        if (e.new_top > e.stack.len) return error.StackOverflow;
        try dispatch.step(e);
        e.top = e.new_top;
        e.ip += e.length;
        if (e.ip >= e.code[@backingInt(e.cur)].len) {
            if (e.call_top > 0) return error.CodeOverflow;
            return;
        }
    }
}

/// Continues in `program` at `ip` (`Ins_Goto_CodeRange`).
pub fn goTo(e: *Engine, program: Program, ip: usize) Error!void {
    const code = e.code[@backingInt(program)];
    if (code.len == 0) return error.InvalidCodeRange;
    if (ip > code.len) return error.CodeOverflow;
    e.cur = program;
    e.ip = ip;
    e.length = 0;
}

/// The next instruction past the current one, as the skips walk
/// (`SkipCode`): sets `op` and `length`.
pub fn skipCode(e: *Engine) Error!void {
    const code = e.code[@backingInt(e.cur)];
    e.ip += e.length;
    if (e.ip >= code.len) return error.CodeOverflow;
    e.op = code[e.ip];
    const l = opcode.lengths[e.op];
    if (l >= 0) {
        e.length = @intCast(l);
        return;
    }
    if (e.ip + 1 >= code.len) return error.CodeOverflow;
    e.length = 2 + @as(usize, @abs(l)) * code[e.ip + 1];
}

/// `d` rounded by the current mode (`func_round`).
pub fn roundDistance(e: *const Engine, d: i32) i32 {
    return round.round(e.round_mode, e.super_round, d);
}

/// FreeType's `pedantic_hinting` checks: `err` when pedantic.
pub fn pedanticError(e: *const Engine, err: Error) Error!void {
    if (e.pedantic) return err;
}

/// Font units to 26.6 for the outline (`metrics.x_scale`): 1.0 for a
/// composite, whose `orus` are already scaled.
pub fn unitsScale(e: *const Engine) i32 {
    return if (e.composite) 0x10000 else e.metrics.scale;
}
