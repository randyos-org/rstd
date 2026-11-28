//! Type 2 charstrings evaluated in 16.16 as FreeType 2.14.3's Adobe engine
//! evaluates them (`src/psaux/psintrp.c`), after fontations
//! `read-fonts/src/ps/charstring.rs` and `cff/stack.rs` (commit 4f115482,
//! MIT): commands go to a sink with `moveTo`, `lineTo`, `curveTo`,
//! `close`, `hstem`, `vstem`, `hintMask`, `counterMask` and
//! `closeOpenPath`. A charstring ending in an implied SEAC returns its
//! components for the caller to draw (`cf2_cmdENDCHAR`).

const std = @import("std");
const fixed = @import("../fixed.zig");
const Index = @import("../../sfnt/tables/cff/Index.zig");

pub const Error = error{ InvalidCharstring, StackOverflow, StackUnderflow, LimitExceeded };

const nesting_limit = 10;
const max_operations = 1_000_000;
const max_stack = 513;
const storage_len = 32;

/// The subroutines a glyph can call.
pub const Context = struct {
    global: Index,
    local: ?Index,
};

/// The bias added to a subroutine number (Type 2 spec, section 4.7).
fn bias(count: usize) i32 {
    return if (count < 1240) 107 else if (count < 33900) 1131 else 32768;
}

/// The operand stack: integers or 16.16 values (`stack.rs`).
const Stack = struct {
    values: [max_stack]i32 = undefined,
    is_fixed: [max_stack]bool = undefined,
    top: usize = 0,

    fn push(s: *Stack, v: i32, is_fixed: bool) Error!void {
        if (s.top >= max_stack) return error.StackOverflow;
        s.values[s.top] = v;
        s.is_fixed[s.top] = is_fixed;
        s.top += 1;
    }

    fn get(s: *const Stack, i: usize) i32 {
        if (i >= s.top) return 0;
        return if (s.is_fixed[i]) s.values[i] else s.values[i] *% 65536;
    }

    fn popFixed(s: *Stack) Error!i32 {
        if (s.top == 0) return error.StackUnderflow;
        s.top -= 1;
        const v = s.values[s.top];
        return if (s.is_fixed[s.top]) v else v *% 65536;
    }

    /// ROLL: the top `count` values move `shift` places toward the top.
    fn roll(s: *Stack, count: i32, shift: i32) Error!void {
        if (count < 2) return;
        if (count > s.top) return error.StackUnderflow;
        const n: usize = @intCast(count);
        const r: usize = @intCast(@mod(@as(i64, shift), count));
        if (r == 0) return;
        std.mem.rotate(i32, s.values[s.top - n .. s.top], n - r);
        std.mem.rotate(bool, s.is_fixed[s.top - n .. s.top], n - r);
    }

    fn popInt(s: *Stack) Error!i32 {
        if (s.top == 0) return error.StackUnderflow;
        s.top -= 1;
        if (s.is_fixed[s.top]) return error.InvalidCharstring;
        return s.values[s.top];
    }
};

const Mode = enum { dxdy, xdy, dxy, dx_initial_y, larger_coord_dist, dx_maybe_dy, maybe_dx_dy };

/// An implied SEAC: an accent at (`adx`, `ady`), 16.16, over a base, both
/// by Standard Encoding code.
pub const Seac = struct { adx: i32, ady: i32, base: i32, accent: i32 };

pub const Result = struct {
    /// The width the charstring gives, if any (before `nominalWidthX`).
    width: ?i32,
    seac: ?Seac = null,
};

/// Evaluates `charstring` into `sink`, from `start` (16.16).
pub fn evaluate(ctx: Context, charstring: []const u8, sink: anytype, start: [2]i32) Error!Result {
    var e: Evaluator(@TypeOf(sink)) = .{ .ctx = ctx, .sink = sink, .x = start[0], .y = start[1] };
    const seen_endchar = try e.run(charstring, 0);
    if (!seen_endchar) _ = try e.endchar();
    if (e.is_open) e.sink.close();
    e.sink.closeOpenPath();
    return .{ .width = if (e.have_width) e.width else null, .seac = e.seac };
}

fn Evaluator(comptime Sink: type) type {
    return struct {
        const Self = @This();

        ctx: Context,
        sink: Sink,
        stack: Stack = .{},
        stack_ix: usize = 0,
        x: i32 = 0,
        y: i32 = 0,
        is_open: bool = false,
        seen_width: bool = false,
        have_width: bool = false,
        width: i32 = 0,
        stem_count: usize = 0,
        ops: u32 = 0,
        storage: [storage_len]i32 = @splat(0),
        seac: ?Seac = null,

        fn reset(e: *Self) void {
            e.stack.top = 0;
            e.stack_ix = 0;
        }

        fn readWidth(e: *Self) void {
            e.width = e.stack.get(0);
            e.seen_width = true;
            e.have_width = true;
        }

        fn ensureOpen(e: *Self) void {
            if (e.is_open) return;
            e.sink.moveTo(0, 0);
            e.is_open = true;
        }

        fn line(e: *Self, x: i32, y: i32) void {
            e.ensureOpen();
            e.sink.lineTo(x, y);
        }

        fn moveTo(e: *Self) void {
            if (e.is_open) e.sink.close() else e.is_open = true;
            e.sink.moveTo(e.x, e.y);
        }

        fn curves(e: *Self, modes: []const Mode) void {
            const x0 = e.x;
            const y0 = e.y;
            var count: usize = 0;
            var pts: [2][2]i32 = undefined;
            e.ensureOpen();
            for (modes) |m| {
                const s = &e.stack;
                const i = e.stack_ix;
                switch (m) {
                    .dxdy => {
                        e.x +%= s.get(i);
                        e.y +%= s.get(i + 1);
                        e.stack_ix += 2;
                    },
                    .xdy => {
                        e.y +%= s.get(i);
                        e.stack_ix += 1;
                    },
                    .dxy => {
                        e.x +%= s.get(i);
                        e.stack_ix += 1;
                    },
                    .dx_initial_y => {
                        e.x +%= s.get(i);
                        e.y = y0;
                        e.stack_ix += 1;
                    },
                    .larger_coord_dist => {
                        const d = s.get(i);
                        if (@abs(e.x -% x0) > @abs(e.y -% y0)) {
                            e.x +%= d;
                            e.y = y0;
                        } else {
                            e.y +%= d;
                            e.x = x0;
                        }
                        e.stack_ix += 1;
                    },
                    .dx_maybe_dy, .maybe_dx_dy => unreachable,
                }
                if (count == 2) {
                    e.sink.curveTo(pts[0][0], pts[0][1], pts[1][0], pts[1][1], e.x, e.y);
                    count = 0;
                } else {
                    pts[count] = .{ e.x, e.y };
                    count += 1;
                }
            }
        }

        /// HVCURVETO and VHCURVETO's curve: its last delta in the other
        /// coordinate too when `last`.
        fn alternating(e: *Self, horizontal: bool, last: bool) void {
            const s = &e.stack;
            e.ensureOpen();
            const i = e.stack_ix;
            if (horizontal) e.x +%= s.get(i) else e.y +%= s.get(i);
            const c0: [2]i32 = .{ e.x, e.y };
            e.x +%= s.get(i + 1);
            e.y +%= s.get(i + 2);
            const c1: [2]i32 = .{ e.x, e.y };
            if (horizontal) {
                e.y +%= s.get(i + 3);
                if (last) e.x +%= s.get(i + 4);
            } else {
                e.x +%= s.get(i + 3);
                if (last) e.y +%= s.get(i + 4);
            }
            e.stack_ix += if (last) 5 else 4;
            e.sink.curveTo(c0[0], c0[1], c1[0], c1[1], e.x, e.y);
        }

        fn stems(e: *Self, horizontal: bool) void {
            var i: usize = 0;
            var len = e.stack.top;
            if (len & 1 != 0 and !e.seen_width) {
                e.readWidth();
                i = 1;
                len -= 1;
            }
            e.seen_width = true;
            var u: i32 = 0;
            while (i < e.stack.top) : (i += 2) {
                u +%= e.stack.get(i);
                const v = u +% e.stack.get(i + 1);
                if (horizontal) e.sink.hstem(u, v) else e.sink.vstem(u, v);
                u = v;
            }
            e.stem_count += len / 2;
            e.reset();
        }

        fn endchar(e: *Self) Error!bool {
            const n = e.stack.top;
            if ((n == 1 or n == 5) and !e.seen_width) e.readWidth();
            e.seen_width = true;
            // an implied SEAC builds the glyph from two others
            if (n > 1) {
                const accent = try e.stack.popInt();
                const base = try e.stack.popInt();
                const ady = try e.stack.popFixed();
                const adx = try e.stack.popFixed();
                e.seac = .{ .adx = adx, .ady = ady, .base = base, .accent = accent };
            }
            return false;
        }

        /// Runs `code`; true if it reached ENDCHAR.
        fn run(e: *Self, code: []const u8, depth: u32) Error!bool {
            if (depth > nesting_limit) return error.LimitExceeded;
            var at: usize = 0;
            var seen_endchar = false;
            while (at < code.len) {
                e.ops += 1;
                if (e.ops > max_operations) return error.LimitExceeded;
                const b0 = code[at];
                at += 1;
                switch (b0) {
                    28, 32...254 => try e.stack.push(try int(code, &at, b0), false),
                    255 => {
                        if (at + 4 > code.len) return error.InvalidCharstring;
                        try e.stack.push(std.mem.readInt(i32, code[at..][0..4], .big), true);
                        at += 4;
                    },
                    else => {
                        var op: u16 = b0;
                        if (b0 == 12) {
                            if (at >= code.len) return error.InvalidCharstring;
                            op = 0x0C00 | @as(u16, code[at]);
                            at += 1;
                        }
                        if (op == 14) seen_endchar = true;
                        if (!try e.operator(op, code, &at, depth)) break;
                    },
                }
            }
            return seen_endchar;
        }

        /// One operator; false to stop the current charstring.
        fn operator(e: *Self, op: u16, code: []const u8, at: *usize, depth: u32) Error!bool {
            const s = &e.stack;
            switch (op) {
                1, 3, 18, 23 => e.stems(op == 1 or op == 18),
                19, 20 => {
                    // HINTMASK and CNTRMASK: implied vstems first
                    var i: usize = 0;
                    var len = s.top;
                    if (len & 1 != 0 and !e.seen_width) {
                        e.readWidth();
                        i = 1;
                        len -= 1;
                    }
                    e.seen_width = true;
                    var u: i32 = 0;
                    while (i < s.top) : (i += 2) {
                        u +%= s.get(i);
                        const v = u +% s.get(i + 1);
                        e.sink.vstem(u, v);
                        u = v;
                    }
                    e.stem_count += len / 2;
                    const count = (e.stem_count + 7) / 8;
                    if (at.* + count > code.len) return error.InvalidCharstring;
                    const mask = code[at.*..][0..count];
                    at.* += count;
                    if (op == 19) e.sink.hintMask(mask) else e.sink.counterMask(mask);
                    e.reset();
                },
                21 => {
                    if (s.top > 2 and !e.seen_width) e.readWidth();
                    e.seen_width = true;
                    const dy = try s.popFixed();
                    const dx = try s.popFixed();
                    e.x +%= dx;
                    e.y +%= dy;
                    e.moveTo();
                    e.reset();
                },
                22, 4 => {
                    if (s.top > 1 and !e.seen_width) e.readWidth();
                    e.seen_width = true;
                    const d = try s.popFixed();
                    if (op == 22) e.x +%= d else e.y +%= d;
                    e.moveTo();
                    e.reset();
                },
                5 => {
                    var i: usize = 0;
                    while (i < s.top) : (i += 2) {
                        e.x +%= s.get(i);
                        e.y +%= s.get(i + 1);
                        e.line(e.x, e.y);
                    }
                    e.reset();
                },
                6, 7 => {
                    var horizontal = op == 6;
                    for (0..s.top) |i| {
                        if (horizontal) e.x +%= s.get(i) else e.y +%= s.get(i);
                        horizontal = !horizontal;
                        e.line(e.x, e.y);
                    }
                    e.reset();
                },
                27, 26 => {
                    // HHCURVETO, VVCURVETO: an odd first value moves the other way
                    const count1 = s.top;
                    const count = count1 & ~@as(usize, 2);
                    e.stack_ix = count1 - count;
                    while (e.stack_ix < count) {
                        if ((count - e.stack_ix) & 1 != 0) {
                            if (op == 27) e.y +%= s.get(e.stack_ix) else e.x +%= s.get(e.stack_ix);
                            e.stack_ix += 1;
                        }
                        if (op == 27) e.curves(&.{ .dxy, .dxdy, .dxy }) else e.curves(&.{ .xdy, .dxdy, .xdy });
                    }
                    e.reset();
                },
                31, 30 => {
                    const count1 = s.top;
                    const count = count1 & ~@as(usize, 2);
                    var horizontal = op == 31;
                    e.stack_ix = count1 - count;
                    while (e.stack_ix < count) {
                        e.alternating(horizontal, count - e.stack_ix == 5);
                        horizontal = !horizontal;
                    }
                    e.reset();
                },
                8, 24 => {
                    while (s.top -| e.stack_ix >= 6) e.curves(&.{ .dxdy, .dxdy, .dxdy });
                    if (op == 24) {
                        e.x +%= s.get(e.stack_ix);
                        e.y +%= s.get(e.stack_ix + 1);
                        e.line(e.x, e.y);
                    }
                    e.reset();
                },
                25 => {
                    while (s.top -| e.stack_ix > 6) {
                        e.x +%= s.get(e.stack_ix);
                        e.y +%= s.get(e.stack_ix + 1);
                        e.line(e.x, e.y);
                        e.stack_ix += 2;
                    }
                    while (s.top -| e.stack_ix >= 6) e.curves(&.{ .dxdy, .dxdy, .dxdy });
                    e.reset();
                },
                0x0C23 => {
                    e.curves(&.{ .dxdy, .dxdy, .dxdy, .dxdy, .dxdy, .dxdy });
                    e.reset();
                },
                0x0C22 => {
                    e.curves(&.{ .dxy, .dxdy, .dxy, .dxy, .dx_initial_y, .dxy });
                    e.reset();
                },
                0x0C24 => {
                    e.curves(&.{ .dxdy, .dxdy, .dxy, .dxy, .dxdy, .dx_initial_y });
                    e.reset();
                },
                0x0C25 => {
                    e.curves(&.{ .dxdy, .dxdy, .dxdy, .dxdy, .dxdy, .larger_coord_dist });
                    e.reset();
                },
                10, 29 => {
                    const n = try s.popInt();
                    const idx: ?Index = if (op == 10) e.ctx.local else e.ctx.global;
                    const subrs = idx orelse return error.InvalidCharstring;
                    const k = n +% bias(subrs.count());
                    if (k < 0 or k >= subrs.count()) return error.InvalidCharstring;
                    _ = try e.run(subrs.get(@intCast(k)), depth + 1);
                },
                11 => return false,
                14 => return e.endchar(),
                0x0C03, 0x0C04, 0x0C0F => {
                    const b = try s.popFixed();
                    const a = try s.popFixed();
                    const r = switch (op) {
                        0x0C03 => a != 0 and b != 0,
                        0x0C04 => a != 0 or b != 0,
                        else => a == b,
                    };
                    try s.push(@intFromBool(r), false);
                },
                0x0C05 => try s.push(@intFromBool(try s.popFixed() == 0), false),
                0x0C09, 0x0C0E => {
                    const a = try s.popFixed();
                    const r: i32 = if (a < -std.math.maxInt(i32)) std.math.maxInt(i32) else if (op == 0x0C09) @intCast(@abs(a)) else -a;
                    try s.push(r, true);
                },
                0x0C0A, 0x0C0B => {
                    const b = try s.popFixed();
                    const a = try s.popFixed();
                    try s.push(if (op == 0x0C0A) a +% b else a -% b, true);
                },
                0x0C0C => {
                    const b = try s.popFixed();
                    const a = try s.popFixed();
                    try s.push(fixed.div(a, b), true);
                },
                0x0C1A => try s.push(sqrt(try s.popFixed()), true),
                0x0C1E => {
                    const shift = try s.popInt();
                    const count = try s.popInt();
                    try s.roll(count, shift);
                },
                0x0C18 => {
                    const b = try s.popFixed();
                    const a = try s.popFixed();
                    try s.push(fixed.mul(a, b), true);
                },
                0x0C12 => _ = try s.popFixed(),
                0x0C1B => {
                    const a = try s.popFixed();
                    try s.push(a, true);
                    try s.push(a, true);
                },
                0x0C1C => {
                    if (s.top < 2) return error.StackUnderflow;
                    std.mem.swap(i32, &s.values[s.top - 1], &s.values[s.top - 2]);
                    std.mem.swap(bool, &s.is_fixed[s.top - 1], &s.is_fixed[s.top - 2]);
                },
                0x0C1D => {
                    const i = try s.popInt();
                    if (s.top > 0) {
                        const n = s.top;
                        const from: usize = if (i < 0) n - 1 else if (i >= n) 0 else n - 1 - @as(usize, @intCast(i));
                        try s.push(s.get(from), true);
                    }
                },
                0x0C14 => {
                    const i = try s.popInt();
                    const v = try s.popFixed();
                    if (i >= 0 and i < storage_len) e.storage[@intCast(i)] = v;
                },
                0x0C15 => {
                    const i = try s.popInt();
                    if (i >= 0 and i < storage_len) try s.push(e.storage[@intCast(i)], true);
                },
                0x0C16 => {
                    const c2 = try s.popFixed();
                    const c1 = try s.popFixed();
                    const a2 = try s.popFixed();
                    const a1 = try s.popFixed();
                    try s.push(if (c1 <= c2) a1 else a2, true);
                },
                // HSBW, DOTSECTION, SBW, POP, SETCURRENTPOINT: nothing in Type 2
                13, 0x0C00, 0x0C07, 0x0C11, 0x0C21 => {},
                // VSINDEX and BLEND belong to CFF2
                15, 16 => return error.InvalidCharstring,
                // the reserved operators, HSTEM3, VSTEM3, RANDOM and Type 1's SEAC clear the stack
                else => e.reset(),
            }
            return true;
        }
    };
}

/// The 16.16 square root, Babylonian, as `cf2_interpT2CharString` takes it.
fn sqrt(v: i32) i32 {
    if (v <= 0) return 0;
    var root = if (v < 10) v else v >> 1;
    while (true) {
        const next = (root +% fixed.div(v, root) +% 1) >> 1;
        if (next == root) return next;
        root = next;
    }
}

fn int(code: []const u8, at: *usize, b0: u8) Error!i32 {
    switch (b0) {
        32...246 => return @as(i32, b0) - 139,
        247...250 => {
            if (at.* >= code.len) return error.InvalidCharstring;
            defer at.* += 1;
            return (@as(i32, b0) - 247) * 256 + code[at.*] + 108;
        },
        251...254 => {
            if (at.* >= code.len) return error.InvalidCharstring;
            defer at.* += 1;
            return -(@as(i32, b0) - 251) * 256 - @as(i32, code[at.*]) - 108;
        },
        else => {
            if (at.* + 2 > code.len) return error.InvalidCharstring;
            defer at.* += 2;
            return std.mem.readInt(i16, code[at.*..][0..2], .big);
        },
    }
}
