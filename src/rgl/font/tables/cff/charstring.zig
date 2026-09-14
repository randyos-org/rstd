const std = @import("std");
const assert = std.debug.assert;

const zigstd = @import("../../../../zigstd.zig");
const types = @import("../../stream_types.zig");
const aux_types = @import("../../aux_types.zig");
const CffTable = @import("CffTable.zig");
const Index = @import("Index.zig");

const Allocator = zigstd.Allocator;
const Vertex = aux_types.Vertex;
const GlyphIndex = types.GlyphIndex;

pub const CharstringCtx = struct {
    first_x: f32,
    first_y: f32,
    x: f32,
    y: f32,
    min_x: i32,
    min_y: i32,
    max_x: i32,
    max_y: i32,
    num_vertices: u32,
    vertices: [*]Vertex,
    flags: Flags,

    pub const Flags = packed struct(u8) {
        started: bool = false,
        mode: enum(u1) {
            /// set min/max and num_vertices
            bounds,
            /// set vertices and num_vertices
            verts,
        },
        _padding: u6 = undefined,
    };

    pub fn init(flags: Flags, vertices: [*]Vertex) CharstringCtx {
        return .{
            .flags = flags,
            .vertices = vertices,
            .first_x = 0,
            .first_y = 0,
            .x = 0,
            .y = 0,
            .min_x = 0,
            .min_y = 0,
            .max_x = 0,
            .max_y = 0,
            .num_vertices = 0,
        };
    }
    pub fn deinit(ctx: *CharstringCtx, alloc: Allocator) void {
        if (ctx.flags.mode == .verts)
            alloc.free(ctx.allVertices());
    }

    fn trackVertex(ctx: *CharstringCtx, x: i32, y: i32) void {
        if (x > ctx.max_x or !ctx.flags.started) ctx.max_x = x;
        if (y > ctx.max_y or !ctx.flags.started) ctx.max_y = y;
        if (x < ctx.min_x or !ctx.flags.started) ctx.min_x = x;
        if (y < ctx.min_y or !ctx.flags.started) ctx.min_y = y;
        ctx.flags.started = true;
    }

    fn v(ctx: *CharstringCtx, ty: Vertex.Type, x: i32, y: i32, cx: i32, cy: i32, cx1: i32, cy1: i32) !void {
        if (ctx.flags.mode == .bounds) {
            trackVertex(ctx, x, y);
            if (ty == .vcubic) {
                trackVertex(ctx, cx, cy);
                trackVertex(ctx, cx1, cy1);
            }
        } else {
            ctx.vertices[ctx.num_vertices].set(ty, x, y, cx, cy);
            ctx.vertices[ctx.num_vertices].cx1 = @truncate(cx1);
            ctx.vertices[ctx.num_vertices].cy1 = @truncate(cy1);
        }
        ctx.num_vertices += 1;
    }

    fn closeShape(ctx: *CharstringCtx) !void {
        if (ctx.first_x != ctx.x or ctx.first_y != ctx.y)
            try ctx.v(.vline, @trunc(ctx.first_x), @trunc(ctx.first_y), 0, 0, 0, 0);
    }

    fn rmoveTo(ctx: *CharstringCtx, dx: f32, dy: f32) !void {
        try ctx.closeShape();
        ctx.first_x = ctx.x + dx;
        ctx.x = ctx.first_x;
        ctx.first_y = ctx.y + dy;
        ctx.y = ctx.first_y;
        try ctx.v(.vmove, @trunc(ctx.x), @trunc(ctx.y), 0, 0, 0, 0);
    }

    fn rlineTo(ctx: *CharstringCtx, dx: f32, dy: f32) !void {
        ctx.x += dx;
        ctx.y += dy;
        try ctx.v(.vline, @trunc(ctx.x), @trunc(ctx.y), 0, 0, 0, 0);
    }

    fn rccurveTo(ctx: *CharstringCtx, dx1: f32, dy1: f32, dx2: f32, dy2: f32, dx3: f32, dy3: f32) !void {
        const cx1 = ctx.x + dx1;
        const cy1 = ctx.y + dy1;
        const cx2 = cx1 + dx2;
        const cy2 = cy1 + dy2;
        ctx.x = cx2 + dx3;
        ctx.y = cy2 + dy3;
        try ctx.v(
            .vcubic,
            @trunc(ctx.x),
            @trunc(ctx.y),
            @trunc(cx1),
            @trunc(cy1),
            @trunc(cx2),
            @trunc(cy2),
        );
    }
};

/// Errors from the Type 2 charstring interpreter itself (malformed operand
/// stacks, unimplemented operators, a missing subroutine, ...) - distinct
/// from CffTable.CffError (which is about resolving the table structure
/// itself) and from Allocator.Error, which callers union in separately.
pub const CharstringError = error{
    Unimplemented,
    RMoveToStack,
    VMoveToStack,
    HMoveToStack,
    RLineToStack,
    VLineToStack,
    HLineToStack,
    HCurveToStack,
    RCurveToStack,
    RCurveLineStack,
    CurveLineStack,
    RLineCurveStack,
    CallGSubRStack,
    RecursionLimit,
    SubRNotFound,
    ReturnOutsideSubR,
    HFlexStack,
    FlexStack,
    HFlex1Stack,
    Flex1Stack,
    CurveToStack,
    ReservedOperator,
    PushStackOverflow,
    NoEndChar,
};

pub fn glyphBoxT2(cff: CffTable, glyph: GlyphIndex) ?aux_types.GlyphBoundingBox {
    var ctx = CharstringCtx.init(.{ .mode = .bounds }, undefined);
    runCharstring(cff, glyph, &ctx) catch return null;

    return .{
        .x_min = ctx.min_x,
        .y_min = ctx.min_y,
        .x_max = ctx.max_x,
        .y_max = ctx.max_y,
    };
}

pub fn glyphShapeT2(cff: CffTable, gpa: Allocator, glyph: GlyphIndex) (Allocator.Error || CharstringError)![]Vertex {
    // mode=bounds to get bounds and num_vertices
    var count_ctx = CharstringCtx.init(.{ .mode = .bounds }, undefined);
    try runCharstring(cff, glyph, &count_ctx);
    const vertices = try gpa.alloc(Vertex, count_ctx.num_vertices);
    errdefer gpa.free(vertices);
    // mode=verts to assign vertices
    var out_ctx = CharstringCtx.init(.{ .mode = .verts }, vertices.ptr);
    try runCharstring(cff, glyph, &out_ctx);
    assert(out_ctx.num_vertices == count_ctx.num_vertices);

    return out_ctx.vertices[0..out_ctx.num_vertices];
}

const Instruction = enum(u8) {
    hintmask = 0x13,
    cntrmask = 0x14,
    hstem = 0x01,
    vstem = 0x03,
    hstemhm = 0x12,
    vstemhm = 0x17,
    rmoveto = 0x15,
    vmoveto = 0x04,
    hmoveto = 0x16,
    rlineto = 0x05,
    vlineto = 0x07,
    hlineto = 0x06,
    hvcurveto = 0x1F,
    vhcurveto = 0x1E,
    rrcurveto = 0x08,
    rcurveline = 0x18,
    rlinecurve = 0x19,
    vvcurveto = 0x1A,
    hhcurveto = 0x1B,
    callsubr = 0x0A,
    callgsubr = 0x1D,
    /// return
    ret = 0x0B,
    endchar = 0x0E,
    twoByteEscape = 0x0C,
    hflex = 0x22,
    flex = 0x23,
    hflex1 = 0x24,
    flex1 = 0x25,

    pub fn toInt(i: Instruction) u16 {
        return @backingInt(i);
    }
};

/// A cursor into one charstring or subroutine's bytes, for the interpreter's
/// own bytecode walk. Index/Dict/CffTable already own INDEX- and DICT-
/// structure parsing (with real bounds checking); this is deliberately
/// slimmer than a general-purpose CFF reader - just enough to execute one
/// program.
const Buf = struct {
    data: []const u8,
    cursor: u32 = 0,

    pub const empty: Buf = .{ .data = &.{} };

    pub fn init(data: []const u8) Buf {
        return .{ .data = data };
    }

    pub inline fn size(b: Buf) u32 {
        return @intCast(b.data.len);
    }

    pub fn skip(b: *Buf, o: u32) void {
        b.seek(b.cursor + o);
    }

    pub fn seek(b: *Buf, o: u32) void {
        b.cursor = @min(o, b.size());
    }

    pub fn peek8(b: *Buf) u8 {
        if (b.cursor >= b.size()) return 0;
        return b.data[b.cursor];
    }

    pub fn get8(b: *Buf) u8 {
        if (b.cursor >= b.size()) return 0;
        defer b.cursor += 1;
        return b.data[b.cursor];
    }

    pub fn get16(b: *Buf) u16 {
        return @truncate(b.get(2));
    }

    pub fn get32(b: *Buf) u32 {
        return b.get(4);
    }

    pub fn get(b: *Buf, n: u32) u32 {
        var v: u32 = 0;
        assert(n >= 1 and n <= 4);
        for (0..n) |_|
            v = (v << 8) | b.get8();
        return v;
    }

    /// Decodes one Type 2 charstring number operand at the cursor.
    pub fn cffInt(b: *Buf) u32 {
        const b0: i32 = b.get8();
        return switch (b0) {
            32...246 => @bitCast(b0 - 139),
            247...250 => @bitCast((b0 - 247) * 256 + b.get8() + 108),
            251...254 => @bitCast(-(b0 - 251) * 256 - b.get8() - 108),
            28 => b.get16(),
            29 => b.get32(),
            else => @panic("invalid instruction"),
        };
    }
};

/// Subroutine numbers are biased by a count-dependent amount before
/// indexing (Type 2 charstring spec, "Subroutine Operators").
fn getBiasedSubr(idx: Index, n: i32) []const u8 {
    const count = idx.count();
    var un: u32 = @bitCast(n);
    un +%= if (count >= 33900)
        32768
    else if (count >= 1240)
        1131
    else
        107;
    if (un >= count) return &.{};
    return idx.get(un);
}

fn runCharstring(cff: CffTable, glyph: GlyphIndex, ctx: *CharstringCtx) !void {
    var maskbits: u32 = 0;
    var in_header = true;
    var has_subrs = false;
    var clear_stack = false;
    var s: [48]f32 = @splat(0); // stack
    var sp: u32 = 0; // stack pointer
    var subr_buf: [10]Buf = undefined;
    var subr_stack: std.ArrayList(Buf) = .initBuffer(&subr_buf);
    var local_subrs: ?Index = null;
    // this currently ignores the initial width value, which isn't needed if we have hmtx
    var b: Buf = .init(cff.getCharString(glyph));

    while (b.cursor < b.size()) {
        var i: u32 = 0;
        clear_stack = true;
        const b0: u16 = b.get8();

        sw: switch (b0) {
            // @TODO implement hinting
            Instruction.hintmask.toInt(), // 0x13
            Instruction.cntrmask.toInt(), // 0x14
            => {
                if (in_header) maskbits += (sp / 2); // implicit "vstem"
                in_header = false;
                b.skip((maskbits + 7) / 8);
            },
            Instruction.hstem.toInt(), // 0x01
            Instruction.vstem.toInt(), // 0x03
            Instruction.hstemhm.toInt(), // 0x12
            Instruction.vstemhm.toInt(), // 0x17
            => {
                maskbits += (sp / 2);
            },
            Instruction.rmoveto.toInt() => { // 0x15
                in_header = false;
                if (sp < 2) return error.RMoveToStack;
                try ctx.rmoveTo(s[sp - 2], s[sp - 1]);
            },
            Instruction.vmoveto.toInt() => { // 0x04
                in_header = false;
                if (sp < 1) return error.VMoveToStack;
                try ctx.rmoveTo(0, s[sp - 1]);
            },
            Instruction.hmoveto.toInt() => { // 0x16
                in_header = false;
                if (sp < 1) return error.HMoveToStack;
                try ctx.rmoveTo(s[sp - 1], 0);
            },
            Instruction.rlineto.toInt() => { // 0x05
                if (sp < 2) return error.RLineToStack;
                while (i + 1 < sp) : (i += 2)
                    try ctx.rlineTo(s[i], s[i + 1]);
            },
            // hlineto/vlineto and vhcurveto/hvcurveto alternate horizontal and vertical
            // starting from a different place.
            Instruction.vlineto.toInt() => { // 0x07
                if (sp < 1) return error.VLineToStack;
                while (true) {
                    if (i >= sp) break;
                    try ctx.rlineTo(0, s[i]);
                    i += 1;
                    if (i >= sp) break;
                    try ctx.rlineTo(s[i], 0);
                    i += 1;
                }
            },
            Instruction.hlineto.toInt() => { // 0x06
                if (sp < 1) return error.HLineToStack;
                while (true) {
                    if (i >= sp) break;
                    try ctx.rlineTo(s[i], 0);
                    i += 1;
                    if (i >= sp) break;
                    try ctx.rlineTo(0, s[i]);
                    i += 1;
                }
            },
            Instruction.hvcurveto.toInt() => { // 0x1F
                if (sp < 4) return error.HCurveToStack;
                while (true) {
                    if (i + 3 >= sp) break;
                    try ctx.rccurveTo(s[i], 0, s[i + 1], s[i + 2], if (sp - i == 5) s[i + 4] else 0.0, s[i + 3]);
                    i += 4;
                    if (i + 3 >= sp) break;
                    try ctx.rccurveTo(0, s[i], s[i + 1], s[i + 2], s[i + 3], if (sp - i == 5) s[i + 4] else 0.0);
                    i += 4;
                }
            },
            Instruction.vhcurveto.toInt() => { // 0x1E
                if (sp < 4) return error.HCurveToStack;
                while (true) {
                    if (i + 3 >= sp) break;
                    try ctx.rccurveTo(0, s[i], s[i + 1], s[i + 2], s[i + 3], if (sp - i == 5) s[i + 4] else 0.0);
                    i += 4;
                    if (i + 3 >= sp) break;
                    try ctx.rccurveTo(s[i], 0, s[i + 1], s[i + 2], if (sp - i == 5) s[i + 4] else 0.0, s[i + 3]);
                    i += 4;
                }
            },
            Instruction.rrcurveto.toInt() => { // 0x08
                if (sp < 6) return error.RCurveToStack;
                while (i + 5 < sp) : (i += 6)
                    try ctx.rccurveTo(s[i], s[i + 1], s[i + 2], s[i + 3], s[i + 4], s[i + 5]);
            },
            Instruction.rcurveline.toInt() => { // 0x18
                if (sp < 8) return error.RCurveLineStack;
                while (i + 5 < sp - 2) : (i += 6)
                    try ctx.rccurveTo(s[i], s[i + 1], s[i + 2], s[i + 3], s[i + 4], s[i + 5]);
                if (i + 1 >= sp) return error.CurveLineStack;
                try ctx.rlineTo(s[i], s[i + 1]);
            },
            Instruction.rlinecurve.toInt() => { // 0x19
                if (sp < 8) return error.RLineCurveStack;
                while (i + 1 < sp - 6) : (i += 2)
                    try ctx.rlineTo(s[i], s[i + 1]);
                if (i + 5 >= sp) return error.RLineCurveStack;
                try ctx.rccurveTo(s[i], s[i + 1], s[i + 2], s[i + 3], s[i + 4], s[i + 5]);
            },
            Instruction.vvcurveto.toInt(), // 0x1A
            Instruction.hhcurveto.toInt(), // 0x1B
            => {
                if (sp < 4) return error.CurveToStack;
                var f: f32 = 0.0;
                if (sp & 1 != 0) {
                    f = s[i];
                    i += 1;
                }
                while (i + 3 < sp) : (i += 4) {
                    if (b0 == Instruction.hhcurveto.toInt()) //  0x1B
                        try ctx.rccurveTo(s[i], f, s[i + 1], s[i + 2], s[i + 3], 0.0)
                    else
                        try ctx.rccurveTo(f, s[i], s[i + 1], s[i + 2], 0.0, s[i + 3]);
                    f = 0.0;
                }
            },
            Instruction.callsubr.toInt() => { // 0x0A
                if (!has_subrs) {
                    local_subrs = cff.getLocalSubrIndexForGlyph(glyph);
                    has_subrs = true;
                }
                continue :sw Instruction.callgsubr.toInt();
                // FALLTHROUGH
            },
            Instruction.callgsubr.toInt() => { // 0x1D
                sp = std.math.sub(u32, sp, 1) catch return error.CallGSubRStack;
                const v: i32 = @trunc(s[sp]);
                subr_stack.appendBounded(b) catch return error.RecursionLimit;
                const subr_bytes = if (b0 == Instruction.callsubr.toInt()) // 0x0A
                    getBiasedSubr(local_subrs orelse return error.SubRNotFound, v)
                else
                    getBiasedSubr(cff.global_subrs_index, v);
                if (subr_bytes.len == 0) return error.SubRNotFound;
                b = .init(subr_bytes);
                clear_stack = false;
            },
            Instruction.ret.toInt() => { // 0x0B
                b = subr_stack.pop() orelse return error.ReturnOutsideSubR;
                clear_stack = false;
            },
            Instruction.endchar.toInt() => { // 0x0E
                try ctx.closeShape();
                return;
            },
            Instruction.twoByteEscape.toInt() => { // 0x0C
                const b1 = b.get8();
                switch (b1) {
                    // @TODO These "flex" implementations ignore the flex-depth and resolution,
                    // and always draw beziers.
                    Instruction.hflex.toInt() => { // 0x22
                        if (sp < 7) return error.HFlexStack;
                        const dx1 = s[0];
                        const dx2 = s[1];
                        const dy2 = s[2];
                        const dx3 = s[3];
                        const dx4 = s[4];
                        const dx5 = s[5];
                        const dx6 = s[6];
                        try ctx.rccurveTo(dx1, 0, dx2, dy2, dx3, 0);
                        try ctx.rccurveTo(dx4, 0, dx5, -dy2, dx6, 0);
                    },
                    Instruction.flex.toInt() => { // 0x23
                        if (sp < 13) return error.FlexStack;
                        const dx1 = s[0];
                        const dy1 = s[1];
                        const dx2 = s[2];
                        const dy2 = s[3];
                        const dx3 = s[4];
                        const dy3 = s[5];
                        const dx4 = s[6];
                        const dy4 = s[7];
                        const dx5 = s[8];
                        const dy5 = s[9];
                        const dx6 = s[10];
                        const dy6 = s[11];
                        //fd is s[12]
                        try ctx.rccurveTo(dx1, dy1, dx2, dy2, dx3, dy3);
                        try ctx.rccurveTo(dx4, dy4, dx5, dy5, dx6, dy6);
                    },
                    Instruction.hflex1.toInt() => { // 0x24
                        if (sp < 9) return error.HFlex1Stack;
                        const dx1 = s[0];
                        const dy1 = s[1];
                        const dx2 = s[2];
                        const dy2 = s[3];
                        const dx3 = s[4];
                        const dx4 = s[5];
                        const dx5 = s[6];
                        const dy5 = s[7];
                        const dx6 = s[8];
                        try ctx.rccurveTo(dx1, dy1, dx2, dy2, dx3, 0);
                        try ctx.rccurveTo(dx4, 0, dx5, dy5, dx6, -(dy1 + dy2 + dy5));
                    },
                    Instruction.flex1.toInt() => { // 0x25
                        if (sp < 11) return error.Flex1Stack;
                        const dx1 = s[0];
                        const dy1 = s[1];
                        const dx2 = s[2];
                        const dy2 = s[3];
                        const dx3 = s[4];
                        const dy3 = s[5];
                        const dx4 = s[6];
                        const dy4 = s[7];
                        const dx5 = s[8];
                        const dy5 = s[9];
                        var dx6 = s[10];
                        var dy6 = s[10];
                        const dx = dx1 + dx2 + dx3 + dx4 + dx5;
                        const dy = dy1 + dy2 + dy3 + dy4 + dy5;
                        if (@abs(dx) > @abs(dy))
                            dy6 = -dy
                        else
                            dx6 = -dx;
                        try ctx.rccurveTo(dx1, dy1, dx2, dy2, dx3, dy3);
                        try ctx.rccurveTo(dx4, dy4, dx5, dy5, dx6, dy6);
                    },

                    else => return error.Unimplemented,
                }
            },
            else => {
                if (b0 != 255 and b0 != 28 and b0 < 32)
                    return error.ReservedOperator;

                // push immediate
                const f: f32 = if (b0 == 255)
                    @floatFromInt(@as(i32, @intCast(b.get32() / 0x10000)))
                else blk: {
                    b.cursor -= 1;
                    break :blk @floatFromInt(@as(i16, @truncate(@as(i32, @bitCast(b.cffInt())))));
                };
                if (sp >= 48) return error.PushStackOverflow;
                s[sp] = f;
                sp += 1;
                clear_stack = false;
            },
        }
        if (clear_stack) sp = 0;
    }
    return error.NoEndChar;
}
