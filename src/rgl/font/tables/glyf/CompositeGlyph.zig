const std = @import("std");
const types = @import("../../stream_types.zig");
const bytes = @import("../../../../bytes.zig");
const zigstd = @import("../../../../zigstd.zig");
const GlyfTable = @import("GlyfTable.zig");

const Allocator = zigstd.Allocator;
const GlyphIndex = types.GlyphIndex;
const F2Dot14 = types.F2Dot14;

pub const Component = struct {
    glyph_index: GlyphIndex,
    /// Offset added to the component's points, in the parent's design
    /// units (once scaled_component_offset/unscaled_component_offset is
    /// accounted for by the caller).
    dx: f32,
    dy: f32,
    /// x' = xscale*x + scale10*y; y' = scale01*x + yscale*y
    xscale: f32,
    scale01: f32,
    scale10: f32,
    yscale: f32,

    round_xy_to_grid: bool,
    use_my_metrics: bool,
    overlap_compound: bool,
    scaled_component_offset: bool,
    unscaled_component_offset: bool,
};

const Self = @This();

header: GlyfTable.GlyphHeader,
components: []Component,
/// Instructions for the composite as a whole (separate from any
/// instructions belonging to individual component glyphs).
instructions: []const u8,

pub fn destroy(self: *Self, alloc: Allocator) void {
    alloc.free(self.instructions);
    alloc.free(self.components);
}

const ComponentFlags = packed struct(u16) {
    args_are_words: bool,
    args_are_xy_values: bool,
    round_xy_to_grid: bool,
    we_have_a_scale: bool,
    _reserved0: bool = false,
    more_components: bool,
    we_have_an_x_and_y_scale: bool,
    we_have_a_two_by_two: bool,
    we_have_instructions: bool,
    use_my_metrics: bool,
    overlap_compound: bool,
    scaled_component_offset: bool,
    unscaled_component_offset: bool,
    _reserved1: u3 = 0,
};

fn readF2Dot14(glyph_data: []const u8, pos: usize) f32 {
    return F2Dot14.fromBytes(glyph_data[pos..]).toFloat(f32);
}

pub fn parse(alloc: Allocator, glyph_data: []const u8, header: GlyfTable.GlyphHeader, check: bool) !Self {
    var components: std.ArrayList(Component) = .empty;
    errdefer components.deinit(alloc);

    var pos: usize = @sizeOf(GlyfTable.GlyphHeader);
    var last_flags: ComponentFlags = undefined;
    while (true) {
        if (check) try bytes.checkSize(glyph_data, pos, u16);
        const flags: ComponentFlags = @bitCast(bytes.readIntBE(u16, glyph_data[pos..]));
        pos += 2;

        if (check) try bytes.checkSize(glyph_data, pos, u16);
        const glyph_index: GlyphIndex = @fromBackingInt(bytes.readIntBE(u16, glyph_data[pos..]));
        pos += 2;

        // Point-matching (argument1/2 as point numbers to align) isn't
        // supported; every composite's first component must use xy values
        // per spec, so this only affects rare multi-component glyphs.
        if (!flags.args_are_xy_values) return GlyfTable.Error.PointMatchingNotSupported;

        var dx: f32 = undefined;
        var dy: f32 = undefined;
        if (flags.args_are_words) {
            if (check) try bytes.checkSize(glyph_data, pos, i16);
            dx = @floatFromInt(bytes.readIntBE(i16, glyph_data[pos..]));
            pos += 2;
            if (check) try bytes.checkSize(glyph_data, pos, i16);
            dy = @floatFromInt(bytes.readIntBE(i16, glyph_data[pos..]));
            pos += 2;
        } else {
            if (check) try bytes.checkSize(glyph_data, pos, i8);
            dx = @floatFromInt(@as(i8, @bitCast(glyph_data[pos])));
            pos += 1;
            if (check) try bytes.checkSize(glyph_data, pos, i8);
            dy = @floatFromInt(@as(i8, @bitCast(glyph_data[pos])));
            pos += 1;
        }

        var xscale: f32 = 1;
        var scale01: f32 = 0;
        var scale10: f32 = 0;
        var yscale: f32 = 1;

        if (flags.we_have_a_scale) {
            if (check) try bytes.checkSize(glyph_data, pos, F2Dot14);
            const s = readF2Dot14(glyph_data, pos);
            pos += 2;
            xscale = s;
            yscale = s;
        } else if (flags.we_have_an_x_and_y_scale) {
            if (check) try bytes.checkArray(glyph_data, pos, 2, F2Dot14);
            xscale = readF2Dot14(glyph_data, pos);
            pos += 2;
            yscale = readF2Dot14(glyph_data, pos);
            pos += 2;
        } else if (flags.we_have_a_two_by_two) {
            if (check) try bytes.checkArray(glyph_data, pos, 4, F2Dot14);
            xscale = readF2Dot14(glyph_data, pos);
            pos += 2;
            scale01 = readF2Dot14(glyph_data, pos);
            pos += 2;
            scale10 = readF2Dot14(glyph_data, pos);
            pos += 2;
            yscale = readF2Dot14(glyph_data, pos);
            pos += 2;
        }

        try components.append(alloc, .{
            .glyph_index = glyph_index,
            .dx = dx,
            .dy = dy,
            .xscale = xscale,
            .scale01 = scale01,
            .scale10 = scale10,
            .yscale = yscale,
            .round_xy_to_grid = flags.round_xy_to_grid,
            .use_my_metrics = flags.use_my_metrics,
            .overlap_compound = flags.overlap_compound,
            .scaled_component_offset = flags.scaled_component_offset,
            .unscaled_component_offset = flags.unscaled_component_offset,
        });

        last_flags = flags;
        if (!flags.more_components) break;
    }

    var instructions: []const u8 = &.{};
    if (last_flags.we_have_instructions) {
        if (check) try bytes.checkSize(glyph_data, pos, u16);
        const num_instr = bytes.readIntBE(u16, glyph_data[pos..]);
        pos += 2;
        if (check) try bytes.checkSlice(glyph_data, pos, num_instr);
        instructions = try alloc.dupe(u8, glyph_data[pos..][0..num_instr]);
    }
    errdefer alloc.free(instructions);

    return .{
        .header = header,
        .components = try components.toOwnedSlice(alloc),
        .instructions = instructions,
    };
}
