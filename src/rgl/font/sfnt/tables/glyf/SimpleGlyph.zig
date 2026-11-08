const bytes = @import("../../../../../bytes.zig");
const zigstd = @import("../../../../../zigstd.zig");
const GlyfTable = @import("GlyfTable.zig");

const Allocator = zigstd.Allocator;

const Self = @This();

header: GlyfTable.GlyphHeader,
/// point index of the last point of each contour.
end_pts_of_contours: []u16,
instructions: []const u8,
/// one per point (indices 0..numPoints()).
flags: []Flag,
/// absolute coordinates, already delta-accumulated from the stream.
x_coordinates: []i16,
y_coordinates: []i16,

pub const Flag = packed struct(u8) {
    on_curve_point: bool,
    x_short_vector: bool,
    y_short_vector: bool,
    repeat_flag: bool,
    x_is_same_or_positive: bool,
    y_is_same_or_positive: bool,
    overlap_simple: bool,
    reserved: bool,
};

pub fn numPoints(self: *const Self) usize {
    const n = self.end_pts_of_contours.len;
    if (n == 0) return 0;
    return @as(usize, self.end_pts_of_contours[n - 1]) + 1;
}

pub fn deinit(self: *Self, alloc: Allocator) void {
    alloc.free(self.y_coordinates);
    alloc.free(self.x_coordinates);
    alloc.free(self.flags);
    alloc.free(self.instructions);
    alloc.free(self.end_pts_of_contours);
}

pub fn parse(alloc: Allocator, glyph_data: []const u8, header: GlyfTable.GlyphHeader, check: bool) !Self {
    const num_contours: usize = @intCast(header.number_of_contours);
    const header_end = @sizeOf(GlyfTable.GlyphHeader);

    if (check) try bytes.checkArray(glyph_data, header_end, num_contours, u16);
    const end_pts_of_contours = try bytes.allocArrayBEToNative(alloc, u16, num_contours, glyph_data[header_end..]);
    errdefer alloc.free(end_pts_of_contours);

    const num_points: usize = if (num_contours == 0) 0 else @as(usize, end_pts_of_contours[num_contours - 1]) + 1;

    const instr_len_pos = header_end + num_contours * @sizeOf(u16);
    if (check) try bytes.checkSize(glyph_data, instr_len_pos, u16);
    const instruction_length = bytes.readIntBE(u16, glyph_data[instr_len_pos..]);

    const instructions_start = instr_len_pos + @sizeOf(u16);
    if (check) try bytes.checkSlice(glyph_data, instructions_start, instruction_length);
    const instructions = try alloc.dupe(u8, glyph_data[instructions_start..][0..instruction_length]);
    errdefer alloc.free(instructions);

    const flags = try alloc.alloc(Flag, num_points);
    errdefer alloc.free(flags);

    var pos = instructions_start + instruction_length;
    var i: usize = 0;
    while (i < num_points) {
        if (check) try bytes.checkSize(glyph_data, pos, u8);
        const flag: Flag = @bitCast(glyph_data[pos]);
        pos += 1;
        flags[i] = flag;
        i += 1;

        if (flag.repeat_flag) {
            if (check) try bytes.checkSize(glyph_data, pos, u8);
            const repeat_count = glyph_data[pos];
            pos += 1;
            const n = @min(repeat_count, num_points - i);
            @memset(flags[i..][0..n], flag);
            i += n;
        }
    }

    const x_coordinates = try readCoordinates(alloc, glyph_data, &pos, flags, .x, check);
    errdefer alloc.free(x_coordinates);
    const y_coordinates = try readCoordinates(alloc, glyph_data, &pos, flags, .y, check);
    errdefer alloc.free(y_coordinates);

    return .{
        .header = header,
        .end_pts_of_contours = end_pts_of_contours,
        .instructions = instructions,
        .flags = flags,
        .x_coordinates = x_coordinates,
        .y_coordinates = y_coordinates,
    };
}

const Axis = enum { x, y };

/// One axis's coordinates from `pos` on, delta-accumulated to absolute values.
fn readCoordinates(alloc: Allocator, glyph_data: []const u8, pos: *usize, flags: []const Flag, comptime axis: Axis, check: bool) ![]i16 {
    const coordinates = try alloc.alloc(i16, flags.len);
    errdefer alloc.free(coordinates);
    var v: i32 = 0;
    for (flags, 0..) |flag, idx| {
        const short_vector, const same_or_positive = switch (axis) {
            .x => .{ flag.x_short_vector, flag.x_is_same_or_positive },
            .y => .{ flag.y_short_vector, flag.y_is_same_or_positive },
        };
        if (short_vector) {
            if (check) try bytes.checkSize(glyph_data, pos.*, u8);
            const delta: i32 = glyph_data[pos.*];
            pos.* += 1;
            v += if (same_or_positive) delta else -delta;
        } else if (!same_or_positive) {
            if (check) try bytes.checkSize(glyph_data, pos.*, i16);
            v += bytes.readIntBE(i16, glyph_data[pos.*..]);
            pos.* += 2;
        }
        // else: short vector clear and same-or-positive set means "same as
        // previous", a zero delta with no bytes consumed.
        coordinates[idx] = @intCast(v);
    }
    return coordinates;
}
