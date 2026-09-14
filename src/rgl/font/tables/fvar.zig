const types = @import("../stream_types.zig");
const bytes = @import("../../../bytes.zig");
const zigstd = @import("../../../zigstd.zig");

/// Schema for the fvar header, used only to compute field offsets/sizes at
/// comptime; never instantiated.
const StreamLayout = extern struct {
    major_version: u16,
    minor_version: u16,
    axes_array_offset: types.Offset16,
    reserved: u16,
    axis_count: u16,
    axis_size: u16,
    instance_count: u16,
    instance_size: u16,
};

/// One VariationAxisRecord. Callers should index the axes array using
/// axisSize() (see FvarTable.getAxis()), not @sizeOf, so a future
/// minor-version-extended record still walks at the right stride.
pub const VariationAxisRecord = extern struct {
    axis_tag: types.TagStr,
    min_value: types.Fixed,
    default_value: types.Fixed,
    max_value: types.Fixed,
    flags: u16,
    axis_name_id: u16,
};

pub const AxisFlags = packed struct(u16) {
    /// Font developer's recommendation that this axis not be exposed
    /// directly in end-user interfaces.
    hidden_axis: bool,
    _reserved: u15 = 0,
};

/// Font variations: the axes of variation a variable font supports, and any
/// designer-provided named instances (e.g. "Bold"). Resolving a coordinate
/// set into actual outline/metric deltas is gvar/HVAR/MVAR's job - fvar
/// only describes the variation space itself.
pub const FvarTable = struct {
    const Self = @This();

    table_data: []const u8,

    pub fn init(table_data: []const u8, check: bool) !Self {
        if (check) try bytes.checkSize(table_data, 0, StreamLayout);
        const self: Self = .{ .table_data = table_data };

        if (check) {
            try bytes.checkSlice(table_data, self.axesArrayOffset(), @as(usize, self.axisCount()) * self.axisSize());
            try bytes.checkSlice(table_data, self.instancesStart(), @as(usize, self.instanceCount()) * self.instanceSize());
        }

        return self;
    }

    pub inline fn majorVersion(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "major_version")..]);
    }

    pub inline fn minorVersion(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "minor_version")..]);
    }

    pub inline fn axesArrayOffset(self: Self) types.Offset16 {
        return bytes.readIntBE(types.Offset16, self.table_data[@offsetOf(StreamLayout, "axes_array_offset")..]);
    }

    pub inline fn axisCount(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "axis_count")..]);
    }

    /// Size in bytes of each VariationAxisRecord in this font; use this
    /// (not @sizeOf(VariationAxisRecord)) as the stride when walking the
    /// axes array, per spec.
    pub inline fn axisSize(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "axis_size")..]);
    }

    pub inline fn instanceCount(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "instance_count")..]);
    }

    /// Size in bytes of each InstanceRecord in this font; either
    /// axisCount*sizeof(Fixed)+4 (no postScriptNameID) or +6 (with it) -
    /// use this (not a computed size) as the stride when walking the
    /// instances array, per spec.
    pub inline fn instanceSize(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "instance_size")..]);
    }

    pub fn getAxis(self: Self, index: usize) VariationAxisRecord {
        const off = self.axesArrayOffset() + index * self.axisSize();
        return bytes.typeBEToNative(VariationAxisRecord, self.table_data[off..]);
    }

    inline fn instancesStart(self: Self) usize {
        return self.axesArrayOffset() + @as(usize, self.axisCount()) * self.axisSize();
    }

    pub fn getInstanceSubfamilyNameId(self: Self, index: usize) u16 {
        return bytes.getU16BE(self.table_data[self.instancesStart() + index * self.instanceSize() ..]);
    }

    pub fn getInstanceFlags(self: Self, index: usize) u16 {
        const off = self.instancesStart() + index * self.instanceSize() + @sizeOf(u16);
        return bytes.getU16BE(self.table_data[off..]);
    }

    /// coordinates[axis], one per axis defined in this font (0..axisCount()).
    pub fn getInstanceCoordinate(self: Self, index: usize, axis: usize) types.Fixed {
        const off = self.instancesStart() + index * self.instanceSize() + 2 * @sizeOf(u16) + axis * @sizeOf(types.Fixed);
        return types.Fixed.fromBytes(self.table_data[off..]);
    }

    inline fn coordsEnd(self: Self) usize {
        return 2 * @sizeOf(u16) + @as(usize, self.axisCount()) * @sizeOf(types.Fixed);
    }

    /// Optional per spec: present only if instanceSize is large enough to
    /// include it. Null if this font's InstanceRecords omit it.
    pub fn getInstancePostScriptNameId(self: Self, index: usize) ?u16 {
        if (self.instanceSize() < self.coordsEnd() + @sizeOf(u16)) return null;
        const off = self.instancesStart() + index * self.instanceSize() + self.coordsEnd();
        return bytes.getU16BE(self.table_data[off..]);
    }
};

test "fvar StreamLayout/VariationAxisRecord match the OpenType stream layout" {
    try zigstd.testExpectEqual(@as(usize, 16), @sizeOf(StreamLayout));
    try zigstd.testExpectEqual(@as(usize, 6), @offsetOf(StreamLayout, "reserved"));
    try zigstd.testExpectEqual(@as(usize, 20), @sizeOf(VariationAxisRecord));
    try zigstd.testExpectEqual(@as(usize, 16), @offsetOf(VariationAxisRecord, "flags"));
}

// Builds the exact example from the OpenType spec's own 'fvar' worked
// example: two axes (wght, wdth), four named instances, 112 bytes total.
fn buildFvarExample(alloc: zigstd.Allocator) ![]u8 {
    const header_size = 16;
    const axis_size = 20;
    const axis_count = 2;
    const instance_size = 14; // 2 (subfamilyNameID) + 2 (flags) + 2*4 (coords) + 2 (postScriptNameID)
    const instance_count = 4;
    const axes_start = header_size;
    const instances_start = axes_start + axis_count * axis_size;
    const total = instances_start + instance_count * instance_size;

    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);

    bytes.setU16BE(buf, 1); // majorVersion
    bytes.setU16BE(buf[2..], 0); // minorVersion
    bytes.setU16BE(buf[4..], axes_start); // axesArrayOffset
    bytes.setU16BE(buf[6..], 2); // reserved
    bytes.setU16BE(buf[8..], axis_count);
    bytes.setU16BE(buf[10..], axis_size);
    bytes.setU16BE(buf[12..], instance_count);
    bytes.setU16BE(buf[14..], instance_size);

    const Axis = struct { tag: []const u8, min: i32, default: i32, max: i32, name_id: u16 };
    const axes = [_]Axis{
        .{ .tag = "wght", .min = 0x012C0000, .default = 0x01900000, .max = 0x02BC0000, .name_id = 0x0100 },
        .{ .tag = "wdth", .min = 0x003E8000, .default = 0x00640000, .max = 0x00960000, .name_id = 0x0101 },
    };
    for (axes, 0..) |a, i| {
        const off = axes_start + i * axis_size;
        @memcpy(buf[off..][0..4], a.tag);
        bytes.setI32BE(buf[off + 4 ..], a.min);
        bytes.setI32BE(buf[off + 8 ..], a.default);
        bytes.setI32BE(buf[off + 12 ..], a.max);
        bytes.setU16BE(buf[off + 16 ..], 0); // flags
        bytes.setU16BE(buf[off + 18 ..], a.name_id);
    }

    const Instance = struct { subfamily_id: u16, wght: i32, wdth: i32, ps_id: u16 };
    const instances = [_]Instance{
        .{ .subfamily_id = 0x0102, .wght = 0x01900000, .wdth = 0x00640000, .ps_id = 0x0106 },
        .{ .subfamily_id = 0x0103, .wght = 0x02BC0000, .wdth = 0x00640000, .ps_id = 0x0107 },
        .{ .subfamily_id = 0x0104, .wght = 0x01900000, .wdth = 0x004B0000, .ps_id = 0x0108 },
        .{ .subfamily_id = 0x0105, .wght = 0x02BC0000, .wdth = 0x004B0000, .ps_id = 0x0109 },
    };
    for (instances, 0..) |inst, i| {
        const off = instances_start + i * instance_size;
        bytes.setU16BE(buf[off..], inst.subfamily_id);
        bytes.setU16BE(buf[off + 2 ..], 0); // flags
        bytes.setI32BE(buf[off + 4 ..], inst.wght);
        bytes.setI32BE(buf[off + 8 ..], inst.wdth);
        bytes.setU16BE(buf[off + 12 ..], inst.ps_id);
    }

    return buf;
}

test "FvarTable reads the spec's own worked example (112 bytes, 2 axes, 4 instances)" {
    const alloc = zigstd.testAllocator;
    const buf = try buildFvarExample(alloc);
    defer alloc.free(buf);

    try zigstd.testExpectEqual(@as(usize, 112), buf.len);
    const fvar = try FvarTable.init(buf, true);

    try zigstd.testExpectEqual(@as(u16, 2), fvar.axisCount());
    try zigstd.testExpectEqual(@as(u16, 4), fvar.instanceCount());

    const wght = fvar.getAxis(0);
    try zigstd.testExpectEqualSlices(u8, "wght", &wght.axis_tag.data);
    try zigstd.testExpectEqual(@as(f64, 300), wght.min_value.toFloat(f64));
    try zigstd.testExpectEqual(@as(f64, 400), wght.default_value.toFloat(f64));
    try zigstd.testExpectEqual(@as(f64, 700), wght.max_value.toFloat(f64));
    try zigstd.testExpectEqual(@as(u16, 256), wght.axis_name_id);

    const wdth = fvar.getAxis(1);
    try zigstd.testExpectEqualSlices(u8, "wdth", &wdth.axis_tag.data);
    try zigstd.testExpectEqual(@as(f64, 62.5), wdth.min_value.toFloat(f64));

    try zigstd.testExpectEqual(@as(u16, 258), fvar.getInstanceSubfamilyNameId(0));
    try zigstd.testExpectEqual(@as(f64, 400), fvar.getInstanceCoordinate(0, 0).toFloat(f64));
    try zigstd.testExpectEqual(@as(f64, 100), fvar.getInstanceCoordinate(0, 1).toFloat(f64));
    try zigstd.testExpectEqual(@as(?u16, 262), fvar.getInstancePostScriptNameId(0));

    try zigstd.testExpectEqual(@as(u16, 261), fvar.getInstanceSubfamilyNameId(3));
    try zigstd.testExpectEqual(@as(f64, 700), fvar.getInstanceCoordinate(3, 0).toFloat(f64));
    try zigstd.testExpectEqual(@as(f64, 75), fvar.getInstanceCoordinate(3, 1).toFloat(f64));
    try zigstd.testExpectEqual(@as(?u16, 265), fvar.getInstancePostScriptNameId(3));
}

test "FvarTable.getInstancePostScriptNameId returns null when instanceSize omits it" {
    const alloc = zigstd.testAllocator;
    const header_size = 16;
    const axis_size = 20;
    const instance_size = 6; // subfamilyNameID + flags + 1 coordinate, no postScriptNameID
    const total = header_size + axis_size + instance_size;

    const buf = try alloc.alloc(u8, total);
    defer alloc.free(buf);
    @memset(buf, 0);

    bytes.setU16BE(buf, 1);
    bytes.setU16BE(buf[4..], header_size); // axesArrayOffset
    bytes.setU16BE(buf[6..], 2); // reserved
    bytes.setU16BE(buf[8..], 1); // axisCount
    bytes.setU16BE(buf[10..], axis_size);
    bytes.setU16BE(buf[12..], 1); // instanceCount
    bytes.setU16BE(buf[14..], instance_size);
    @memcpy(buf[header_size..][0..4], "wght");

    const fvar = try FvarTable.init(buf, true);
    try zigstd.testExpectEqual(@as(?u16, null), fvar.getInstancePostScriptNameId(0));
}
