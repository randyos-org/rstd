//! The OpenType spec's worked `fvar` example: two axes (wght, wdth), four
//! named instances, 112 bytes. Shared by the fvar and coordinate tests.

const bytes = @import("../../../../../bytes.zig");
const zigstd = @import("../../../../../zigstd.zig");

pub fn build(alloc: zigstd.Allocator) ![]u8 {
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
