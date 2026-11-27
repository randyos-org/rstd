//! Tests: `devpath.covers`, which decides whether stderr also goes to
//! `std_err` (#740).

const std = @import("std");
const devpath = @import("../devpath.zig");

const testing = std.testing;

/// `parts`, one after another; comptime only.
fn cat(comptime parts: []const []const u8) []const u8 {
    var len: usize = 0;
    for (parts) |p| len += p.len;
    var out: [len]u8 = undefined;
    var at: usize = 0;
    for (parts) |p| {
        @memcpy(out[at..][0..p.len], p);
        at += p.len;
    }
    const final = out;
    return &final;
}

/// A device path node: `kind`, `sub`, then `body`.
fn node(comptime kind: u8, comptime sub: u8, comptime body: []const u8) []const u8 {
    const len = 4 + body.len;
    return cat(&.{ &.{ kind, sub, len & 0xFF, len >> 8 }, body });
}

const end_instance = node(0x7F, 0x01, "");
const end_entire = node(0x7F, 0xFF, "");
// a PCI node and a UART node, as a serial console's path has them
const serial: []const u8 = cat(&.{ node(0x01, 0x01, &.{ 0, 0x1F }), node(0x03, 0x0E, "uart") });
const video: []const u8 = cat(&.{ node(0x01, 0x01, &.{ 0, 0x02 }), node(0x02, 0x01, "acpi") });
const serial_only: []const u8 = cat(&.{ serial, end_entire });
const video_only: []const u8 = cat(&.{ video, end_entire });
const both: []const u8 = cat(&.{ serial, end_instance, video, end_entire });

test "UEFI device paths: the same devices, in any order, cover each other" {
    const swapped: []const u8 = comptime cat(&.{ video, end_instance, serial, end_entire });
    try testing.expect(devpath.covers(both, both));
    try testing.expect(devpath.covers(both, swapped));
    try testing.expect(devpath.covers(swapped, both));
    try testing.expect(devpath.covers(serial_only, serial_only));
}

test "UEFI device paths: a superset covers, a subset or other devices do not" {
    try testing.expect(devpath.covers(both, serial_only));
    try testing.expect(devpath.covers(both, video_only));
    try testing.expect(!devpath.covers(serial_only, both));
    try testing.expect(!devpath.covers(serial_only, video_only));
}

test "UEFI device paths: empty or malformed lists never cover" {
    try testing.expect(!devpath.covers("", ""));
    try testing.expect(!devpath.covers(end_entire, end_entire));
    try testing.expect(!devpath.covers(both, end_entire));
    // no end-entire node
    try testing.expect(!devpath.covers(serial, serial));
    // a node longer than the list, after a good instance
    const long: []const u8 = comptime cat(&.{ serial, end_instance, &.{ 0x01, 0x01, 0x40, 0x00 } });
    try testing.expect(!devpath.covers(long, serial_only));
    try testing.expect(!devpath.covers(both, long));
    // a node shorter than its header
    const short: []const u8 = comptime cat(&.{ &.{ 0x01, 0x01, 0x02, 0x00 }, end_entire });
    try testing.expect(!devpath.covers(short, short));
}
