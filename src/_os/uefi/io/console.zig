//! Bytes to a UEFI `SimpleTextOutput` console, as UCS-2 with "\n" as "\r\n"; ASCII only. Shared by `stderr.zig`'s drain and
//! `operate.zig`'s stdout and stderr writes, which land on the same console.

const std = @import("std");
const uefi = std.os.uefi;

pub fn writeConsole(protocol: ?*uefi.protocol.SimpleTextOutput, bytes: []const u8) void {
    const out = protocol orelse return;
    // batched: outputString is a firmware call, one per char would be slow
    var buf: [64:0]u16 = undefined;
    var len: usize = 0;
    for (bytes) |byte| {
        if (len + 2 > buf.len) {
            buf[len] = 0;
            _ = out.outputString(buf[0..len :0]) catch false;
            len = 0;
        }
        if (byte == '\n') {
            buf[len] = '\r';
            len += 1;
        }
        buf[len] = byte;
        len += 1;
    }
    if (len > 0) {
        buf[len] = 0;
        _ = out.outputString(buf[0..len :0]) catch false;
    }
}

/// `header`, then `data` with its last slice `splat` times; the bytes written.
pub fn writeSplat(protocol: ?*uefi.protocol.SimpleTextOutput, header: []const u8, data: []const []const u8, splat: usize) usize {
    var n: usize = header.len;
    writeConsole(protocol, header);
    if (data.len != 0) {
        for (data[0 .. data.len - 1]) |chunk| {
            writeConsole(protocol, chunk);
            n += chunk.len;
        }
        const last = data[data.len - 1];
        for (0..splat) |_| writeConsole(protocol, last);
        n += last.len * splat;
    }
    return n;
}
