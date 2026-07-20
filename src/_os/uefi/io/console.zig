//! Write bytes to a UEFI `SimpleTextOutput` console: UCS-2 encoding,
//! "\n" -> "\r\n" (LF doesn't return column). Assumes ASCII -- no UTF-8
//! decoding, multi-byte chars break. Shared by `stderr.zig`'s `debug_io`
//! drain and `operate.zig`'s real stdout/stderr streaming writes -- both
//! ultimately write to the same firmware console.

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
