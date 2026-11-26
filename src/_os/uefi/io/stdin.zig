//! fd 0 from the UEFI `SimpleTextInput` console (`state.con_in`): keystrokes
//! as UTF-8, Enter as "\n", with no echo. Keys with no character (arrows,
//! function keys) are dropped.

const std = @import("std");
const Io = std.Io;

const state = @import("state.zig");

const Error = Io.Operation.FileReadStreaming.Error;

/// A key's UTF-8 bytes that did not fit the caller's buffers.
var pending: [3]u8 = undefined;
var pending_len: usize = 0;

/// Waits for a key, then reads every key already typed into `data`; the
/// bytes read. `EndOfStream` with no console (none, or after `state.stop`).
pub fn read(data: []const []u8) Error!usize {
    const in = state.con_in orelse return error.EndOfStream;
    const bs = state.bootServices() orelse return error.EndOfStream;
    var total: usize = 0;
    for (data) |buffer| {
        var at: usize = 0;
        while (at < buffer.len) {
            if (pending_len == 0) {
                const blocking = total + at == 0;
                if (!try nextKey(bs, in, blocking)) break;
            }
            const n = @min(pending_len, buffer.len - at);
            @memcpy(buffer[at..][0..n], pending[0..n]);
            @memmove(pending[0 .. pending_len - n], pending[n..pending_len]);
            pending_len -= n;
            at += n;
        }
        total += at;
        if (at < buffer.len) break;
    }
    return total;
}

/// Encodes the next keystroke into `pending`; false when none is waiting and
/// `blocking` is false.
fn nextKey(bs: *std.os.uefi.tables.BootServices, in: *std.os.uefi.protocol.SimpleTextInput, blocking: bool) Error!bool {
    while (true) {
        const key = in.readKeyStroke() catch |err| switch (err) {
            error.NotReady => {
                if (!blocking) return false;
                _ = bs.waitForEvent(&[_]std.os.uefi.Event{in.wait_for_key}) catch return error.InputOutput;
                continue;
            },
            else => return error.InputOutput,
        };
        const c: u21 = switch (key.unicode_char) {
            0 => continue,
            '\r' => '\n',
            else => |u| u,
        };
        pending_len = std.unicode.utf8Encode(c, &pending) catch continue;
        return true;
    }
}
