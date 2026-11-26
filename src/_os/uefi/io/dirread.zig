//! `dirRead` over the boot volume's root (`state.root_dir`, the only `Io.Dir`
//! here; see `dir.zig`): one `File.Info.File` per `File.read`, names from
//! UCS-2 to WTF-8 in the reader's buffer, `.` and `..` skipped.

const std = @import("std");
const uefi = std.os.uefi;
const Io = std.Io;

const state = @import("state.zig");

/// The longest name, as WTF-8 from at most 255 UCS-2 units.
pub const max_name_bytes = 255 * 3;

const FileInfo = uefi.protocol.File.Info.File;
/// One entry's `File.Info.File` with a name of up to 255 units and its NUL.
var info_buf: [@sizeOf(FileInfo) + 256 * 2]u8 align(@alignOf(FileInfo)) = undefined;

pub fn dirRead(userdata: ?*anyopaque, r: *Io.Dir.Reader, entries: []Io.Dir.Entry) Io.Dir.Reader.Error!usize {
    _ = userdata;
    const root = state.root_dir orelse return error.Unexpected; // openRootDir was never called
    switch (r.state) {
        .finished => return 0,
        .reset => {
            root.setPosition(0) catch return error.Unexpected;
            r.state = .reading;
        },
        .reading => {},
    }
    // names live in `r.buffer` until the next call
    r.end = 0;
    var n: usize = 0;
    while (n < entries.len and r.buffer.len - r.end >= max_name_bytes) {
        const len = root.read(&info_buf) catch |err| return switch (err) {
            error.BufferTooSmall => error.SystemResources,
            else => error.Unexpected,
        };
        if (len == 0) {
            r.state = .finished;
            break;
        }
        const info: *const FileInfo = @ptrCast(&info_buf);
        const name = std.mem.span(info.getFileName());
        if (isDot(name)) continue;
        const out = r.buffer[r.end..];
        const name_len = std.unicode.wtf16LeToWtf8(out, name);
        entries[n] = .{
            .name = out[0..name_len],
            .kind = if (info.attribute.directory) .directory else .file,
            .inode = if (Io.File.INode == void) {} else 0,
        };
        r.end += name_len;
        n += 1;
    }
    return n;
}

fn isDot(name: []const u16) bool {
    return (name.len == 1 and name[0] == '.') or (name.len == 2 and name[0] == '.' and name[1] == '.');
}
