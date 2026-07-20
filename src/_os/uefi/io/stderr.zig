const std = @import("std");
const uefi = std.os.uefi;
const Io = std.Io;

const state = @import("state.zig");
const console = @import("console.zig");

/// `Io.Writer.VTable` backing `Threaded.stderr_writer` (see state.zig).
/// `file_writer` must point at Io.File.Writer but never touches `.file` --
/// writes go through `.interface`, this vtable. state.zig's `File` is a
/// dummy; `stderrDrain` does the real work.
pub const stderr_writer_vtable: Io.Writer.VTable = .{
    .drain = stderrDrain,
};

fn stderrDrain(w: *Io.Writer, data: []const []const u8, splat: usize) Io.Writer.Error!usize {
    const buffered = w.buffered();
    console.writeConsole(state.con_out, buffered);
    var n: usize = buffered.len;
    if (data.len != 0) {
        for (data[0 .. data.len - 1]) |chunk| {
            console.writeConsole(state.con_out, chunk);
            n += chunk.len;
        }
        const last = data[data.len - 1];
        for (0..splat) |_| console.writeConsole(state.con_out, last);
        n += last.len * splat;
    }
    return w.consume(n);
}

pub fn lockStderr(userdata: ?*anyopaque, terminal_mode: ?Io.Terminal.Mode) Io.Cancelable!Io.LockedStderr {
    _ = userdata;
    state.stderr_lock_count += 1;
    return .{
        .file_writer = &state.stderr_writer,
        // UEFI console has no ANSI escapes; `.no_color` no-ops setColor.
        // (SetAttribute color needs a custom Terminal mode, not supported)
        .terminal_mode = terminal_mode orelse .no_color,
    };
}

pub fn tryLockStderr(userdata: ?*anyopaque, terminal_mode: ?Io.Terminal.Mode) Io.Cancelable!?Io.LockedStderr {
    // single-threaded: lock always available (recursively)
    return try lockStderr(userdata, terminal_mode);
}

pub fn unlockStderr(userdata: ?*anyopaque) void {
    _ = userdata;
    // mirrors Io.Threaded.unlockStderr: lock owns final flush, clear caller buffer
    if (state.stderr_writer.err == null) state.stderr_writer.interface.flush() catch {};
    state.stderr_writer.err = null;
    state.stderr_writer.interface.end = 0;
    state.stderr_writer.interface.buffer = &.{};
    state.stderr_lock_count -= 1;
}
