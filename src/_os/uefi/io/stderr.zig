const std = @import("std");
const Io = std.Io;

const state = @import("state.zig");
const console = @import("console.zig");

/// `Io.Writer.VTable` of `state.stderr_writer`, the `Io.File.Writer` that
/// `lockStderr` hands out; its `.file` is never written, `stderrDrain` is.
pub const stderr_writer_vtable: Io.Writer.VTable = .{
    .drain = stderrDrain,
};

fn stderrDrain(w: *Io.Writer, data: []const []const u8, splat: usize) Io.Writer.Error!usize {
    return w.consume(console.writeSplat(state.con_out, w.buffered(), data, splat));
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
