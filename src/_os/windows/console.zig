//! Console setup for GUI-subsystem Windows programs, which start without a console and, unless redirected, without valid std handles,
//! which Zig reads from the PEB on every print. `setup` points them at a console or NUL first, leaving redirected ones be.
const std = @import("std");
const sysinfo = @import("builtin");
const win32 = @import("win32");
const bytes = @import("../../bytes.zig");
const zigstd = @import("../../zigstd.zig");

const L = bytes.L;
const console = win32.system.console;
const fs = win32.storage.file_system;
const kernel32 = win32.kernel32;

/// `setup`, with a visible console when the command line has `--debug` (the
/// parent's, else a new one), else the NUL device.
pub fn maybeAllocDebugConsole(init: zigstd.MainInit) bool {
    // `toSlice` documents that its result may reinterpret/point into the
    // allocator's memory in ways that aren't a single well-formed `alloc`
    // call; freeing it back piecewise (e.g. via a GeneralPurposeAllocator)
    // trips its bucket/canary accounting. It requires an arena, whose
    // backing memory `init.arena` already cleans up on process exit, so no
    // manual free is needed here.
    const args = init.minimal.args.toSlice(init.arena.allocator()) catch return false;
    if (sysinfo.target.os.tag != .windows) return false;
    return setup(debugRequested(args));
}

/// What `maybeAllocDebugConsole` does, for a caller that has decided for
/// itself whether it wants a console. With `want_console`, output goes to the
/// launching console if there is one, else to a newly allocated console
/// window; otherwise it is discarded. Returns whether output reaches a
/// console. Call before anything prints.
pub fn setup(want_console: bool) bool {
    if (sysinfo.target.os.tag != .windows) return true;
    if (redirected()) return false;
    if (want_console and setupDebugConsole()) return true;
    setStdHandlesToNul();
    return false;
}

/// Whether `args` (as from `std.process.Args.toSlice`) contain "--debug".
pub fn debugRequested(args: []const [:0]const u8) bool {
    if (args.len == 0) return false;
    for (args[1..]) |arg| {
        if (zigstd.equal(u8, arg, "--debug")) return true;
    }
    return false;
}

/// GetFileType results.
const file_type_disk: u32 = 1;
const file_type_pipe: u32 = 3;

/// Whether the launcher redirected stderr to a file or pipe, which must not
/// be replaced. Read from the PEB, where Zig's own prints read it: zigwin32's
/// GetStdHandle and std's PEB type both declare it non-optional, but a
/// GUI-subsystem process without redirection has 0 there, so the field is
/// read as a plain integer.
///
/// A console handle doesn't count: a launcher such as `zig build run` can
/// pass its own console handles, which a GUI-subsystem process (not attached
/// to that console) cannot write to. Writing to them fails, so those get the
/// same treatment as no handle at all.
fn redirected() bool {
    const raw: *const usize = @ptrCast(&std.os.windows.peb().ProcessParameters.hStdError);
    const value = raw.*;
    if (value == 0 or value == @intFromPtr(win32.foundation.INVALID_HANDLE_VALUE)) return false;
    const file_type = kernel32.GetFileType(@ptrFromInt(value));
    return file_type == file_type_disk or file_type == file_type_pipe;
}

/// Attempts to attach/allocate a real console and point
/// STD_OUTPUT_HANDLE/STD_ERROR_HANDLE at it. Returns whether it succeeded.
pub fn setupDebugConsole() bool {
    // Prefer reusing whatever console launched us (e.g. `zig build run --
    // --debug` from an existing terminal) so output lands there instead of
    // popping a separate window. Falls back to allocating our own when
    // there's no parent console to attach to (double-click/shortcut launch).
    const attached = kernel32.AttachConsole(console.ATTACH_PARENT_PROCESS) != 0;
    if (!attached and kernel32.AllocConsole() == 0) return false;

    const conout = openHandle(L("CONOUT$"), true, false);
    if (invalidHandle(conout)) return false;

    setStdHandles(conout);
    return true;
}

fn invalidHandle(handle: win32.foundation.HANDLE) bool {
    return handle == win32.foundation.INVALID_HANDLE_VALUE;
}

/// No console at all (or "--debug" wasn't passed, or setup above failed):
/// point std handles at the NUL device instead of leaving them unset, so
/// any live-PEB-reading print always has somewhere safe to write.
fn setStdHandlesToNul() void {
    const nul = openHandle(L("NUL"), false, false);
    if (invalidHandle(nul)) return;
    setStdHandles(nul);
}

/// Temporary smoke-test diagnostic: points STD_OUTPUT_HANDLE/STD_ERROR_HANDLE
/// at a real log file, same mechanism as console.zig's NUL redirect, so
/// std.debug.print and error return traces are visible regardless of
/// whether a console is attached.
pub fn redirectStdioToFile(alloc: zigstd.Allocator, name: []const u8) !void {
    if (sysinfo.target.os.tag != .windows) return;
    const name_w = try bytes.utf8ToWinStrAllocZ(alloc, name);
    defer alloc.free(name_w);
    const handle = openHandle(name_w, false, true);
    if (invalidHandle(handle)) return error.InvalidHandle;
    setStdHandles(handle);
}

/// `SetStdHandle` updates the PEB fields Zig prints through; `AllocConsole` and
/// `AttachConsole` alone leave a GUI-subsystem process's handles invalid.
fn setStdHandles(handle: win32.foundation.HANDLE) void {
    _ = kernel32.SetStdHandle(console.STD_OUTPUT_HANDLE, handle);
    _ = kernel32.SetStdHandle(console.STD_ERROR_HANDLE, handle);
}

fn openHandle(name: [:0]const u16, read_too: bool, create: bool) win32.foundation.HANDLE {
    // GENERIC_WRITE (plus GENERIC_READ for the console case), spelled out
    // field-by-field since zigwin32 flag structs are packed structs, not
    // plain ints to `|`.
    var access = fs.FILE_ACCESS_FLAGS{
        .FILE_WRITE_DATA = 1,
        .FILE_APPEND_DATA = 1,
        .FILE_WRITE_EA = 1,
        .FILE_WRITE_ATTRIBUTES = 1,
        .READ_CONTROL = 1,
        .SYNCHRONIZE = 1,
    };
    if (read_too) {
        access.FILE_READ_DATA = 1;
        access.FILE_READ_EA = 1;
        access.FILE_READ_ATTRIBUTES = 1;
    }
    const share_all = fs.FILE_SHARE_MODE{ .READ = 1, .WRITE = 1 };

    return kernel32.CreateFileW(name, access, share_all, null, if (create) fs.CREATE_ALWAYS else fs.OPEN_EXISTING, .{}, null);
}
