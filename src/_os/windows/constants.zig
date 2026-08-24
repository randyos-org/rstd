const std = @import("std");
const win32 = @import("win32").everything;

// zigwin32 skips generating this constant since its value is NULL
pub const HWND_TOP = null;

pub const skip_msg_debug_output = .{
    win32.WM_MOUSEMOVE,
    win32.WM_NCHITTEST,
    win32.WM_GETMINMAXINFO,
    win32.WM_WINDOWPOSCHANGED,
    win32.WM_WINDOWPOSCHANGING,
    win32.WM_SETCURSOR,
    win32.WM_NCMOUSEMOVE,
};

fn errorCheck(return_value: anytype, _error: anyerror) !void {
    const failed = switch (@TypeOf(return_value)) {
        u16, win32.LRESULT => return_value == 0,
        win32.BOOL => return_value == 0,
        else => @compileError("Cannot check this type"),
    };
    if (failed) {
        std.debug.print("win32 error: {d}\n", .{win32.GetLastError()});
        return _error;
    }
}
