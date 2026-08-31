const std = @import("std");
const sysinfo = @import("builtin");
const os_tag = sysinfo.os.tag;
const is_windows = os_tag == .windows;
const DynLib = std.DynLib;
const win32 = @import("win32");
const Allocator = @import("zigstd.zig").Allocator;

pub const WinDll = struct {
    const Self = @This();
    handle: win32.foundation.HINSTANCE,

    pub const Error = error{FileNotFound};

    /// Trusts the file. Malicious file will be able to execute arbitrary code.
    pub fn open(path: []const u8) Error!Self {
        var buf: [std.fs.max_path_bytes:0]u8 = undefined;
        const path_z = std.mem.printSentinel(&buf, "{s}", .{path}, 0) catch return error.FileNotFound;
        return openZ(path_z);
    }

    /// Trusts the file. Malicious file will be able to execute arbitrary code.
    pub fn openZ(path_c: [*:0]const u8) Error!Self {
        const handle = win32.system.library_loader.LoadLibraryA(path_c) orelse
            return error.FileNotFound;
        return .{ .handle = handle };
    }

    /// Trusts the file.
    pub fn close(self: *Self) void {
        _ = win32.system.library_loader.FreeLibrary(self.handle);
        self.* = undefined;
    }

    pub fn lookup(self: *Self, comptime T: type, name: [:0]const u8) ?T {
        const proc = win32.system.library_loader.GetProcAddress(self.handle, name.ptr) orelse
            return null;
        return @as(T, @ptrCast(proc));
    }
};

pub const Dll = switch (os_tag) {
    .windows => WinDll,
    else => DynLib,
};

/// Returns a handle to `name` only if it is already loaded into this
/// process.
///
/// Do not call `.close()` on the result: freeing it would drop a
/// refcount this code never held, risking an early unload other code is still
/// depending on. Treat the returned handle as a permanent, borrowed reference
/// for the rest of the process's lifetime.
pub fn findLoaded(name: [:0]const u8) ?Dll {
    if (is_windows) {
        const handle = win32.system.library_loader.GetModuleHandleA(name.ptr) orelse
            return null;
        return .{ .handle = handle };
    }
    const handle = std.c.dlopen(name.ptr, .{ .LAZY = true, .NOLOAD = true }) orelse
        return null;
    return .{ .inner = .{ .handle = handle } };
}

pub fn loadDll(alloc: Allocator, lib_name: []const u8) !Dll {
    const lib_name_fmt = if (is_windows) "{s}.dll" else "lib{s}.so";
    const lib = try alloc.print(lib_name_fmt, .{lib_name});
    defer alloc.free(lib);
    return try Dll.open(lib);
}
