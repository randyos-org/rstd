const std = @import("std");
const uefi = std.os.uefi;
const Io = std.Io;
const log = std.log.scoped(.bootfs);

const state = @import("state.zig");
const fdtable = @import("fdtable.zig");

/// Locate Simple File System protocol, open root volume into `t.root_dir`,
/// return handle-less `Io.Dir` token (see state.zig). Idempotent.
pub fn openRootDir() Io.Dir {
    if (state.root_dir == null) {
        const bs = state.bootServices() orelse {
            log.err("boot services unavailable", .{});
            return .{ .handle = {} };
        };

        log.debug("locating simple file system protocol", .{});
        const res = bs.locateProtocol(uefi.protocol.SimpleFileSystem, null) catch |err| {
            log.err("locating simple file system protocol failed: {s}", .{@errorName(err)});
            return .{ .handle = {} };
        };
        const file_system = res orelse {
            log.err("simple file system protocol not found!", .{});
            return .{ .handle = {} };
        };

        log.debug("opening root volume", .{});
        state.root_dir = file_system.openVolume() catch |err| {
            log.err("opening root volume failed: {s}", .{@errorName(err)});
            return .{ .handle = {} };
        };
    }
    // return .{ .handle = state.root_dir };
    return .{ .handle = {} };
}

// ------------------------------
// io vtable callbacks below here

pub const dirCreateDir = Io.failingDirCreateDir;
pub const dirCreateDirPath = Io.failingDirCreateDirPath;
pub const dirCreateDirPathOpen = Io.failingDirCreateDirPathOpen;
pub const dirOpenDir = Io.failingDirOpenDir;
pub const dirStat = Io.failingDirStat;
pub const dirStatFile = Io.failingDirStatFile;
pub const dirAccess = Io.failingDirAccess;
pub const dirCreateFile = Io.failingDirCreateFile;
pub const dirCreateFileAtomic = Io.failingDirCreateFileAtomic;

pub fn dirOpenFile(userdata: ?*anyopaque, _: Io.Dir, sub_path: []const u8, options: Io.Dir.OpenFileOptions) Io.File.OpenError!Io.File {
    _ = userdata;
    // only dir that exists (see state.zig); read-only
    if (options.mode != .read_only) return error.ReadOnlyFileSystem;
    const root = state.root_dir orelse return error.Unexpected; // openRootDir was never called

    // UEFI wants UCS-2 with backslash separators.
    var path_utf16: [128:0]u16 = undefined;
    if (sub_path.len >= path_utf16.len) return error.NameTooLong;
    const len = std.unicode.utf8ToUtf16Le(&path_utf16, sub_path) catch return error.BadPathName;
    for (path_utf16[0..len]) |*unit| {
        if (unit.* == '/') unit.* = '\\';
    }
    path_utf16[len] = 0;

    const opened = root.open(path_utf16[0..len :0], .read, .{}) catch |err| return switch (err) {
        error.NotFound => error.FileNotFound,
        error.AccessDenied => error.AccessDenied,
        error.OutOfResources => error.SystemResources,
        error.WriteProtected => error.ReadOnlyFileSystem,
        error.VolumeFull => error.NoSpaceLeft,
        error.InvalidParameter => error.BadPathName,
        error.NoMedia, error.MediaChanged, error.DeviceError, error.VolumeCorrupted => error.NoDevice,
        error.Unexpected => error.Unexpected,
    };
    const fd = fdtable.alloc(opened) orelse {
        opened.close() catch {};
        return error.SystemFdQuotaExceeded;
    };
    return .{ .handle = fd, .flags = .{ .nonblocking = false } };
}

pub fn dirClose(userdata: ?*anyopaque, _: []const Io.Dir) void {
    _ = userdata;
    // every Io.Dir is boot volume root (see state.zig)
    if (state.root_dir) |root| {
        root.close() catch {};
        state.root_dir = null;
    }
}

pub const dirRead = Io.noDirRead;
pub const dirRealPath = Io.failingDirRealPath;
pub const dirRealPathFile = Io.failingDirRealPathFile;
pub const dirDeleteFile = Io.failingDirDeleteFile;
pub const dirDeleteDir = Io.failingDirDeleteDir;
pub const dirRename = Io.failingDirRename;
pub const dirRenamePreserve = Io.failingDirRenamePreserve;
pub const dirSymLink = Io.failingDirSymLink;
pub const dirReadLink = Io.failingDirReadLink;
pub const dirSetOwner = Io.failingDirSetOwner;
pub const dirSetFileOwner = Io.failingDirSetFileOwner;
pub const dirSetPermissions = Io.failingDirSetPermissions;
pub const dirSetFilePermissions = Io.failingDirSetFilePermissions;
pub const dirSetTimestamps = Io.noDirSetTimestamps;
pub const dirHardLink = Io.failingDirHardLink;
