const std = @import("std");
const uefi = std.os.uefi;
const Io = std.Io;
const log = std.log.scoped(.iostate);

const time = @import("time.zig");
const stderr_mod = @import("stderr.zig");
const random = @import("random.zig");
const dir = @import("dir.zig");
const net = @import("net.zig");
const file = @import("file.zig");
const progress = @import("progress.zig");
const async = @import("async.zig");
const proc = @import("proc.zig");
const operate = @import("operate.zig");
const stream = @import("stream.zig");
const fdtable = @import("fdtable.zig");

/// set by init(), cleared by stop(); null = console output dropped
pub var con_out: ?*uefi.protocol.SimpleTextOutput = null;

/// set by stop(); guards all Boot Services use so late print/panic after
/// exitBootServices() drops output instead of touching reclaimed memory.
/// Not derived from `boot_services == null` -- don't trust firmware for that.
pub var stopped: bool = false;

/// opened boot-volume root dir (see dir.openRootDir). single slot, not a
/// table: `Io.Dir` has no handle bits on UEFI, so values can't be told apart.
pub var root_dir: ?*uefi.protocol.File = null;

/// per-task cancel-protection state; only one task here
pub var cancel_protection: Io.CancelProtection = .unblocked;

/// recursion depth not a mutex: single thread, but lock is recursive so
/// panic handler can re-enter while a log line holds it.
pub var stderr_lock_count: u32 = 0;

/// wire UEFI console + tick clock in. call once, before anything logs.
pub fn init() void {
    if (con_out == null) {
        if (uefi.system_table.con_out) |out| {
            out.reset(false) catch {};
            con_out = out;
        }
    }
    time.init() catch |err| {
        // not fatal: log lines just lack a monotonic clock; console's already up
        log.err("starting bootloader timer failed: {s}", .{@errorName(err)});
    };
}

/// sever from Boot Services. call right after exitBootServices() succeeds --
/// con_out/stall/RNG live in memory firmware may reclaim then.
pub fn stop() void {
    con_out = null;
    // no close() calls: protocols already gone, just drop dangling pointers
    root_dir = null;
    fdtable.reset();
    stopped = true;
}

pub fn bootServices() ?*uefi.tables.BootServices {
    if (stopped) return null;
    return uefi.system_table.boot_services;
}

/// Comptime-known `Io` value -- `userdata` is a fixed global pointer and
/// `vtable` a fixed function table, so this is safe to use directly as
/// `std_options_debug_io` (evaluated at comptime; see `std.Options.debug_io`
/// and `src/boot/uefi/__root__.zig`). Runtime state (the stderr writer,
/// `state.io`) is wired separately by `io()`, which must run before this is
/// used for real output.
pub const io_inst: Io = .{
    .userdata = null, //Io.Threaded.global_single_threaded,
    .vtable = &.{
        .crashHandler = async.crashHandler,
        .async = async.async,
        .concurrent = async.concurrent,
        .await = async.await,
        .cancel = async.cancel,
        .groupAsync = async.groupAsync,
        .groupConcurrent = async.groupConcurrent,
        .groupAwait = async.groupAwait,
        .groupCancel = async.groupCancel,
        .recancel = async.recancel,
        .swapCancelProtection = async.swapCancelProtection,
        .checkCancel = async.checkCancel,
        .futexWait = async.futexWait,
        .futexWaitUncancelable = async.futexWaitUncancelable,
        .futexWake = async.futexWake,

        .operate = operate.operate,

        .batchAwaitAsync = async.batchAwaitAsync,
        .batchAwaitConcurrent = async.batchAwaitConcurrent,
        .batchCancel = async.batchCancel,

        .dirCreateDir = dir.dirCreateDir,
        .dirCreateDirPath = dir.dirCreateDirPath,
        .dirCreateDirPathOpen = dir.dirCreateDirPathOpen,
        .dirOpenDir = dir.dirOpenDir,
        .dirStat = dir.dirStat,
        .dirStatFile = dir.dirStatFile,
        .dirAccess = dir.dirAccess,
        .dirCreateFile = dir.dirCreateFile,
        .dirCreateFileAtomic = dir.dirCreateFileAtomic,
        .dirOpenFile = dir.dirOpenFile,
        .dirClose = dir.dirClose,
        .dirRead = dir.dirRead,
        .dirRealPath = dir.dirRealPath,
        .dirRealPathFile = dir.dirRealPathFile,
        .dirDeleteFile = dir.dirDeleteFile,
        .dirDeleteDir = dir.dirDeleteDir,
        .dirRename = dir.dirRename,
        .dirRenamePreserve = dir.dirRenamePreserve,
        .dirSymLink = dir.dirSymLink,
        .dirReadLink = dir.dirReadLink,
        .dirSetOwner = dir.dirSetOwner,
        .dirSetFileOwner = dir.dirSetFileOwner,
        .dirSetPermissions = dir.dirSetPermissions,
        .dirSetFilePermissions = dir.dirSetFilePermissions,
        .dirSetTimestamps = dir.dirSetTimestamps,
        .dirHardLink = dir.dirHardLink,

        .fileStat = file.fileStat,
        .fileLength = file.fileLength,
        .fileClose = file.fileClose,
        .fileWritePositional = file.fileWritePositional,
        .fileWriteFileStreaming = file.fileWriteFileStreaming,
        .fileWriteFilePositional = file.fileWriteFilePositional,
        .fileReadPositional = file.fileReadPositional,
        .fileSeekBy = file.fileSeekBy,
        .fileSeekTo = file.fileSeekTo,
        .fileSync = file.fileSync,
        .fileIsTty = file.fileIsTty,
        .fileEnableAnsiEscapeCodes = file.fileEnableAnsiEscapeCodes,
        .fileSupportsAnsiEscapeCodes = file.fileSupportsAnsiEscapeCodes,
        .fileSetLength = file.fileSetLength,
        .fileSetOwner = file.fileSetOwner,
        .fileSetPermissions = file.fileSetPermissions,
        .fileSetTimestamps = file.fileSetTimestamps,
        .fileLock = file.fileLock,
        .fileTryLock = file.fileTryLock,
        .fileUnlock = file.fileUnlock,
        .fileDowngradeLock = file.fileDowngradeLock,
        .fileRealPath = file.fileRealPath,
        .fileHardLink = file.fileHardLink,

        .fileMemoryMapCreate = file.fileMemoryMapCreate,
        .fileMemoryMapDestroy = file.fileMemoryMapDestroy,
        .fileMemoryMapSetLength = file.fileMemoryMapSetLength,
        .fileMemoryMapRead = file.fileMemoryMapRead,
        .fileMemoryMapWrite = file.fileMemoryMapWrite,

        .processExecutableOpen = proc.processExecutableOpen,
        .processExecutablePath = proc.processExecutablePath,

        .lockStderr = stderr_mod.lockStderr,
        .tryLockStderr = stderr_mod.tryLockStderr,
        .unlockStderr = stderr_mod.unlockStderr,

        .processCurrentPath = proc.processCurrentPath,
        .processSetCurrentDir = proc.processSetCurrentDir,
        .processSetCurrentPath = proc.processSetCurrentPath,
        .processReplace = proc.processReplace,
        .processReplacePath = proc.processReplacePath,
        .processSpawn = proc.processSpawn,
        .processSpawnPath = proc.processSpawnPath,
        .childWait = proc.childWait,
        .childKill = proc.childKill,

        .progressParentFile = progress.progressParentFile,

        .now = time.now,
        .clockResolution = time.clockResolution,
        .sleep = time.sleep,

        .random = random.random,
        .randomSecure = random.randomSecure,

        .netListenIp = net.netListenIp,
        .netAccept = net.netAccept,
        .netBindIp = net.netBindIp,
        .netConnectIp = net.netConnectIp,
        .netListenUnix = net.netListenUnix,
        .netConnectUnix = net.netConnectUnix,
        .netSocketCreatePair = net.netSocketCreatePair,
        .netSend = net.netSend,
        .netWrite = net.netWrite,
        .netWriteFile = net.netWriteFile,
        .netClose = net.netClose,
        .netShutdown = net.netShutdown,
        .netInterfaceNameResolve = net.netInterfaceNameResolve,
        .netInterfaceName = net.netInterfaceName,
        .netLookup = net.netLookup,
    },
};

pub var stderr_writer: Io.File.Writer = .{
    .io = io_inst,
    .file = stream.stderr(),
    .mode = .streaming_simple,
    .interface = .{
        .vtable = &@import("stderr.zig").stderr_writer_vtable,
        .buffer = &.{},
    },
};
