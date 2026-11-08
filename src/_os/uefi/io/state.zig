const std = @import("std");
const uefi = std.os.uefi;
const Io = std.Io;

const logging = @import("../../../log.zig");
const log = logging.Logger(@This());

const time = @import("time.zig");
const stderr_mod = @import("stderr.zig");
const random = @import("random.zig");
const dir = @import("dir.zig");
const file = @import("file.zig");
const async = @import("async.zig");
const stubs = @import("stubs.zig");
const operate = @import("operate.zig");
const fdtable = @import("fdtable.zig");

/// set by init(), cleared by stop(); null = console output dropped
pub var con_out: ?*uefi.protocol.SimpleTextOutput = null;

/// set by stop(); guards all Boot Services use so late print/panic after
/// exitBootServices() drops output instead of touching reclaimed memory.
/// Not derived from `boot_services == null`; don't trust firmware for that.
pub var stopped: bool = false;

/// opened boot-volume root dir (see dir.openRootDir). single slot, not a
/// table: there's only ever one `Io.Dir` here, so its handle value is just a
/// placeholder (see dir.zig's `no_handle`), not a lookup key.
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
        log.err(@src(), "starting bootloader timer failed: {s}", .{@errorName(err)});
    };
}

/// sever from Boot Services. call right after exitBootServices() succeeds;
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

/// Comptime-known `Io`: a null `userdata` and a fixed vtable, so it can be
/// `std.Options.debug_io` (RandyOS's `src/boot/uefi/__root__.zig`). `init`
/// runs before any real output.
pub const io_inst: Io = .{
    .userdata = null, //Io.Threaded.global_single_threaded,
    .vtable = &vtable,
};

/// Each `Io.VTable` entry, from the one part that declares its name.
const vtable_parts = .{ async, operate, dir, file, stubs, stderr_mod, time, random };

const vtable: Io.VTable = blk: {
    var vt: Io.VTable = undefined;
    for (@typeInfo(Io.VTable).@"struct".field_names) |name| {
        @field(vt, name) = for (vtable_parts) |part| {
            if (@hasDecl(part, name)) break @field(part, name);
        } else @compileError("no UEFI Io entry for " ++ name);
    }
    break :blk vt;
};

pub var stderr_writer: Io.File.Writer = .{
    .io = io_inst,
    .file = fdtable.stderr(),
    .mode = .streaming_simple,
    .interface = .{
        .vtable = &stderr_mod.stderr_writer_vtable,
        .buffer = &.{},
    },
};
