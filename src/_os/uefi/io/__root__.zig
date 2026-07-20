//! `std.Io` impl backed by UEFI protocols.
//!
//! std.debug.print/panic/std.log route through `std.Options.debug_io`.
//! Default `std.Io.Threaded` needs POSIX/Windows syscalls, doesn't compile
//! for UEFI. This replaces it; exported as `std_options_debug_io`.
//!
//! `Io.VTable` has no optional fields -- everything must be populated even
//! with no UEFI equivalent (process, net...). Real impls:
//!
//! * `lockStderr`/`tryLockStderr`/`unlockStderr` -- Simple Text Output
//!   Protocol (`con_out`).
//! * `dirOpenFile`/`fileReadPositional`/`fileLength`/`fileClose`/`dirClose`
//!   -- Simple File System protocol on boot volume (see `openRootDir`). Only
//!   ONE dir + ONE file open at a time: `Io.Dir`/`Io.File` have no handle
//!   bits on UEFI (`posix.fd_t` is `void`) -- see state.zig.
//! * `now`/`clockResolution` -- `RuntimeServices.getTime()` for `.real`
//!   (whole-second RTC), tick counter (time.zig) otherwise.
//! * `sleep` -- `BootServices.stall()`.
//! * `random`/`randomSecure` -- EFI RNG protocol, timer-seeded PRNG fallback
//!   for infallible `random`.
//! * `async`/`groupAsync` -- run eagerly on calling "thread" (only option,
//!   single-threaded).
//! * Cancellation/futex -- no-ops; one task, nothing to cancel/wake.
//!
//! Everything else panics with the missing op's name -- loud failure beats
//! silent garbage.
//!
//! Boot Services lifetime: call `init()` early (wires `con_out`, tick
//! clock), `stop()` right after `exitBootServices()` -- after that, console
//! output drops and stall/RNG degrade to no-ops (backing protocols gone).

const std = @import("std");
const Io = std.Io;

pub const time = @import("time.zig");

const state = @import("state.zig");
const dir = @import("dir.zig");
const fdtable = @import("fdtable.zig");
// const stderr_mod = @import("stderr.zig");
// const random = @import("random.zig");
// const net = @import("net.zig");
// const file = @import("file.zig");
// const progress = @import("progress.zig");
// const async = @import("async.zig");
// const proc = @import("proc.zig");
// const operate = @import("operate.zig");

// test "hello" {
//     try std.debug.print("hello\n", .{});
// }

pub const stdout = fdtable.stdout;
pub const stderr = fdtable.stderr;
pub const stdin = fdtable.stdin;
pub const cwdFn: ?fn () std.Io.Dir = dir.openRootDir;
pub const init: ?fn () void = state.init;
pub const stop: ?fn () void = state.stop;
pub const io_inst: Io = state.io_inst;

pub fn ioFactory() Io {
    return io_inst;
}
