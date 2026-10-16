//! Test runner for `rstd-check-uefi`: links the OS-free tests (the qemu
//! rigs' root, `_test_qemu.zig`: simd and BLAS) for UEFI without running
//! them, so the code only UEFI builds compile (`types.exprecisionOn`, the
//! arm64 ID-register probe) is checked. No machine runs them yet.

const std = @import("std");
const builtin = @import("builtin");

pub fn main() void {
    for (builtin.test_functions) |t| std.mem.doNotOptimizeAway(t.func);
}

/// Nothing here runs, so failing stand-ins do for what std would otherwise
/// take from an OS: the testing allocator and Io, pages, and (as in
/// `_test_qemu_runner.zig`) one CPU and no threads.
pub const std_os_options: std.Options.OperatingSystem = .{
    .testing = struct {
        pub const allocator = std.testing.failing_allocator;
        pub const io = std.Io.failing;
    },
    .heap = struct {
        pub const page_size_min = 4096;
        pub const page_size_max = 4096;
        pub const page_allocator = std.testing.failing_allocator;
        pub fn queryPageSize() usize {
            return 4096;
        }
    },
    .Thread = struct {
        pub const Id = usize;
        pub const ThreadHandle = void;
        pub fn getCurrentId() Id {
            return 0;
        }
        pub fn getCpuCount() std.Thread.CpuCountError!usize {
            return 1;
        }
        pub fn yield() std.Thread.YieldError!void {}
        pub fn spawn(_: std.Thread.SpawnConfig, comptime _: anytype, _: anytype) std.Thread.SpawnError!@This() {
            return error.SystemResources;
        }
        pub fn getHandle(_: @This()) ThreadHandle {}
        pub fn detach(_: @This()) void {}
        pub fn join(_: @This()) void {}
    },
};

/// Nothing is printed: the default Io would bring in the threaded one.
pub const std_options_debug_io: std.Io = std.Io.failing;
