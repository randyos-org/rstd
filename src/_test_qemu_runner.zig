//! Freestanding test runner for `qemu-system-*`. Prints results through the
//! architecture's semihosting-style console and exits with the failure count
//! (0 = all passed), which the run step maps to the process exit code.
//!
//! Use with `zig test -target <arch>-freestanding --test-runner runner.zig`
//! and the matching `_qemutest/<arch>.ld`; `<arch>.zig` supplies the boot code.

const std = @import("std");
const sysinfo = @import("builtin");
const mem = @import("mem.zig");

const arch = switch (sysinfo.cpu.arch) {
    .aarch64 => @import("_qemutest/aarch64.zig"),
    .riscv64 => @import("_qemutest/riscv64.zig"),
    .arm => @import("_qemutest/arm.zig"),
    .loongarch64 => @import("_qemutest/loongarch64.zig"),
    .mips64el => @import("_qemutest/mips64.zig"),
    .s390x => @import("_qemutest/s390x.zig"),
    .powerpc64le => @import("_qemutest/ppc64le.zig"),
    else => @compileError("qemutest: no boot code for " ++ @tagName(sysinfo.cpu.arch)),
};

pub const panic = std.debug.no_panic;

/// Freestanding targets have no page size to infer.
pub const std_os_options: std.Options.OperatingSystem = .{
    .heap = struct {
        pub const page_size_min = mem.page_size;
        pub const page_size_max = mem.page_size;
    },
};

/// `std.testing` failure messages go through `std.debug.print`; there is no OS
/// to print to, and the default Io would drag in the threaded implementation.
pub const std_options_debug_io: std.Io = std.Io.failing;

comptime {
    _ = arch;
}

fn print(s: []const u8) void {
    for (s) |c| arch.writeByte(c);
}

fn printInt(v: usize) void {
    var buf: [20]u8 = undefined;
    var i: usize = buf.len;
    var n = v;
    while (true) {
        i -= 1;
        buf[i] = '0' + @as(u8, @intCast(n % 10));
        n /= 10;
        if (n == 0) break;
    }
    print(buf[i..]);
}

export fn qemutest_main() callconv(.c) noreturn {
    print("qemutest ");
    print(@tagName(sysinfo.cpu.arch));
    print("\n");
    var pass: usize = 0;
    var fail: usize = 0;
    var skip: usize = 0;
    for (sysinfo.test_functions) |t| {
        print(t.name);
        print("... ");
        if (t.func()) |_| {
            pass += 1;
            print("ok\n");
        } else |err| switch (err) {
            error.SkipZigTest => {
                skip += 1;
                print("skip\n");
            },
            else => {
                fail += 1;
                print("FAIL ");
                print(@errorName(err));
                print("\n");
            },
        }
    }
    printInt(pass);
    print(" passed, ");
    printInt(skip);
    print(" skipped, ");
    printInt(fail);
    print(" failed\n");
    arch.exit(if (fail == 0) 0 else 1);
}
