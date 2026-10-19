//! PID 1 of the Linux qemu rigs (`build/qemutest.zig`): mounts what tests
//! read (`/proc`, `/sys`, `/dev`), runs `/test`, reports how it ended, and
//! restarts the machine, which `-no-reboot` turns into qemu exiting.

const std = @import("std");
const linux = std.os.linux;

/// What the rig's run step looks for; anything else fails it.
pub const passed = "linuxtest: /test exited 0";

pub fn main() void {
    mount("proc", "/proc", "proc");
    mount("sysfs", "/sys", "sysfs");
    mount("devtmpfs", "/dev", "devtmpfs");
    mount("tmpfs", "/tmp", "tmpfs");

    // Through a pipe, not the console: on a terminal the test runner redraws
    // one progress line, which the serial log barely shows; on a pipe it
    // prints a line per test.
    var fds: [2]i32 = undefined;
    if (linux.errno(linux.pipe2(&fds, .{})) != .SUCCESS) say("linuxtest: no pipe; output stays on the console\n");
    const pid = linux.fork();
    if (pid == 0) {
        _ = linux.dup3(fds[1], 1, 0);
        _ = linux.dup3(fds[1], 2, 0);
        _ = linux.close(fds[0]);
        _ = linux.close(fds[1]);
        const argv = [_:null]?[*:0]const u8{"/test"};
        // RSTD_TEST_PROGRESS: tests that run every core name each one (`blas1/check.zig`)
        const envp = [_:null]?[*:0]const u8{ "PATH=/", "HOME=/", "TMPDIR=/tmp", "TERM=dumb", "RSTD_TEST_PROGRESS=1" };
        _ = linux.execve("/test", &argv, &envp);
        say("linuxtest: cannot run /test\n");
        linux.exit(127);
    }
    _ = linux.close(fds[1]);
    var buf: [4096]u8 = undefined;
    while (true) {
        const n = linux.read(fds[0], &buf, buf.len);
        if (linux.errno(n) == .INTR) continue;
        if (linux.errno(n) != .SUCCESS or n == 0) break;
        say(buf[0..n]);
    }
    var status: i32 = 0;
    while (true) {
        const r = linux.wait4(@intCast(pid), &status, 0, null);
        if (linux.errno(r) != .INTR) break;
    }
    const s: u32 = @bitCast(status);
    var line: [64]u8 = undefined;
    say(if (linux.W.IFEXITED(s))
        std.fmt.bufPrint(&line, "linuxtest: /test exited {d}\n", .{linux.W.EXITSTATUS(s)}) catch unreachable
    else
        std.fmt.bufPrint(&line, "linuxtest: /test killed by signal {d}\n", .{s & 0x7f}) catch unreachable);
    linux.sync();
    _ = linux.reboot(.MAGIC1, .MAGIC2, .RESTART, null);
    while (true) {}
}

fn mount(special: [*:0]const u8, dir: [*:0]const u8, fstype: [*:0]const u8) void {
    if (linux.errno(linux.mount(special, dir, fstype, 0, 0)) != .SUCCESS) {
        say("linuxtest: mount failed: ");
        say(std.mem.span(dir));
        say("\n");
    }
}

fn say(s: []const u8) void {
    _ = linux.write(1, s.ptr, s.len);
}
