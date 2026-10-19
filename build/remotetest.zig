//! `remote-test-<target>` steps: rstd's full test suite, cross-compiled here
//! and run on a real machine over SSH (`tools/remote_test.zig`), for targets
//! no emulator covers. The machine is an SSH alias, so `ssh <alias>` must
//! work without a password from the shell the build runs in; its address
//! and login live in that SSH config, not here.

const std = @import("std");
const rstdbuild = @import("../src/buildutils/__root__.zig");
const Build = rstdbuild.Build;
const Optimize = rstdbuild.Optimize;
const buildroot = @import("__root__.zig");
const qemutest = @import("qemutest.zig");

/// Adds the remote steps; `filters` is `-Drstd-test-filter`.
pub fn addRemoteTests(b: *Build, optimize: Optimize, filters: []const []const u8) !void {
    // Apple silicon: the cores only real hardware runs (vortex, vortexm4, SME)
    const mac = b.option([]const u8, "remote-mac", "SSH alias of the Apple silicon Mac for remote-test-aarch64-macos (default macmini-m6)") orelse "macmini-m6";
    const target = b.resolveTargetQuery(.{ .cpu_arch = .aarch64, .os_tag = .macos });
    const mod = try buildroot.addRstdModule(b, null, target, optimize, .all, false, false);
    const test_exe = try buildroot.addRstdTest(b, mod, target, filters);

    const tool = b.addExecutable(.{ .name = "remote_test", .root_module = b.createModule(.{
        .root_source_file = b.path("tools/remote_test.zig"),
        .target = b.graph.host,
    }) });
    // under the qemu rigs' watchdog: an hour with no output ends it
    const run = b.addRunArtifact(qemutest.watchTool(b));
    run.addArgs(&.{ "60", "--" });
    run.addArtifactArg(tool);
    run.addArg(mac);
    run.addArtifactArg(test_exe);
    run.stdio = .inherit;
    run.has_side_effects = true;
    b.step("remote-test-aarch64-macos", b.fmt("Run rstd's tests on the Mac at SSH alias {s} (-Dremote-mac)", .{mac})).dependOn(&run.step);
}
