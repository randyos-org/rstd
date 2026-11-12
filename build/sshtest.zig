//! `ssh-test-<target>` steps: rstd's full suite cross-compiled here and run over SSH (`tools/ssh_test.zig`) where no emulator covers
//! the target; the machine is an SSH alias whose address stays in the SSH config. The bench steps are `bench.zig`'s.

const std = @import("std");
const rstdbuild = @import("../src/buildutils/__root__.zig");
const Build = rstdbuild.Build;
const Optimize = rstdbuild.Optimize;
const buildroot = @import("__root__.zig");
const qemutest = @import("qemu/__root__.zig");
const bench = @import("bench.zig");
const tools = @import("tools.zig");

/// Adds the SSH steps and `bench.zig`'s; `filters` is `-Drstd-test-filter`.
pub fn addSshTests(b: *Build, optimize: Optimize, filters: []const []const u8) !void {
    const mac = b.option([]const u8, "ssh-mac", "SSH alias of the Apple silicon Mac for the ssh-test steps (default " ++ tools.mac_alias ++ ")") orelse tools.mac_alias;
    const tool = tools.hostTool(b, .ssh_test);
    // Apple silicon: the cores only real hardware runs (vortex, vortexm4, SME)
    const arm = b.resolveTargetQuery(.{ .cpu_arch = .aarch64, .os_tag = .macos });
    try addSshTest(b, optimize, filters, mac, tool, arm, &.{.{ .where = "on the Mac" }});
    // x86_64 macOS under Rosetta 2, which the Mac must have installed. Its
    // cpuid hides AVX unless asked, so rstd's AVX cores run only in `-avx2`.
    // `hw.optional.avx2_0` stays 0 there, so Accelerate likely keeps its SSE
    // path: its AVX2 path needs an Intel Mac.
    const x86 = b.resolveTargetQuery(.{ .cpu_arch = .x86_64, .os_tag = .macos });
    try addSshTest(b, optimize, filters, mac, tool, x86, &.{
        .{ .where = "on the Mac under Rosetta 2" },
        .{ .suffix = "-avx2", .env = "ROSETTA_ADVERTISE_AVX=1", .where = "on the Mac under Rosetta 2 with AVX2 and FMA in cpuid (rstd's AVX cores too)" },
    });
    addGroupRepro(b, optimize, mac, tool, arm);
    try bench.addBenches(b, mac, arm, tool);
}

/// One `ssh-test-<arch>-macos<suffix>` step: where it runs, and a variable
/// set for the run.
const SshVariant = struct { suffix: []const u8 = "", env: ?[]const u8 = null, where: []const u8 };

/// `ssh-test-<arch>-macos`: the full suite for `target`, run on `mac`, once
/// per variant.
fn addSshTest(b: *Build, optimize: Optimize, filters: []const []const u8, mac: []const u8, tool: *rstdbuild.CompileStep, target: rstdbuild.ResolvedTarget, variants: []const SshVariant) !void {
    const mod = try buildroot.addRstdModule(b, null, target, optimize, .all, .{});
    const test_exe = try buildroot.addRstdTest(b, mod, target, optimize, filters);
    for (variants) |v| {
        // under the qemu rigs' watchdog: an hour with no output ends it
        const run = qemutest.watchCommand(b);
        run.addArg("--");
        run.addArtifactArg2(tool, .{});
        run.addArg(mac);
        run.addArtifactArg2(test_exe, .{});
        if (v.env) |e| run.addArg(b.fmt("env:{s}", .{e}));
        run.stdio = .inherit;
        run.has_side_effects = true;
        const name = b.fmt("ssh-test-{t}-macos{s}", .{ target.result.cpu.arch, v.suffix });
        b.step(name, b.fmt("Run rstd's {t} tests {s} at SSH alias {s} (-Dssh-mac)", .{ target.result.cpu.arch, v.where, mac })).dependOn(&run.step);
        // what that runs, built alone, so a build queue frees its slot before
        // the Mac runs the tests
        const check = b.step(b.fmt("ssh-check-{t}-macos{s}", .{ target.result.cpu.arch, v.suffix }), b.fmt("Build what {s} runs, without running it", .{name}));
        check.dependOn(&test_exe.step);
        check.dependOn(&tool.step);
    }
}

/// `ssh-repro-group-aarch64-macos`: `tools/probe/evented_group_repro.zig` on the
/// Mac, `futures` then `group`; fails while std's `Group` still crashes on
/// `Dispatch` (issue #208). Rerun on each toolchain bump.
fn addGroupRepro(b: *Build, optimize: Optimize, mac: []const u8, tool: *rstdbuild.CompileStep, target: rstdbuild.ResolvedTarget) void {
    const exe = b.addExecutable(.{ .name = "evented_group_repro", .root_module = b.createModule(.{
        .root_source_file = b.path("tools/probe/evented_group_repro.zig"),
        .target = target,
        .optimize = optimize,
    }) });
    const step = b.step("ssh-repro-group-aarch64-macos", b.fmt("Run std's Io.Group crash repro on Io.Dispatch at SSH alias {s}; passing means rstd can drop its futures fallback", .{mac}));
    var prev: ?*Build.Step = null;
    for ([_][]const u8{ "futures", "group" }) |mode| {
        const run = qemutest.watchCommand(b);
        run.addArg("--");
        run.addArtifactArg2(tool, .{});
        run.addArg(mac);
        run.addArtifactArg2(exe, .{});
        run.addArg(mode);
        run.stdio = .inherit;
        run.has_side_effects = true;
        if (prev) |p| run.step.dependOn(p);
        prev = &run.step;
    }
    step.dependOn(prev.?);
    const check = b.step("ssh-repro-group-aarch64-macos-build", "Build what ssh-repro-group-aarch64-macos runs, without running it");
    check.dependOn(&exe.step);
    check.dependOn(&tool.step);
}
