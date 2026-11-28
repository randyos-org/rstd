//! The host tools the build runs, each built once per build, and the Mac they reach by default.

const std = @import("std");
const rstdbuild = @import("../src/buildutils/__root__.zig");
const Build = rstdbuild.Build;
const Compile = rstdbuild.CompileStep;

/// The SSH alias of the Mac the rigs and SSH steps use by default (`-Dtest-host`, `-Dssh-mac`).
pub const mac_alias = "macmini-m6";
/// The SSH alias of the Intel Mac the x86_64 SSH step uses by default (`-Dssh-mac-intel`).
pub const intel_mac_alias = "macmini-intel";

pub const Tool = enum { ssh_test, cfffonts, test_host, qemu_watch, load_window, initramfs };

var built: std.EnumArray(Tool, ?*Compile) = .initFill(null);

/// Where `tool`'s source is: the rigs' tools in `tools/rig/`.
fn dir(tool: Tool) []const u8 {
    return switch (tool) {
        .ssh_test, .cfffonts => "tools",
        .test_host, .qemu_watch, .load_window, .initramfs => "tools/rig",
    };
}

/// `<dir>/<tool>.zig` for the host, cached per build, with `tools/remote.zig`
/// as the `remote` module (a file in a subfolder cannot import it by path).
pub fn hostTool(b: *Build, comptime tool: Tool) *Compile {
    const slot = built.getPtr(tool);
    if (slot.*) |t| if (t.step.owner == b) return t;
    const mod = b.createModule(.{
        .root_source_file = b.path(comptime dir(tool) ++ "/" ++ @tagName(tool) ++ ".zig"),
        .target = b.graph.host,
    });
    mod.addImport("remote", b.createModule(.{ .root_source_file = b.path("tools/remote.zig"), .target = b.graph.host }));
    slot.* = b.addExecutable(.{ .name = @tagName(tool), .root_module = mod });
    return slot.*.?;
}

/// `tools/nix/blas.sh <system>` on `host` (`ssh_test`), or with none on
/// this machine's Nix (CI's runners, issue #408), at the nixpkgs of
/// `tools/nix/flake.lock`: nixpkgs' OpenBLAS and BLIS for `system`, and with
/// `reference` the bench's reference binaries (`tools/blasbench/nix`).
pub fn nixBlas(b: *Build, host: ?[]const u8, system: []const u8, reference: bool) !*Build.Step.Run {
    const run = try nixScript(b, host, "blas.sh", system);
    if (reference) {
        // `ssh_test` copies a `local:` directory over
        run.addDirectoryArg2(b.path("tools/blasbench/nix"), .{ .prefix = if (host != null) "local:" else "" });
        // a directory argument's contents are not hashed
        for ([_][]const u8{ "eigen.nix", "eigen_shim.cpp" }) |f| run.addFileInput(b.path(b.fmt("tools/blasbench/nix/{s}", .{f})));
    }
    return run;
}

/// nixpkgs' prebuilt Linux kernel for `system` (`tools/nix/kernel.sh`), from
/// `host`'s Nix or this machine's, at `tools/nix/flake.lock`'s nixpkgs: the
/// Linux rigs' and `qemu-bench-linux` steps' kernel (issue #656).
pub fn nixKernel(b: *Build, host: ?[]const u8, system: []const u8) !Build.LazyPath {
    const run = try nixScript(b, host, "kernel.sh", system);
    // the nixpkgs an x86_64-darwin Nix machine evaluates instead
    run.addArg((try lockRevs(b)).darwin);
    return run.captureStdOut(.{});
}

/// `tools/nix/<script> <system> <nixpkgs rev>` on `host` (`ssh_test`), or
/// with none on this machine's Nix (CI's runners, issue #408).
fn nixScript(b: *Build, host: ?[]const u8, comptime script: []const u8, system: []const u8) !*Build.Step.Run {
    const run = if (host) |h| blk: {
        const r = b.addRunArtifact(hostTool(b, .ssh_test));
        r.addArg(h);
        break :blk r;
    } else b.addSystemCommand(&.{"bash"});
    run.addFileArg(b.path("tools/nix/" ++ script));
    run.addArgs(&.{ system, (try lockRevs(b)).nixpkgs });
    return run;
}

/// `tools/nix/flake.lock`'s nixpkgs revisions: unstable, and the release an
/// x86_64-darwin host still evaluates (`flake.nix`).
pub fn lockRevs(b: *Build) !struct { nixpkgs: []const u8, darwin: []const u8 } {
    const lock = try std.Io.Dir.cwd().readFileAlloc(b.graph.io, try b.root.joinString(b.allocator, "tools/nix/flake.lock"), b.allocator, .limited(1 << 16));
    const Rev = struct { locked: struct { rev: []const u8 } };
    const Lock = struct { nodes: struct { nixpkgs: Rev, @"nixpkgs-x86_64-darwin": Rev } };
    const parsed = try std.json.parseFromSliceLeaky(Lock, b.allocator, lock, .{ .ignore_unknown_fields = true });
    return .{ .nixpkgs = parsed.nodes.nixpkgs.locked.rev, .darwin = parsed.nodes.@"nixpkgs-x86_64-darwin".locked.rev };
}
