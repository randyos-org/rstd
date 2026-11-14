//! The host tools the build runs, each built once per build, and the Mac they reach by default.

const std = @import("std");
const rstdbuild = @import("../src/buildutils/__root__.zig");
const Build = rstdbuild.Build;
const Compile = rstdbuild.CompileStep;

/// The SSH alias of the Mac the rigs and SSH steps use by default (`-Dtest-host`, `-Dssh-mac`).
pub const mac_alias = "macmini-m6";
/// The SSH alias of the Intel Mac the x86_64 SSH step uses by default (`-Dssh-mac-intel`).
pub const intel_mac_alias = "macmini-intel";

pub const Tool = enum { ssh_test, test_host, qemu_watch, load_window, initramfs, macos_vm };

var built: std.EnumArray(Tool, ?*Compile) = .initFill(null);

/// Where `tool`'s source is: the rigs' tools in `tools/rig/`, the macOS VM's
/// in `tools/macos/`.
fn dir(tool: Tool) []const u8 {
    return switch (tool) {
        .ssh_test => "tools",
        .test_host, .qemu_watch, .load_window, .initramfs => "tools/rig",
        .macos_vm => "tools/macos",
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

/// `tools/nix/blas.sh <system>` on `host` (`ssh_test`), at the nixpkgs of
/// `tools/nix/flake.lock`: nixpkgs' OpenBLAS and BLIS for `system`, and with
/// `reference` the bench's reference binaries (`tools/blasbench/nix`).
pub fn nixBlas(b: *Build, host: []const u8, system: []const u8, reference: bool) !*Build.Step.Run {
    const lock = try std.Io.Dir.cwd().readFileAlloc(b.graph.io, try b.root.joinString(b.allocator, "tools/nix/flake.lock"), b.allocator, .limited(1 << 16));
    const Lock = struct { nodes: struct { nixpkgs: struct { locked: struct { rev: []const u8 } } } };
    const parsed = try std.json.parseFromSliceLeaky(Lock, b.allocator, lock, .{ .ignore_unknown_fields = true });
    const run = b.addRunArtifact(hostTool(b, .ssh_test));
    run.addArg(host);
    run.addFileArg(b.path("tools/nix/blas.sh"));
    run.addArgs(&.{ system, parsed.nodes.nixpkgs.locked.rev });
    if (reference) {
        run.addPrefixedDirectoryArg("local:", b.path("tools/blasbench/nix"));
        // a directory argument's contents are not hashed
        for ([_][]const u8{ "eigen.nix", "eigen_shim.cpp" }) |f| run.addFileInput(b.path(b.fmt("tools/blasbench/nix/{s}", .{f})));
    }
    return run;
}
