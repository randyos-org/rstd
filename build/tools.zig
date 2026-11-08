//! The host tools the build runs, each built once per build, and the Mac they reach by default.

const std = @import("std");
const rstdbuild = @import("../src/buildutils/__root__.zig");
const Build = rstdbuild.Build;
const Compile = rstdbuild.CompileStep;

/// The SSH alias of the Mac the rigs and SSH steps use by default (`-Dtest-host`, `-Dssh-mac`).
pub const mac_alias = "macmini-m6";

pub const Tool = enum { ssh_test, test_host, qemu_watch, load_window, initramfs, macos_vm };

var built: std.EnumArray(Tool, ?*Compile) = .initFill(null);

/// `tools/<tool>.zig` for the host, cached per build.
pub fn hostTool(b: *Build, comptime tool: Tool) *Compile {
    const slot = built.getPtr(tool);
    if (slot.*) |t| if (t.step.owner == b) return t;
    slot.* = b.addExecutable(.{ .name = @tagName(tool), .root_module = b.createModule(.{
        .root_source_file = b.path("tools/" ++ @tagName(tool) ++ ".zig"),
        .target = b.graph.host,
    }) });
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
