//! `hint-refs`: FreeType's outlines of each test font, from nixpkgs' FreeType on a Mac (`tools/nix/hintrefs.sh`), written to
//! `src/rgl/font/hint/_test/refs/` for the hinting tests (issue #805).

const rstdbuild = @import("../src/buildutils/__root__.zig");
const Build = rstdbuild.Build;
const testfonts = @import("testfonts.zig");
const tools = @import("tools.zig");

pub fn add(b: *Build, host: []const u8) !void {
    const step = b.step("hint-refs", "Regenerate the hinting tests' FreeType outlines (src/rgl/font/hint/_test/refs/) with Nix on -Dssh-mac-intel; needs -Dtest-fonts");
    if (!testfonts.enabled) {
        step.dependOn(&b.addFail("hint-refs needs -Dtest-fonts").step);
        return;
    }
    const update = b.addUpdateSourceFiles();
    step.dependOn(&update.step);
    // 26.05-darwin's nixpkgs on either Mac, so either writes the same files
    const rev = (try tools.lockRevs(b)).darwin;
    for (testfonts.fonts) |f| {
        const dep = b.lazyDependency(f.dep, .{}) orelse continue;
        const run = b.addRunArtifact(tools.hostTool(b, .ssh_test));
        run.addArg(host);
        run.addFileArg(b.path("tools/nix/hintrefs.sh"));
        run.addArg(rev);
        // `ssh_test` copies a `local:` path over; a directory's contents are not hashed
        run.addDirectoryArg2(b.path("tools/hintrefs"), .{ .prefix = "local:" });
        for ([_][]const u8{ "hintrefs.c", "hintrefs.nix" }) |file| run.addFileInput(b.path(b.fmt("tools/hintrefs/{s}", .{file})));
        run.addFileArg2(dep.path(f.path), .{ .prefix = "local:" });
        update.addCopyFileToSource(run.captureStdOut(.{}), b.fmt("src/rgl/font/hint/_test/refs/{s}.bin.gz", .{f.name}));
    }
}
