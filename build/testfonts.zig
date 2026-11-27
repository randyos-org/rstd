//! The fonts of rstd's font-backed tests and `hint-refs` (issue #805):
//! rgui's, `build.zig.zon`'s lazy `font_*` dependencies, fetched only with
//! `-Dtest-fonts`, as a lazy dependency is fetched whenever it is asked for.

const std = @import("std");
const rstdbuild = @import("../src/buildutils/__root__.zig");
const Build = rstdbuild.Build;
const Module = Build.Module;

pub const Font = struct {
    /// Its field in `rstd_test_fonts`, and its `hint-refs` file's name.
    name: []const u8,
    /// Its `build.zig.zon` dependency.
    dep: []const u8,
    /// The file in that dependency.
    path: []const u8,
};

pub const fonts = [_]Font{
    .{ .name = "arimo", .dep = "font_arimo", .path = "ArimoNerdFont-Regular.ttf" },
    .{ .name = "jetbrains_mono", .dep = "font_jetbrains_mono", .path = "JetBrainsMonoNerdFontMono-Regular.ttf" },
    .{ .name = "ubuntu", .dep = "font_ubuntu", .path = "UbuntuNerdFont-Regular.ttf" },
    .{ .name = "ibm_plex_serif", .dep = "font_ibm_plex_serif", .path = "fonts/complete/ttf/IBMPlexSerif-Regular.ttf" },
    .{ .name = "commit_mono", .dep = "font_commit_mono", .path = "CommitMonoNerdFontMono-Regular.otf" },
    // the release's own name for Atkinson Hyperlegible Mono
    .{ .name = "atkynson_mono", .dep = "font_atkinson_hyperlegible_mono", .path = "AtkynsonMonoNerdFontMono-Regular.otf" },
    // CJK, for the autohinter's CJK writing system
    .{ .name = "cica", .dep = "font_cica", .path = "Cica-Regular.ttf" },
};

/// `-Dtest-fonts`, which `build` reads before any test is added.
pub var enabled = false;

var cached: ?*Module = null;

/// `rstd_test_fonts`: each font's bytes as a `?[]const u8` named for it,
/// null without `-Dtest-fonts` (or until the fonts are fetched).
pub fn module(b: *Build) !*Module {
    if (cached) |m| if (m.owner == b) return m;
    const wf = b.addWriteFiles();
    var src: std.ArrayList(u8) = .empty;
    try src.appendSlice(b.allocator, "//! The test fonts (build/testfonts.zig).\n\n");
    for (fonts) |f| {
        const dep = if (enabled) b.lazyDependency(f.dep, .{}) else null;
        if (dep) |d| {
            const file = std.Io.Dir.path.basename(f.path);
            _ = wf.addCopyFile(d.path(f.path), file);
            try src.print(b.allocator, "pub const {s}: ?[]const u8 = @embedFile(\"{s}\");\n", .{ f.name, file });
        } else try src.print(b.allocator, "pub const {s}: ?[]const u8 = null;\n", .{f.name});
    }
    const m = b.createModule(.{ .root_source_file = wf.add("test_fonts.zig", src.items) });
    cached = m;
    return m;
}
