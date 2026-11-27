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
    // Limbu, for the autohinter's Indic writing system
    .{ .name = "noto_sans_limbu", .dep = "font_noto_sans_limbu", .path = "NotoSansLimbu/hinted/ttf/NotoSansLimbu-Regular.ttf" },
    // variable (wght 200 to 700, avar), VTT bytecode with cvar, for variable hinting
    .{ .name = "cascadia_code", .dep = "font_cascadia_code", .path = "ttf/CascadiaCode.ttf" },
};

/// A `hint-refs` file: a font at an instance, its axis values in 16.16
/// (`FT_Set_Var_Design_Coordinates`), none for the default.
pub const Ref = struct {
    /// The file's name.
    name: []const u8,
    font: []const u8,
    design: []const i32 = &.{},
};

/// Each font at its default, and Cascadia Code at named and unnamed weights.
pub const refs = blk: {
    var r: [fonts.len + 3]Ref = undefined;
    for (fonts, 0..) |f, i| r[i] = .{ .name = f.name, .font = f.name };
    r[fonts.len] = .{ .name = "cascadia_code_wght300", .font = "cascadia_code", .design = &.{300 << 16} };
    r[fonts.len + 1] = .{ .name = "cascadia_code_wght550", .font = "cascadia_code", .design = &.{550 << 16} };
    r[fonts.len + 2] = .{ .name = "cascadia_code_wght700", .font = "cascadia_code", .design = &.{700 << 16} };
    break :blk r;
};

pub fn find(name: []const u8) Font {
    for (fonts) |f| if (std.mem.eql(u8, f.name, name)) return f;
    unreachable;
}

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
