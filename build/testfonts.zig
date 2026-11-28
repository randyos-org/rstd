//! The fonts of rstd's font-backed tests and `hint-refs` (issue #805):
//! rgui's, `build.zig.zon`'s lazy `font_*` dependencies, fetched only with
//! `-Dtest-fonts`, as a lazy dependency is fetched whenever it is asked for;
//! and synthetic CFF fonts `tools/cfffonts.zig` writes (issue #832).

const std = @import("std");
const rstdbuild = @import("../src/buildutils/__root__.zig");
const Build = rstdbuild.Build;
const Module = Build.Module;
const LazyPath = Build.LazyPath;
const tools = @import("tools.zig");

pub const Font = struct {
    /// Its field in `rstd_test_fonts`, and its `hint-refs` file's name.
    name: []const u8,
    /// Its `build.zig.zon` dependency.
    dep: []const u8 = "",
    /// The file in that dependency, or the name of the one generated.
    path: []const u8,
    /// `tools/cfffonts.zig`'s font of this name instead of a dependency.
    generated: ?[]const u8 = null,
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
    // variable CFF2 (wght 200 to 900, avar, HVAR), for blended charstrings and Private DICTs
    .{ .name = "source_sans_3", .dep = "font_source_sans_3_vf", .path = "VF/SourceSans3VF-Upright.otf" },
    // CID-keyed with FontMatrix shears and units per em, and SEAC accents
    .{ .name = "cff_cid", .path = "rstd-cff-cid.otf", .generated = "cid" },
    .{ .name = "cff_seac", .path = "rstd-cff-seac.otf", .generated = "seac" },
};

/// A `hint-refs` file: a font at an instance, its axis values in 16.16
/// (`FT_Set_Var_Design_Coordinates`), none for the default.
pub const Ref = struct {
    /// The file's name.
    name: []const u8,
    font: []const u8,
    design: []const i32 = &.{},
};

/// Each font at its default, and the variable ones at named and unnamed weights.
pub const refs = blk: {
    var r: [fonts.len + 6]Ref = undefined;
    for (fonts, 0..) |f, i| r[i] = .{ .name = f.name, .font = f.name };
    r[fonts.len] = .{ .name = "cascadia_code_wght300", .font = "cascadia_code", .design = &.{300 << 16} };
    r[fonts.len + 1] = .{ .name = "cascadia_code_wght550", .font = "cascadia_code", .design = &.{550 << 16} };
    r[fonts.len + 2] = .{ .name = "cascadia_code_wght700", .font = "cascadia_code", .design = &.{700 << 16} };
    r[fonts.len + 3] = .{ .name = "source_sans_3_wght400", .font = "source_sans_3", .design = &.{400 << 16} };
    r[fonts.len + 4] = .{ .name = "source_sans_3_wght650", .font = "source_sans_3", .design = &.{650 << 16} };
    r[fonts.len + 5] = .{ .name = "source_sans_3_wght900", .font = "source_sans_3", .design = &.{900 << 16} };
    break :blk r;
};

/// `f`'s file, null until fetched.
pub fn file(b: *Build, f: Font) ?LazyPath {
    if (f.generated) |which| {
        const run = b.addRunArtifact(tools.hostTool(b, .cfffonts));
        run.addArg(which);
        return run.addOutputFileArg(f.path);
    }
    const dep = b.lazyDependency(f.dep, .{}) orelse return null;
    return dep.path(f.path);
}

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
        const path = if (enabled) file(b, f) else null;
        if (path) |p| {
            const base = std.Io.Dir.path.basename(f.path);
            _ = wf.addCopyFile(p, base);
            try src.print(b.allocator, "pub const {s}: ?[]const u8 = @embedFile(\"{s}\");\n", .{ f.name, base });
        } else try src.print(b.allocator, "pub const {s}: ?[]const u8 = null;\n", .{f.name});
    }
    const m = b.createModule(.{ .root_source_file = wf.add("test_fonts.zig", src.items) });
    cached = m;
    return m;
}
