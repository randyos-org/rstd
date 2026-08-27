const rstdbuild = @import("__root__.zig");
const files = @import("files.zig");

const Build = rstdbuild.Build;
const LazyPath = rstdbuild.LazyPath;
const ResolvedTarget = rstdbuild.ResolvedTarget;
const allocPrint = @import("std").fmt.allocPrint;

/// add a FreeDesktop `.desktop` file and icon to the install step for Linux targets
pub fn addLinuxDesktopFile(
    b: *Build,
    target: ResolvedTarget,
    icon_png: LazyPath,
    name: []const u8,
    display_name: []const u8,
    version: []const u8,
) !void {
    if (target.result.os.tag != .linux) return;

    const desktop_name = try allocPrint(b.allocator, "{s}.desktop", .{name});
    const desktop_contents = try allocPrint(b.allocator,
        \\[Desktop Entry]
        \\Type=Application
        \\Version=1.0
        \\Name={s}
        \\Comment={s} v{s}
        \\Exec={s}
        \\Icon={s}
        \\Terminal=false
        \\Categories=Utility;
        \\
    , .{ display_name, display_name, version, name, name });
    const build_dir = files.addBuildDir(b);
    const desktop_path = build_dir.add(desktop_name, desktop_contents);

    const icon_name = try allocPrint(b.allocator, "{s}.png", .{name});
    const install_step = b.getInstallStep();
    install_step.dependOn(&b.addInstallFileWithDir(
        desktop_path,
        .prefix,
        b.pathJoin(&.{ "share/applications", desktop_name }),
    ).step);
    install_step.dependOn(&b.addInstallFileWithDir(
        icon_png,
        .prefix,
        b.pathJoin(&.{ "share/icons/hicolor/256x256/apps", icon_name }),
    ).step);
}
