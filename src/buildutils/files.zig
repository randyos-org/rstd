const rstdbuild = @import("__root__.zig");
const Build = rstdbuild.Build;
const LazyPath = rstdbuild.LazyPath;
const InstallDir = rstdbuild.InstallDir;
const Step = rstdbuild.Step;

pub fn addBuildDir(b: *Build) *Step.WriteFile {
    return b.addWriteFiles();
}

pub fn installDirLazyPath(b: *Build, instdir: *InstallDir) LazyPath {
    const dest_rel_path = instdir.options.install_subdir;
    return switch (instdir.options.install_dir) {
        .prefix => .{ .relative = .{ .base = .install_prefix, .sub_path = dest_rel_path } },
        .lib => .{ .relative = .{ .base = .install_lib, .sub_path = dest_rel_path } },
        .bin => .{ .relative = .{ .base = .install_bin, .sub_path = dest_rel_path } },
        .header => .{ .relative = .{ .base = .install_include, .sub_path = dest_rel_path } },
        .custom => |p| .{ .relative = .{
            .base = .install_prefix,
            .sub_path = if (dest_rel_path.len == 0) p else b.pathJoin(&.{ p, dest_rel_path }),
        } },
    };
}
