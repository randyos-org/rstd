// const std = @import("std");
const rstdbuild = @import("__root__.zig");
const Build = rstdbuild.Build;
const LazyPath = rstdbuild.LazyPath;
const OptimizeMode = rstdbuild.OptimizeMode;
const Module = rstdbuild.Module;

const files = @import("files.zig");
const addRstd = @import("rstd.zig").addRstd;
// const rc = @import("rc");

pub const Icon = struct {
    const Self = @This();
    rc_build_path: LazyPath,
    png_build_path: LazyPath,
    png_name: []const u8,
    rgba_build_path: LazyPath,
    rgba_name: []const u8,

    // `addWin32ResourceFile` is deprecated in favor of the `rc` package
    // (https://codeberg.org/ziglang/rc, already in build.zig.zon), but `rc`'s
    // resinator/aro dependency chain doesn't build yet against this
    // toolchain's std-library churn
    pub fn addToModule(self: *const Self, b: *Build, mod: *Module) void {
        _ = b;
        // const target = mod.resolved_target.?;
        // if (target.result.os.tag == .windows) {
        //     const rc_dep = b.dependency("rc", .{});
        //     var compiler: rc.Compiler = .init(rc_dep, .{
        //         .rc_files = &.{self.rc_build_path},
        //         .target = target,
        //         .optimize = mod.optimize.?,
        //     });
        //     mod.addObjectFile(compiler.getObjectFile());
        // }

        mod.addWin32ResourceFile(.{ .file = self.rc_build_path });
        mod.addAnonymousImport(self.png_name, .{ .root_source_file = self.png_build_path });
        mod.addAnonymousImport(self.rgba_name, .{ .root_source_file = self.rgba_build_path });
    }

    pub fn addGenIconStep(
        b: *Build,
        optimize: OptimizeMode,
        in_png_path: LazyPath,
        png_name: []const u8,
        ico_name: []const u8,
        rc_name: []const u8,
        rgba_name: []const u8,
    ) !Icon {
        const target = b.graph.host;
        const rstd = try addRstd(b, target, null, false);

        const mod = b.createModule(.{
            .root_source_file = rstd.dep.?.path("src/rgl/gen_icon.zig"),
            .target = target,
            .optimize = optimize,
        });
        rstd.addImportToModule(mod);

        const exe = b.addExecutable(.{
            .name = "genicon",
            .root_module = mod,
        });
        const run = b.addRunArtifact(exe);
        run.addFileArg(in_png_path);

        // record the paths of the generated output files so the build graph knows about them
        const png_build_path = run.addOutputFileArg(png_name);
        const ico_build_path = run.addOutputFileArg(ico_name);
        const rgba_build_path = run.addOutputFileArg(rgba_name);

        // The .rc ICON directive references "icon.ico" by a bare name, resolved
        // relative to the .rc file's own directory, so the generated icon.ico
        // needs to sit next to a (likewise generated) icon.rc, not just exist
        // somewhere else in the cache under its own path.
        const rc_build_dir = files.addBuildDir(b);
        const rc_build_path = rc_build_dir.add(rc_name, b.fmt("1 ICON \"{s}\"\n", .{ico_name}));
        _ = rc_build_dir.addCopyFile(ico_build_path, ico_name);

        // `zig build icon` - not needed for the app build itself (that already
        // depends on genicon's outputs above), but installs copies to zig-out/
        // so a human can open them after tweaking the drawing code.
        const icon_step = b.step("icon", "Regenerate icon.png/icon.ico into zig-out/ for inspection");
        icon_step.dependOn(&b.addInstallFileWithDir(png_build_path, .prefix, png_name).step);
        icon_step.dependOn(&b.addInstallFileWithDir(ico_build_path, .prefix, ico_name).step);
        icon_step.dependOn(&b.addInstallFileWithDir(rgba_build_path, .prefix, rgba_name).step);

        return .{
            .rc_build_path = rc_build_path,
            .png_build_path = png_build_path,
            .png_name = png_name,
            .rgba_build_path = rgba_build_path,
            .rgba_name = rgba_name,
        };
    }
};
