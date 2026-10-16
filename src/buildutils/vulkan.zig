const rstdbuild = @import("__root__.zig");

pub fn addVulkan(b: *rstdbuild.Build) rstdbuild.DepModule {
    // get the path to the vulkan registry XML file
    const registry = b.dependency("vulkan_headers", .{}).path("registry/vk.xml");

    // generate vulkan.zig from the registry XML file using vulkan-zig-generator.
    // this returns a compile step for the vulkan-zig-generator exe.
    const gen = b.dependency("vulkan", .{}).artifact("vulkan-zig-generator");

    // add a run step to run the vulkan-zig-generator exe
    // The output file is added to the build graph as a LazyPath, so that the build system knows about it and can track its dependencies.
    const cmd = b.addRunArtifact(gen);

    // pass it the registry XML file and telling it to output vulkan.zig.
    cmd.addFileArg(registry);

    // return a Zig module build object for the output vk.zig as "vulkan"
    return .{ .mod = b.addModule("vulkan-zig", .{
        .root_source_file = cmd.addOutputFileArg2("vk.zig", .{}),
    }), .name = "vulkan" };
}

/// Adds rstd's own generated vulkan-zig module (built above, wired onto
/// rstd's shared module when `-Dvulkan=true`) onto `mod` under the import
/// name "vulkan". `rstd_mod` is the `DepModule` a consumer's own
/// `addRstd(...)` call returned.
///
/// A consumer that wants its own code to see the same `vk.*` types
/// rgl/vulkan/gpu.zig uses (so handles can cross that boundary) should
/// call this instead of generating a second, separate copy via `addVulkan`
/// -- vk-zig's types are nominal, so two independently-generated copies
/// are NOT interchangeable even when byte-identical.
pub fn addVulkanImportToModule(rstd_mod: rstdbuild.DepModule, mod: *rstdbuild.Module) void {
    const dep = rstd_mod.dep orelse
        @panic("addVulkanImportToModule needs rstd_mod.dep -- addRstd() must have fetched rstd via b.dependency(), not reused an in-graph module");
    mod.addImport("vulkan", dep.module("vulkan-zig"));
}
