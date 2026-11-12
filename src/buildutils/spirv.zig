const std = @import("std");

const rstdbuild = @import("__root__.zig");
const Build = rstdbuild.Build;
const LazyPath = rstdbuild.LazyPath;
const Module = rstdbuild.Module;
const ResolvedTarget = rstdbuild.ResolvedTarget;
const Optimize = rstdbuild.Optimize;
const RstdOptions = rstdbuild.RstdOptions;

pub const ShaderStage = @import("shadertypes.zig").Stage;

pub const CompiledShader = struct {
    spv_vert_path: LazyPath,
    spv_frag_path: LazyPath,
    vert_mod_name: []const u8,
    frag_mod_name: []const u8,

    /// Exposes the compiled SPIR-V as `shaders.<name>.vert` /
    /// `shaders.<name>.frag` anonymous imports; consuming Zig code
    /// loads the raw bytes with `@embedFile` under those names.
    pub fn addToModule(self: *const CompiledShader, mod: *Module) void {
        mod.addAnonymousImport(self.vert_mod_name, .{ .root_source_file = self.spv_vert_path });
        mod.addAnonymousImport(self.frag_mod_name, .{ .root_source_file = self.spv_frag_path });
    }
};

fn compileShaderStage(
    b: *Build,
    spirv_target: ResolvedTarget,
    optimize: Optimize,
    src: LazyPath,
    name: []const u8,
    stage: ShaderStage,
    rstd_opts: RstdOptions,
    imports: []const Module.Import,
) !LazyPath {
    const shader_opts = b.addOptions();
    shader_opts.addOption(ShaderStage, "stage", stage);
    // shader_opts.addOption(bool, "add_rstd", add_rstd);

    const shader_mod = b.createModule(.{
        .root_source_file = src,
        .target = spirv_target,
        .optimize = optimize,
        .imports = imports,
    });
    shader_mod.addOptions("shader_opts", shader_opts);
    // // A standalone module giving shaders/lib/ access to specific rstd
    // // functions and types without needing the full `rstd` module. Rooted at
    // // _shader_imports.zig so those two modules' own cross-imports of each other
    // // stay inside the module boundary.
    // const shader_imports_mod = b.createModule(.{ .root_source_file = b.path("src/_shader_imports.zig") });
    // shader_mod.addImport("shader_imports", shader_imports_mod);
    // if (add_rstd) {
    const rstd_mod = try rstdbuild.addRstd(b, spirv_target, null, true, rstd_opts);
    rstd_mod.addImportToModule(shader_mod);
    // }

    const exe = b.addExecutable(.{
        .name = name,
        .root_module = shader_mod,
        // ReleaseFast defaults to trying LLVM for this target, which
        // crashes the compiler outright (SpirV.zig asserts !use_llvm).
        // Debug happens to default to the non-LLVM backend already, so
        // this only bites in optimized builds unless forced off here.
        .use_llvm = false,
    });
    return exe.getEmittedBin();
}

/// Compile `src` to a vert+frag SPIR-V pair named `name`.
pub fn addShader(
    b: *Build,
    optimize: Optimize,
    src: LazyPath,
    name: []const u8,
    rstd_opts: RstdOptions,
) !CompiledShader {
    return addShaderWith(b, optimize, src, name, rstd_opts, &.{});
}

/// `addShader`, with more modules for `src` to import.
pub fn addShaderWith(
    b: *Build,
    optimize: Optimize,
    src: LazyPath,
    name: []const u8,
    rstd_opts: RstdOptions,
    imports: []const Module.Import,
) !CompiledShader {
    const spirv_target = b.resolveTargetQuery(.{
        .cpu_arch = .spirv32,
        .cpu_model = .{ .explicit = &std.Target.spirv.cpu.vulkan_v1_2 },
        .os_tag = .vulkan,
    });
    const vert_name = b.fmt("{s}.vert.spv", .{name});
    const frag_name = b.fmt("{s}.frag.spv", .{name});
    return .{
        .spv_vert_path = try compileShaderStage(b, spirv_target, optimize, src, vert_name, .vert, rstd_opts, imports),
        .spv_frag_path = try compileShaderStage(b, spirv_target, optimize, src, frag_name, .frag, rstd_opts, imports),
        .vert_mod_name = vert_name,
        .frag_mod_name = frag_name,
    };
}

/// Compile rgl's three generic pipelines and wire `rgl/shaders/__root__.zig`'s
/// @embedFiles directly onto `rstd_mod`. Called from rstd's own
/// build.zig (gated behind `-Dvulkan`), not from a consumer's; only the
/// build.zig that created `rstd_mod` can add imports to it.
pub fn wireRglShaders(b: *Build, optimize: Optimize, rstd_mod: *Module) !void {
    const rgl = try addShader(b, optimize, b.path("src/rgl/shaders/shapes.zig"), "rgl.shapes", .{ ._rstd_mod = rstd_mod });
    const glyph = try addShader(b, optimize, b.path("src/rgl/shaders/glyph.zig"), "rgl.glyph", .{ ._rstd_mod = rstd_mod });
    const image = try addShader(b, optimize, b.path("src/rgl/shaders/image.zig"), "rgl.image", .{ ._rstd_mod = rstd_mod });
    rgl.addToModule(rstd_mod);
    glyph.addToModule(rstd_mod);
    image.addToModule(rstd_mod);
}
