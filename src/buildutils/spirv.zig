const std = @import("std");

const rstdbuild = @import("__root__.zig");
const Build = rstdbuild.Build;
const LazyPath = rstdbuild.LazyPath;
const Module = rstdbuild.Module;
const ResolvedTarget = rstdbuild.ResolvedTarget;
const OptimizeMode = rstdbuild.OptimizeMode;

const shaderlib = @import("../rgl/shaders/lib/__root__.zig");
pub const ShaderStage = shaderlib.types.ShaderStage;

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
    optimize: OptimizeMode,
    src: LazyPath,
    name: []const u8,
    stage: ShaderStage,
    add_rstd: bool,
) !LazyPath {
    const shader_opts = b.addOptions();
    shader_opts.addOption(ShaderStage, "stage", stage);

    const shader_mod = b.createModule(.{
        .root_source_file = src,
        .target = spirv_target,
        .optimize = optimize,
    });
    shader_mod.addOptions("shader_opts", shader_opts);
    if (add_rstd) {
        const rstd_mod = try rstdbuild.addRstd(b, spirv_target, null, true);
        rstd_mod.addImportToModule(shader_mod);
    }

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
    optimize: OptimizeMode,
    src: LazyPath,
    name: []const u8,
    add_rstd: bool,
) !CompiledShader {
    const spirv_target = b.resolveTargetQuery(.{
        .cpu_arch = .spirv32,
        .cpu_model = .{ .explicit = &std.Target.spirv.cpu.vulkan_v1_2 },
        .os_tag = .vulkan,
    });
    const vert_name = b.fmt("{s}.vert.spv", .{name});
    const frag_name = b.fmt("{s}.frag.spv", .{name});
    return .{
        .spv_vert_path = try compileShaderStage(b, spirv_target, optimize, src, vert_name, .vert, add_rstd),
        .spv_frag_path = try compileShaderStage(b, spirv_target, optimize, src, frag_name, .frag, add_rstd),
        .vert_mod_name = vert_name,
        .frag_mod_name = frag_name,
    };
}

/// Compile rgl's three generic pipelines and wire `rgl/shaders/__root__.zig`'s
/// @embedFiles directly onto `rstd_mod`. Called from rstd's own
/// build.zig (gated behind `-Dvulkan`), not from a consumer's; only the
/// build.zig that created `rstd_mod` can add imports to it.
pub fn wireRglShaders(b: *Build, optimize: OptimizeMode, rstd_mod: *Module) !void {
    const rgl = try addShader(b, optimize, b.path("src/rgl/shaders/rgl.zig"), "rgl", false);
    const glyph = try addShader(b, optimize, b.path("src/rgl/shaders/glyph.zig"), "rgl.glyph", false);
    const image = try addShader(b, optimize, b.path("src/rgl/shaders/image.zig"), "rgl.image", false);
    rgl.addToModule(rstd_mod);
    glyph.addToModule(rstd_mod);
    image.addToModule(rstd_mod);
}
