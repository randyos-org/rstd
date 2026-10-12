const std = @import("std");

pub const machine = @import("../machine.zig");
pub const io = @import("../io.zig");

pub const Build = std.Build;
pub const Step = Build.Step;
pub const BuildDir = Step.WriteFile;
pub const InstallDir = Step.InstallDir;
pub const Optimize = std.builtin.Optimize;
/// Kept for consumers still on the pre-0.18 name (the top-level randyos
/// build); new code uses `Optimize`.
pub const OptimizeMode = Optimize;
pub const TargetQuery = std.Target.Query;
pub const RunStep = Step.Run;
pub const CompileStep = Step.Compile;
pub const Module = Build.Module;
pub const BuildOptions = Step.Options;
pub const ResolvedTarget = Build.ResolvedTarget;
pub const LazyPath = Build.LazyPath;
pub const Dependency = Build.Dependency;

pub const spirv = @import("spirv.zig");
pub const wireRglShaders = spirv.wireRglShaders;
pub const addShader = spirv.addShader;
pub const ShaderStage = spirv.ShaderStage;

const builder_mod = @import("builder.zig");
pub const getOptimize = builder_mod.getOptimize;
pub const getTarget = builder_mod.getTarget;

const module_mod = @import("module.zig");
pub const cloneModuleWithTarget = module_mod.cloneModuleWithTarget;
pub const cloneRootForTarget = module_mod.cloneRootForTarget;

const files_mod = @import("files.zig");
pub const addBuildDir = files_mod.addBuildDir;
pub const installDirLazyPath = files_mod.installDirLazyPath;

const docs_mod = @import("docs.zig");
pub const Docs = docs_mod.Docs;
pub const addDocs = docs_mod.addDocs;
pub const addCompileStepDocsToStep = docs_mod.addCompileStepDocsToStep;

const options_mod = @import("options.zig");
pub const addBuildOptionsModule = options_mod.addBuildOptionsModule;
pub const LogLevel = options_mod.LogLevel;
pub const ScopeLevel = options_mod.ScopeLevel;
pub const addLogScopeOptions = options_mod.addLogScopeOptions;
pub const addBuildOptsModToModule = options_mod.addBuildOptsModToModule;

const deps_mod = @import("deps.zig");
pub const DepCLink = deps_mod.DepCLink;
pub const DepCBind = deps_mod.DepCBind;
pub const DepModule = deps_mod.DepModule;

const freedesktop_mod = @import("freedesktop.zig");
pub const addLinuxDesktopFile = freedesktop_mod.addLinuxDesktopFile;

const icon_mod = @import("icon.zig");
pub const Icon = icon_mod.Icon;

const rstd_dep_mod = @import("rstd.zig");
pub const addRstd = rstd_dep_mod.addRstd;
pub const RstdOptions = rstd_dep_mod.RstdOptions;

const win32_mod = @import("win32.zig");
pub const addWin32 = win32_mod.addWin32;
pub const addWin32ToModIfWinTgt = win32_mod.addWin32ToModIfWinTgt;

const renderdoc_mod = @import("renderdoc.zig");
pub const addRenderDoc = renderdoc_mod.addRenderDoc;

const vulkan_mod = @import("vulkan.zig");
pub const addVulkan = vulkan_mod.addVulkan;
pub const addVulkanImportToModule = vulkan_mod.addVulkanImportToModule;
