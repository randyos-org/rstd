//! Helpers for build scripts (rstd's own and its consumers'). Import only `std`, files here and `../bootinfo.zig`: a build script
//! recompiles on any edit to a file it can reach, so an edge into the rest of rstd rebuilds every consumer's build script.

const std = @import("std");

pub const bootinfo = @import("../bootinfo.zig");
pub const io = @import("io.zig");

pub const Build = std.Build;
pub const Step = Build.Step;
pub const BuildDir = Step.WriteFile;
pub const InstallDir = Step.InstallDir;
pub const Optimize = std.lang.Optimize;
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
pub const addShaderWith = spirv.addShaderWith;
pub const ShaderStage = spirv.ShaderStage;

const builder_file = @import("builder.zig");
pub const getOptimize = builder_file.getOptimize;
pub const getTarget = builder_file.getTarget;
pub const console = builder_file.console;

const module_file = @import("module.zig");
pub const cloneModuleWithTarget = module_file.cloneModuleWithTarget;
pub const cloneRootForTarget = module_file.cloneRootForTarget;

const files_file = @import("files.zig");
pub const addBuildDir = files_file.addBuildDir;
pub const installDirLazyPath = files_file.installDirLazyPath;

const docs_file = @import("docs.zig");
pub const Docs = docs_file.Docs;
pub const addDocs = docs_file.addDocs;
pub const addCompileStepDocsToStep = docs_file.addCompileStepDocsToStep;

const options_file = @import("options.zig");
pub const addBuildOptionsModule = options_file.addBuildOptionsModule;
pub const LogLevel = options_file.LogLevel;
pub const ScopeLevel = options_file.ScopeLevel;
pub const addLogScopeOptions = options_file.addLogScopeOptions;
pub const addBuildOptsModToModule = options_file.addBuildOptsModToModule;
pub const optionsModule = options_file.optionsModule;

const deps_file = @import("deps.zig");
pub const DepCLink = deps_file.DepCLink;
pub const DepCBind = deps_file.DepCBind;
pub const DepModule = deps_file.DepModule;
pub const declaresDependency = deps_file.declaresDependency;

const freedesktop_file = @import("freedesktop.zig");
pub const addLinuxDesktopFile = freedesktop_file.addLinuxDesktopFile;

const icon_file = @import("icon.zig");
pub const Icon = icon_file.Icon;

const rstd_file = @import("rstd.zig");
pub const addRstd = rstd_file.addRstd;
pub const RstdOptions = rstd_file.RstdOptions;
pub const BlasCores = rstd_file.BlasCores;

const win32_file = @import("win32.zig");
pub const addWin32 = win32_file.addWin32;
pub const addWin32ToModIfWinTgt = win32_file.addWin32ToModIfWinTgt;

const renderdoc_file = @import("renderdoc.zig");
pub const addRenderDoc = renderdoc_file.addRenderDoc;

const vulkan_file = @import("vulkan.zig");
pub const addVulkan = vulkan_file.addVulkan;
pub const addVulkanImportToModule = vulkan_file.addVulkanImportToModule;
