//! Root of a per-core BLAS kernel object (see `build/blaskern.zig`).
//!
//! The build compiles this file once per group of cores that share compile
//! flags, with that group's CPU features, and links the result into anything
//! that uses `rstd`. `blaskernel_options.cores` names the cores whose kernels
//! this object exports; `fuse_fma_set`/`fuse_fma` carry the group's
//! contraction policy where the reference's flags fix one (see `fma.zig`).
//!
//! It is one root for every architecture: the registry it exports through is
//! chosen by the target (`kernels/__root__.zig`).

const std = @import("std");
const sysinfo = @import("builtin");
const opts = @import("blaskernel_options");
const fma = @import("fma.zig");

const kernels = @import("kernels/__root__.zig").ForArch(sysinfo.target.cpu.arch);

/// Read by `fma.zig`.
pub const blas_fuse_fma: bool = if (opts.fuse_fma_set) opts.fuse_fma else fma.hardware;

/// Read by `fma.zig`: a group that fuses on MIPS before release 6 can only
/// mean `madd.fmt` does (the Loongson 3A cores).
pub const blas_fused_madd: bool = opts.fuse_fma_set and opts.fuse_fma;

/// The vector multiply-adds'; the group's policy is for the scalar flags.
pub const blas_fuse_fma_vector: bool = if (opts.fuse_fma_set) opts.fuse_fma else fma.vector_madd;

comptime {
    @setEvalBranchQuota(1_000_000);
    for (opts.cores) |name| kernels.exportCore(std.meta.stringToEnum(kernels.Core, name).?);
}
