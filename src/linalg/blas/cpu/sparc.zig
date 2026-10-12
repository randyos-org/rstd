//! SPARC CPU selection. The reference has two static targets and detects
//! nothing at run time, so the core type follows the target: V9 (any 64-bit
//! target, or the 32-bit one built for V9) or the generic one.

const std = @import("std");
const coretypes = @import("_coretypes.zig");
const sysinfo = @import("builtin");

pub const CoreType = coretypes.SparcCoreType;
pub const kernel = @import("../kernels/sparc.zig");

pub fn detect() CoreType {
    if (comptime sysinfo.cpu.arch.isSPARC() and (sysinfo.cpu.arch == .sparc64 or sysinfo.cpu.has(.sparc, .v9))) return .sparc;
    return .sparcv7;
}

const Vtable = @import("../Vtable.zig");
const blas1 = @import("../blas1/__root__.zig");
const blas = @import("../__root__.zig");
const gemm_params = @import("../blas3/gemm/params.zig");
pub fn vtableFromCoreType(core_type: CoreType) Vtable {
    var vt: Vtable = .{};
    blas1.installGeneric(&vt);
    kernel.install(&vt, core_type);
    gemm_params.installCore(&vt, .sparc, core_type, .fromConfig(blas.config));
    return vt;
}

test "detect runs" {
    _ = detect();
}
