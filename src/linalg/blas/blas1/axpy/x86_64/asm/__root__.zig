//! Which converted file has which kernel, its cores in its `keys`; converted from OpenBLAS `v0.3.34-794-gdbf2c24c9` by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const registry = @import("../../../../registry.zig");

pub const Entry = registry.Entry;

pub const entries = [_]Entry{
    .{ .slot = .saxpy_k, .module = @import("saxpy_sse_atom_sysv/__root__.zig") },
    .{ .slot = .saxpy_k, .module = @import("saxpy_sse_dunnington_sysv/__root__.zig") },
    .{ .slot = .saxpy_k, .module = @import("saxpy_sse_atom_win/__root__.zig") },
    .{ .slot = .saxpy_k, .module = @import("saxpy_sse_dunnington_win/__root__.zig") },
    .{ .slot = .daxpy_k, .module = @import("daxpy_atom_sysv.zig") },
    .{ .slot = .daxpy_k, .module = @import("daxpy_atom_win.zig") },
    .{ .slot = .caxpy_k, .module = @import("caxpy_sse_atom_sysv/__root__.zig") },
    .{ .slot = .caxpy_k, .module = @import("caxpy_sse_dunnington_sysv/__root__.zig") },
    .{ .slot = .caxpy_k, .module = @import("caxpy_sse_atom_win/__root__.zig") },
    .{ .slot = .caxpy_k, .module = @import("caxpy_sse_dunnington_win/__root__.zig") },
    .{ .slot = .zaxpy_k, .module = @import("zaxpy_atom_sysv/__root__.zig") },
    .{ .slot = .zaxpy_k, .module = @import("zaxpy_atom_win/__root__.zig") },
    .{ .slot = .caxpyc_k, .module = @import("caxpyc_sse_atom_sysv/__root__.zig") },
    .{ .slot = .caxpyc_k, .module = @import("caxpyc_sse_dunnington_sysv/__root__.zig") },
    .{ .slot = .caxpyc_k, .module = @import("caxpyc_sse_atom_win/__root__.zig") },
    .{ .slot = .caxpyc_k, .module = @import("caxpyc_sse_dunnington_win/__root__.zig") },
    .{ .slot = .zaxpyc_k, .module = @import("zaxpyc_atom_sysv/__root__.zig") },
    .{ .slot = .zaxpyc_k, .module = @import("zaxpyc_atom_win/__root__.zig") },
    .{ .slot = .daxpy_k, .module = @import("daxpy_sse2_barcelona_sysv/__root__.zig") },
    .{ .slot = .daxpy_k, .module = @import("daxpy_sse2_dunnington_sysv.zig") },
    .{ .slot = .daxpy_k, .module = @import("daxpy_sse2_opteron_sysv.zig") },
    .{ .slot = .daxpy_k, .module = @import("daxpy_sse2_barcelona_win/__root__.zig") },
    .{ .slot = .daxpy_k, .module = @import("daxpy_sse2_dunnington_win/__root__.zig") },
    .{ .slot = .zaxpy_k, .module = @import("zaxpy_sse2_barcelona_sysv/__root__.zig") },
    .{ .slot = .zaxpy_k, .module = @import("zaxpy_sse2_opteron_sysv/__root__.zig") },
    .{ .slot = .zaxpy_k, .module = @import("zaxpy_sse2_barcelona_win/__root__.zig") },
    .{ .slot = .zaxpy_k, .module = @import("zaxpy_sse2_opteron_win/__root__.zig") },
    .{ .slot = .zaxpyc_k, .module = @import("zaxpyc_sse2_barcelona_sysv/__root__.zig") },
    .{ .slot = .zaxpyc_k, .module = @import("zaxpyc_sse2_opteron_sysv/__root__.zig") },
    .{ .slot = .zaxpyc_k, .module = @import("zaxpyc_sse2_barcelona_win/__root__.zig") },
    .{ .slot = .zaxpyc_k, .module = @import("zaxpyc_sse2_opteron_win/__root__.zig") },
};
