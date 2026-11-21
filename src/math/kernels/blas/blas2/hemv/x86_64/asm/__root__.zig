//! Which converted file has which kernel, its cores in its `keys`; converted from OpenBLAS `v0.3.34-794-gdbf2c24c9` by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const registry = @import("../../../../registry.zig");

pub const Entry = registry.Entry;

pub const entries = [_]Entry{
    .{ .slot = .zhemv_U, .module = @import("zhemv_U_sse2_atom_sysv/__root__.zig") },
    .{ .slot = .zhemv_U, .module = @import("zhemv_U_sse2_core2_sysv/__root__.zig") },
    .{ .slot = .zhemv_U, .module = @import("zhemv_U_sse2_opteron_sysv/__root__.zig") },
    .{ .slot = .zhemv_U, .module = @import("zhemv_U_sse2_atom_win/__root__.zig") },
    .{ .slot = .zhemv_U, .module = @import("zhemv_U_sse2_core2_win/__root__.zig") },
    .{ .slot = .zhemv_U, .module = @import("zhemv_U_sse2_opteron_win/__root__.zig") },
    .{ .slot = .zhemv_L, .module = @import("zhemv_L_sse2_atom_sysv/__root__.zig") },
    .{ .slot = .zhemv_L, .module = @import("zhemv_L_sse2_core2_sysv/__root__.zig") },
    .{ .slot = .zhemv_L, .module = @import("zhemv_L_sse2_opteron_sysv/__root__.zig") },
    .{ .slot = .zhemv_L, .module = @import("zhemv_L_sse2_atom_win/__root__.zig") },
    .{ .slot = .zhemv_L, .module = @import("zhemv_L_sse2_core2_win/__root__.zig") },
    .{ .slot = .zhemv_L, .module = @import("zhemv_L_sse2_opteron_win/__root__.zig") },
};
