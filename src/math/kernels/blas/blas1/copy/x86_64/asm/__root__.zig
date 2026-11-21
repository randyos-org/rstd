//! Which converted file has which kernel, its cores in its `keys`; converted from OpenBLAS `v0.3.34-794-gdbf2c24c9` by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const registry = @import("../../../../registry.zig");

pub const Entry = registry.Entry;

pub const entries = [_]Entry{
    .{ .slot = .scopy_k, .module = @import("scopy_sse_sysv/__root__.zig") },
    .{ .slot = .scopy_k, .module = @import("scopy_sse_win/__root__.zig") },
    .{ .slot = .dcopy_k, .module = @import("dcopy_sse2_atom_sysv.zig") },
    .{ .slot = .dcopy_k, .module = @import("dcopy_sse2_cooperlake_sysv.zig") },
    .{ .slot = .dcopy_k, .module = @import("dcopy_sse2_atom_win.zig") },
    .{ .slot = .dcopy_k, .module = @import("dcopy_sse2_cooperlake_win.zig") },
    .{ .slot = .ccopy_k, .module = @import("ccopy_sse_sysv/__root__.zig") },
    .{ .slot = .ccopy_k, .module = @import("ccopy_sse_win/__root__.zig") },
    .{ .slot = .zcopy_k, .module = @import("zcopy_sse2_atom_sysv.zig") },
    .{ .slot = .zcopy_k, .module = @import("zcopy_sse2_cooperlake_sysv.zig") },
    .{ .slot = .zcopy_k, .module = @import("zcopy_sse2_atom_win.zig") },
    .{ .slot = .zcopy_k, .module = @import("zcopy_sse2_cooperlake_win.zig") },
    .{ .slot = .dcopy_k, .module = @import("dcopy_bulldozer_sysv.zig") },
    .{ .slot = .dcopy_k, .module = @import("dcopy_bulldozer_win.zig") },
};
