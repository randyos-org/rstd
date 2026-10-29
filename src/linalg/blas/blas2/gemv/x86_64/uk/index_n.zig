//! Which block kernel each core's C driver uses, read from the `#if`
//! chain at the top of `kernel/x86_64/<p>gemv_n_4.c` (OpenBLAS `v0.3.34-317-ga69062524`). Generated
//! by `tools/blas_asm/uk.zig`; re-derive it rather than editing it.

const uk = @import("../../../../blas1/uk_block.zig");

pub const Entry = uk.Entry;

pub const entries = [_]Entry{
    .{ .slot = .sgemv_n, .cores = &.{
        "bulldozer",
        "excavator",
        "piledriver",
        "steamroller",
    }, .mask = 0, .block = .assembly, .kernel = "sgemv_kernel_4x8", .module = @import("sgemv_n_bulldozer.zig"), .own = @import("sgemv_n_4.zig") },
    .{ .slot = .sgemv_n, .cores = &.{
        "nehalem",
    }, .mask = 0, .block = .assembly, .kernel = "sgemv_kernel_4x8", .module = @import("sgemv_n_nehalem.zig"), .own = @import("sgemv_n_4.zig") },
    .{ .slot = .sgemv_n, .cores = &.{
        "sandybridge",
    }, .mask = 0, .block = .assembly, .kernel = "sgemv_kernel_4x8", .module = @import("sgemv_n_sandy.zig"), .own = @import("sgemv_n_4.zig") },
    .{ .slot = .sgemv_n, .cores = &.{
        "haswell",
        "zen",
    }, .mask = 0, .block = .assembly, .kernel = "sgemv_kernel_4x8", .module = @import("sgemv_n_haswell.zig"), .own = @import("sgemv_n_4.zig") },
    .{ .slot = .sgemv_n, .cores = &.{
        "cooperlake",
        "sapphirerapids",
        "skylakex",
    }, .mask = 0, .block = .assembly, .kernel = "sgemv_kernel_4x8", .module = @import("sgemv_n_haswell.zig"), .own = @import("sgemv_n_4.zig") },
    .{ .slot = .sgemv_n, .cores = &.{
        "cooperlake",
        "sapphirerapids",
        "skylakex",
    }, .mask = 0, .block = .vector, .kernel = "", .module = void, .own = @import("sgemv_n_4.zig") },
    .{ .slot = .dgemv_n, .cores = &.{
        "nehalem",
    }, .mask = 0, .block = .assembly, .kernel = "dgemv_kernel_4x4", .module = @import("dgemv_n_nehalem.zig"), .own = @import("dgemv_n_4.zig") },
    .{ .slot = .dgemv_n, .cores = &.{
        "excavator",
        "haswell",
        "steamroller",
        "zen",
    }, .mask = 0, .block = .assembly, .kernel = "dgemv_kernel_4x4", .module = @import("dgemv_n_haswell.zig"), .own = @import("dgemv_n_4.zig") },
    .{ .slot = .dgemv_n, .cores = &.{
        "cooperlake",
        "sapphirerapids",
        "skylakex",
    }, .mask = 0, .block = .vector, .kernel = "", .module = void, .own = @import("dgemv_n_4.zig") },
    .{ .slot = .cgemv_n, .cores = &.{
        "cooperlake",
        "haswell",
        "sapphirerapids",
        "skylakex",
        "zen",
    }, .mask = 0, .block = .assembly, .kernel = "cgemv_kernel_4x4", .module = @import("cgemv_n_haswell.zig") },
    .{ .slot = .cgemv_r, .cores = &.{
        "cooperlake",
        "haswell",
        "sapphirerapids",
        "skylakex",
        "zen",
    }, .mask = 0, .block = .assembly, .kernel = "cgemv_kernel_4x4", .module = @import("cgemv_n_haswell.zig") },
    .{ .slot = .cgemv_o, .cores = &.{
        "cooperlake",
        "haswell",
        "sapphirerapids",
        "skylakex",
        "zen",
    }, .mask = 0, .block = .assembly, .kernel = "cgemv_kernel_4x4", .module = @import("cgemv_n_haswell.zig") },
    .{ .slot = .cgemv_s, .cores = &.{
        "cooperlake",
        "haswell",
        "sapphirerapids",
        "skylakex",
        "zen",
    }, .mask = 0, .block = .assembly, .kernel = "cgemv_kernel_4x4", .module = @import("cgemv_n_haswell.zig") },
    .{ .slot = .cgemv_n, .cores = &.{
        "bulldozer",
        "excavator",
        "piledriver",
        "steamroller",
    }, .mask = 0, .block = .assembly, .kernel = "cgemv_kernel_4x4", .module = @import("cgemv_n_bulldozer.zig") },
    .{ .slot = .cgemv_r, .cores = &.{
        "bulldozer",
        "excavator",
        "piledriver",
        "steamroller",
    }, .mask = 0, .block = .assembly, .kernel = "cgemv_kernel_4x4", .module = @import("cgemv_n_bulldozer.zig") },
    .{ .slot = .cgemv_o, .cores = &.{
        "bulldozer",
        "excavator",
        "piledriver",
        "steamroller",
    }, .mask = 0, .block = .assembly, .kernel = "cgemv_kernel_4x4", .module = @import("cgemv_n_bulldozer.zig") },
    .{ .slot = .cgemv_s, .cores = &.{
        "bulldozer",
        "excavator",
        "piledriver",
        "steamroller",
    }, .mask = 0, .block = .assembly, .kernel = "cgemv_kernel_4x4", .module = @import("cgemv_n_bulldozer.zig") },
    .{ .slot = .cgemv_n, .cores = &.{
        "atom",
        "barcelona",
        "bobcat",
        "core2",
        "dunnington",
        "nano",
        "nehalem",
        "opteron",
        "opteron_sse3",
        "penryn",
        "prescott",
        "sandybridge",
    }, .mask = 0, .block = .none, .kernel = "", .module = void },
    .{ .slot = .cgemv_r, .cores = &.{
        "atom",
        "barcelona",
        "bobcat",
        "core2",
        "dunnington",
        "nano",
        "nehalem",
        "opteron",
        "opteron_sse3",
        "penryn",
        "prescott",
        "sandybridge",
    }, .mask = 0, .block = .none, .kernel = "", .module = void },
    .{ .slot = .cgemv_o, .cores = &.{
        "atom",
        "barcelona",
        "bobcat",
        "core2",
        "dunnington",
        "nano",
        "nehalem",
        "opteron",
        "opteron_sse3",
        "penryn",
        "prescott",
        "sandybridge",
    }, .mask = 0, .block = .none, .kernel = "", .module = void },
    .{ .slot = .cgemv_s, .cores = &.{
        "atom",
        "barcelona",
        "bobcat",
        "core2",
        "dunnington",
        "nano",
        "nehalem",
        "opteron",
        "opteron_sse3",
        "penryn",
        "prescott",
        "sandybridge",
    }, .mask = 0, .block = .none, .kernel = "", .module = void },
    .{ .slot = .zgemv_n, .cores = &.{
        "cooperlake",
        "haswell",
        "sapphirerapids",
        "skylakex",
        "zen",
    }, .mask = 0, .block = .assembly, .kernel = "zgemv_kernel_4x4", .module = @import("zgemv_n_haswell.zig") },
    .{ .slot = .zgemv_r, .cores = &.{
        "cooperlake",
        "haswell",
        "sapphirerapids",
        "skylakex",
        "zen",
    }, .mask = 0, .block = .assembly, .kernel = "zgemv_kernel_4x4", .module = @import("zgemv_n_haswell.zig") },
    .{ .slot = .zgemv_o, .cores = &.{
        "cooperlake",
        "haswell",
        "sapphirerapids",
        "skylakex",
        "zen",
    }, .mask = 0, .block = .assembly, .kernel = "zgemv_kernel_4x4", .module = @import("zgemv_n_haswell.zig") },
    .{ .slot = .zgemv_s, .cores = &.{
        "cooperlake",
        "haswell",
        "sapphirerapids",
        "skylakex",
        "zen",
    }, .mask = 0, .block = .assembly, .kernel = "zgemv_kernel_4x4", .module = @import("zgemv_n_haswell.zig") },
    .{ .slot = .zgemv_n, .cores = &.{
        "sandybridge",
    }, .mask = 0, .block = .assembly, .kernel = "zgemv_kernel_4x4", .module = @import("zgemv_n_sandy.zig") },
    .{ .slot = .zgemv_r, .cores = &.{
        "sandybridge",
    }, .mask = 0, .block = .assembly, .kernel = "zgemv_kernel_4x4", .module = @import("zgemv_n_sandy.zig") },
    .{ .slot = .zgemv_o, .cores = &.{
        "sandybridge",
    }, .mask = 0, .block = .assembly, .kernel = "zgemv_kernel_4x4", .module = @import("zgemv_n_sandy.zig") },
    .{ .slot = .zgemv_s, .cores = &.{
        "sandybridge",
    }, .mask = 0, .block = .assembly, .kernel = "zgemv_kernel_4x4", .module = @import("zgemv_n_sandy.zig") },
    .{ .slot = .zgemv_n, .cores = &.{
        "bulldozer",
        "excavator",
        "piledriver",
        "steamroller",
    }, .mask = 0, .block = .assembly, .kernel = "zgemv_kernel_4x4", .module = @import("zgemv_n_bulldozer.zig") },
    .{ .slot = .zgemv_r, .cores = &.{
        "bulldozer",
        "excavator",
        "piledriver",
        "steamroller",
    }, .mask = 0, .block = .assembly, .kernel = "zgemv_kernel_4x4", .module = @import("zgemv_n_bulldozer.zig") },
    .{ .slot = .zgemv_o, .cores = &.{
        "bulldozer",
        "excavator",
        "piledriver",
        "steamroller",
    }, .mask = 0, .block = .assembly, .kernel = "zgemv_kernel_4x4", .module = @import("zgemv_n_bulldozer.zig") },
    .{ .slot = .zgemv_s, .cores = &.{
        "bulldozer",
        "excavator",
        "piledriver",
        "steamroller",
    }, .mask = 0, .block = .assembly, .kernel = "zgemv_kernel_4x4", .module = @import("zgemv_n_bulldozer.zig") },
    .{ .slot = .zgemv_n, .cores = &.{
        "barcelona",
        "core2",
        "dunnington",
        "nano",
        "nehalem",
        "opteron",
        "penryn",
    }, .mask = 0, .block = .none, .kernel = "", .module = void },
    .{ .slot = .zgemv_r, .cores = &.{
        "barcelona",
        "core2",
        "dunnington",
        "nano",
        "nehalem",
        "opteron",
        "penryn",
    }, .mask = 0, .block = .none, .kernel = "", .module = void },
    .{ .slot = .zgemv_o, .cores = &.{
        "barcelona",
        "core2",
        "dunnington",
        "nano",
        "nehalem",
        "opteron",
        "penryn",
    }, .mask = 0, .block = .none, .kernel = "", .module = void },
    .{ .slot = .zgemv_s, .cores = &.{
        "barcelona",
        "core2",
        "dunnington",
        "nano",
        "nehalem",
        "opteron",
        "penryn",
    }, .mask = 0, .block = .none, .kernel = "", .module = void },
};

/// The entry for `slot` on `core`, if its driver is C.
pub fn find(comptime slot: uk.Slot, comptime core: []const u8) ?Entry {
    return uk.find(&entries, slot, core);
}
