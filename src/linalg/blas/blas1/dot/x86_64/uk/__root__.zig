//! Which block kernel each core's C driver uses, read from the `#if`
//! chain at the top of `kernel/x86_64/<p>dot.c` (OpenBLAS `v0.3.34-317-ga69062524`). Generated
//! by `tools/blas_asm/uk.zig`; re-derive it rather than editing it.

const uk = @import("../../../uk_block.zig");

pub const Entry = uk.Entry;

pub const entries = [_]Entry{
    .{ .slot = .dsdot_k, .cores = &.{
        "bulldozer",
    }, .mask = 32, .block = .assembly, .kernel = "sdot_kernel_16", .module = @import("sdot_bulldozer.zig") },
    .{ .slot = .dsdot_k, .cores = &.{
        "excavator",
        "piledriver",
        "steamroller",
    }, .mask = 32, .block = .assembly, .kernel = "sdot_kernel_16", .module = @import("sdot_steamroller.zig") },
    .{ .slot = .dsdot_k, .cores = &.{
        "nehalem",
    }, .mask = 32, .block = .assembly, .kernel = "sdot_kernel_16", .module = @import("sdot_nehalem.zig") },
    .{ .slot = .dsdot_k, .cores = &.{
        "haswell",
        "zen",
    }, .mask = 32, .block = .assembly, .kernel = "sdot_kernel_16", .module = @import("sdot_haswell.zig") },
    .{ .slot = .dsdot_k, .cores = &.{
        "cooperlake",
        "sapphirerapids",
        "skylakex",
    }, .mask = 32, .block = .vector, .kernel = "", .module = void },
    .{ .slot = .dsdot_k, .cores = &.{
        "sandybridge",
    }, .mask = 32, .block = .assembly, .kernel = "sdot_kernel_16", .module = @import("sdot_sandy.zig") },
    .{ .slot = .ddot_k, .cores = &.{
        "excavator",
        "steamroller",
    }, .mask = 16, .block = .assembly, .kernel = "ddot_kernel_8", .module = @import("ddot_steamroller.zig") },
    .{ .slot = .ddot_k, .cores = &.{
        "piledriver",
    }, .mask = 16, .block = .assembly, .kernel = "ddot_kernel_8", .module = @import("ddot_piledriver.zig") },
    .{ .slot = .ddot_k, .cores = &.{
        "nehalem",
    }, .mask = 16, .block = .assembly, .kernel = "ddot_kernel_8", .module = @import("ddot_nehalem.zig") },
    .{ .slot = .ddot_k, .cores = &.{
        "haswell",
        "zen",
    }, .mask = 16, .block = .assembly, .kernel = "ddot_kernel_8", .module = @import("ddot_haswell.zig") },
    .{ .slot = .ddot_k, .cores = &.{
        "cooperlake",
        "sapphirerapids",
        "skylakex",
    }, .mask = 16, .block = .vector, .kernel = "", .module = void },
    .{ .slot = .ddot_k, .cores = &.{
        "sandybridge",
    }, .mask = 16, .block = .assembly, .kernel = "ddot_kernel_8", .module = @import("ddot_sandy.zig") },
    .{ .slot = .cdotu_k, .cores = &.{
        "bulldozer",
    }, .mask = 16, .block = .assembly, .kernel = "cdot_kernel_16", .module = @import("cdot_bulldozer.zig") },
    .{ .slot = .cdotc_k, .cores = &.{
        "bulldozer",
    }, .mask = 16, .block = .assembly, .kernel = "cdot_kernel_16", .module = @import("cdot_bulldozer.zig") },
    .{ .slot = .cdotu_k, .cores = &.{
        "excavator",
        "piledriver",
        "steamroller",
    }, .mask = 16, .block = .assembly, .kernel = "cdot_kernel_16", .module = @import("cdot_steamroller.zig") },
    .{ .slot = .cdotc_k, .cores = &.{
        "excavator",
        "piledriver",
        "steamroller",
    }, .mask = 16, .block = .assembly, .kernel = "cdot_kernel_16", .module = @import("cdot_steamroller.zig") },
    .{ .slot = .cdotu_k, .cores = &.{
        "cooperlake",
        "haswell",
        "sapphirerapids",
        "skylakex",
        "zen",
    }, .mask = 16, .block = .assembly, .kernel = "cdot_kernel_16", .module = @import("cdot_haswell.zig") },
    .{ .slot = .cdotc_k, .cores = &.{
        "cooperlake",
        "haswell",
        "sapphirerapids",
        "skylakex",
        "zen",
    }, .mask = 16, .block = .assembly, .kernel = "cdot_kernel_16", .module = @import("cdot_haswell.zig") },
    .{ .slot = .cdotu_k, .cores = &.{
        "sandybridge",
    }, .mask = 16, .block = .assembly, .kernel = "cdot_kernel_16", .module = @import("cdot_sandy.zig") },
    .{ .slot = .cdotc_k, .cores = &.{
        "sandybridge",
    }, .mask = 16, .block = .assembly, .kernel = "cdot_kernel_16", .module = @import("cdot_sandy.zig") },
    .{ .slot = .zdotu_k, .cores = &.{
        "bulldozer",
    }, .mask = 8, .block = .assembly, .kernel = "zdot_kernel_8", .module = @import("zdot_bulldozer.zig") },
    .{ .slot = .zdotc_k, .cores = &.{
        "bulldozer",
    }, .mask = 8, .block = .assembly, .kernel = "zdot_kernel_8", .module = @import("zdot_bulldozer.zig") },
    .{ .slot = .zdotu_k, .cores = &.{
        "excavator",
        "piledriver",
        "steamroller",
    }, .mask = 8, .block = .assembly, .kernel = "zdot_kernel_8", .module = @import("zdot_steamroller.zig") },
    .{ .slot = .zdotc_k, .cores = &.{
        "excavator",
        "piledriver",
        "steamroller",
    }, .mask = 8, .block = .assembly, .kernel = "zdot_kernel_8", .module = @import("zdot_steamroller.zig") },
    .{ .slot = .zdotu_k, .cores = &.{
        "cooperlake",
        "haswell",
        "sapphirerapids",
        "skylakex",
        "zen",
    }, .mask = 8, .block = .assembly, .kernel = "zdot_kernel_8", .module = @import("zdot_haswell.zig") },
    .{ .slot = .zdotc_k, .cores = &.{
        "cooperlake",
        "haswell",
        "sapphirerapids",
        "skylakex",
        "zen",
    }, .mask = 8, .block = .assembly, .kernel = "zdot_kernel_8", .module = @import("zdot_haswell.zig") },
    .{ .slot = .zdotu_k, .cores = &.{
        "sandybridge",
    }, .mask = 8, .block = .assembly, .kernel = "zdot_kernel_8", .module = @import("zdot_sandy.zig") },
    .{ .slot = .zdotc_k, .cores = &.{
        "sandybridge",
    }, .mask = 8, .block = .assembly, .kernel = "zdot_kernel_8", .module = @import("zdot_sandy.zig") },
};

/// The entry for `slot` on `core`, if its driver is C.
pub fn find(comptime slot: uk.Slot, comptime core: []const u8) ?Entry {
    return uk.find(&entries, slot, core);
}
