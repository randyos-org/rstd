//! Which block kernel each core's C driver uses, read from the `#if`
//! chain at the top of `kernel/x86_64/<p>rot.c` (OpenBLAS `v0.3.34-317-ga69062524`). Generated
//! by `tools/blas_asm/uk.zig`; re-derive it rather than editing it.

const uk = @import("../../../uk_block.zig");

pub const Entry = uk.Entry;

pub const entries = [_]Entry{
    .{ .slot = .srot_k, .cores = &.{
        "skylakex",
    }, .mask = 0, .block = .vector, .kernel = "", .module = void },
    .{ .slot = .srot_k, .cores = &.{
        "haswell",
        "zen",
    }, .mask = 0, .block = .vector, .kernel = "", .module = void },
    .{ .slot = .srot_k, .cores = &.{
        "cooperlake",
        "sapphirerapids",
    }, .mask = 0, .block = .none, .kernel = "", .module = void },
    .{ .slot = .drot_k, .cores = &.{
        "skylakex",
    }, .mask = 0, .block = .vector, .kernel = "", .module = void },
    .{ .slot = .drot_k, .cores = &.{
        "haswell",
        "zen",
    }, .mask = 0, .block = .vector, .kernel = "", .module = void },
    .{ .slot = .drot_k, .cores = &.{
        "cooperlake",
        "sapphirerapids",
    }, .mask = 0, .block = .none, .kernel = "", .module = void },
};

/// The entry for `slot` on `core`, if its driver is C.
pub fn find(comptime slot: uk.Slot, comptime core: []const u8) ?Entry {
    return uk.find(&entries, slot, core);
}
