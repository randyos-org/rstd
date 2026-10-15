//! WebAssembly CPU selection. There is nothing to detect at run time: which
//! vector instructions the module may use is fixed when it is compiled, so the
//! core type follows the target's features. The reference has only the one
//! target for it, `WASM128_GENERIC`, in `kernel/wasm/`.

const std = @import("std");
const coretypes = @import("../_coretypes.zig");
const select = @import("../_select.zig");
const sysinfo = @import("builtin");

pub const CoreType = coretypes.WasmCoreType;
pub const Decision = select.Decision(CoreType);

/// Nothing is probed; the module's features are the answer.
pub const Caps = struct {
    simd128: bool = sysinfo.cpu.arch.isWasm() and std.Target.wasm.featureSetHas(sysinfo.cpu.features, .simd128),
};

pub fn probe() Caps {
    return .{};
}

pub fn runnable(core: CoreType, c: Caps) bool {
    return core == .wasm_generic or c.simd128;
}

pub fn decide(c: Caps) Decision {
    return .{ .core = if (c.simd128) .wasm128_generic else .wasm_generic };
}

pub fn detect() CoreType {
    return decide(probe()).core;
}

pub fn detectNoted() Decision {
    return decide(probe());
}

test "SIMD128 modules get the SIMD core" {
    try std.testing.expectEqual(CoreType.wasm128_generic, decide(.{ .simd128 = true }).core);
    try std.testing.expectEqual(CoreType.wasm_generic, decide(.{ .simd128 = false }).core);
    try std.testing.expect(!runnable(.wasm128_generic, .{ .simd128 = false }));
}

test "detect runs" {
    _ = detect();
}
