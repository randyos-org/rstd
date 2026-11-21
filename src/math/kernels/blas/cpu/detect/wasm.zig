//! WebAssembly CPU selection: the module's compiled features decide; the
//! reference's one target is `WASM128_GENERIC` (`kernel/wasm/`).

const std = @import("std");
const coretypes = @import("../coretypes.zig");
const select = @import("../select.zig");
const sysinfo = @import("builtin");

pub const CoreType = coretypes.WasmCoreType;
pub const Decision = select.Decision(CoreType);

/// Nothing is probed; the module's features are the answer.
pub const Caps = struct {
    simd128: bool = sysinfo.target.cpu.arch.isWasm() and std.Target.wasm.featureSetHas(sysinfo.target.cpu.features, .simd128),
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

test "SIMD128 modules get the SIMD core" {
    try std.testing.expectEqual(CoreType.wasm128_generic, decide(.{ .simd128 = true }).core);
    try std.testing.expectEqual(CoreType.wasm_generic, decide(.{ .simd128 = false }).core);
    try std.testing.expect(!runnable(.wasm128_generic, .{ .simd128 = false }));
}

test "detect runs" {
    _ = select.detect(@This(), null);
}
