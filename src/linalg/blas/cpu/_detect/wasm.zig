//! WebAssembly CPU selection. There is nothing to detect at run time: which
//! vector instructions the module may use is fixed when it is compiled, so the
//! core type follows the target's features. The reference has only the one
//! target for it, `WASM128_GENERIC`, in `kernel/wasm/`.

const std = @import("std");
const coretypes = @import("../_coretypes.zig");
const sysinfo = @import("builtin");

pub const CoreType = coretypes.WasmCoreType;

pub fn detect() CoreType {
    if (comptime sysinfo.cpu.arch.isWasm() and std.Target.wasm.featureSetHas(sysinfo.cpu.features, .simd128)) return .wasm128_generic;
    return .wasm_generic;
}

test "detect runs" {
    _ = detect();
}
