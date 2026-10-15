//! Derived from code in OpenBLAS

const std = @import("std");
const zigstd = @import("../../zigstd.zig");

pub const types = @import("types.zig");
pub const cpu = @import("cpu/__root__.zig");
pub const Vtable = @import("Vtable.zig");
pub const Config = @import("Config.zig");
const blas1 = @import("blas1/__root__.zig");
const blas2 = @import("blas2/__root__.zig");
const blas3 = @import("blas3/__root__.zig");
const ext = @import("ext/__root__.zig");
// pub const util = @import("util.zig");

const InitPhase = enum(u8) { uninit, running, ready };

pub var config: Config = .{};
pub var core_type: cpu.CoreType = undefined;
pub var _vtable: Vtable = generic_vtable;

/// The vtable before `init` and after `shutdown`: the portable kernels, where
/// there are any.
const generic_vtable: Vtable = blk: {
    @setEvalBranchQuota(100_000);
    var vt: Vtable = .{};
    blas1.installGeneric(&vt);
    blas2.installGeneric(&vt);
    @import("portable.zig").alias(&vt);
    break :blk vt;
};

/// A vtable slot's name; too many fields for the default branch quota.
pub const Slot = blk: {
    @setEvalBranchQuota(10_000);
    break :blk std.meta.FieldEnum(Vtable);
};

/// The kernel in vtable slot `slot` (`.daxpy_k`), as the reference's `DAXPY_K`
/// is `gotoblas->daxpy_k`. A slot with no kernel for the running core, and no
/// portable one, is a bug in the caller.
pub fn kernel(comptime slot: Slot) @typeInfo(@FieldType(Vtable, @tagName(slot))).optional.child {
    return @field(_vtable, @tagName(slot)).?;
}

var phase = std.atomic.Value(InitPhase).init(.uninit);

/// Runs initialization exactly once.
///
/// A mutex would drag in threading support we do not want on freestanding
/// targets, so this is a three-state handshake: the winner of the
/// compare-exchange initializes while any other thread spins until it publishes.
pub fn init(env: *std.process.Environ.Map) void {
    if (phase.load(.acquire) == .ready) return;
    if (phase.cmpxchgStrong(.uninit, .running, .acquire, .monotonic) == null) {
        config = .init(env);

        if (config.core_type) |corename| {
            if (cpu.fromName(corename)) |ct| {
                core_type = ct;
            } else {
                // OpenBLAS warns and detects (`OPENBLAS_CORETYPE` names no core here)
                blasWarn(1, "unknown OPENBLAS_CORETYPE '{s}', detecting the core", .{corename});
                core_type = cpu.detect();
            }
        } else {
            core_type = cpu.detect();
        }
        blasWarn(2, "core: {s}", .{@tagName(core_type)});

        _vtable = cpu.vtableFromCoreType(core_type);
        if (_vtable.init) |f| f();
        phase.store(.ready, .release);
        return;
    }
    while (phase.load(.acquire) != .ready) std.atomic.spinLoopHint();
}

/// A diagnostic, shown when `OPENBLAS_VERBOSE` is at least `level`. Only where
/// there is somewhere to print to.
fn blasWarn(level: i32, comptime fmt: []const u8, args: anytype) void {
    if (comptime @import("builtin").os.tag == .freestanding) return;
    if (config.verbose < level) return;
    std.debug.print("OpenBLAS: " ++ fmt ++ "\n", args);
}

pub fn shutdown() void {
    if (phase.load(.acquire) != .ready) return;
    // if (vtable.quit) |f| f();
    _vtable = generic_vtable;
    phase.store(.uninit, .release);
}

pub fn coreName() []const u8 {
    return @tagName(core_type);
}

// lives here, not in `Vtable.zig`, so that file stays import-free for the
// per-core kernel objects, which use its type (see `kernels/__root__.zig`)
test "VTable default-initializes with every kernel pointer null" {
    const vt: Vtable = .{};
    try zigstd.testExpectTrue(vt.sgemm_kernel == null);
    try zigstd.testExpectTrue(vt.dgemm_kernel == null);
    try zigstd.testExpectTrue(vt.cgemm_kernel_n == null);
    try zigstd.testExpectTrue(vt.zgemm_kernel_n == null);
    try zigstd.testExpectTrue(vt.init == null);
    try zigstd.testExpectEqual(@as(i32, 0), vt.dtb_entries);
}

test {
    _ = @import("_test_types.zig");
    zigstd.testRefAllDecls(@This());
}
