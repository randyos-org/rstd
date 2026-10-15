//! CPU detection, one module per architecture family -- the same split the
//! reference has across `driver/others/dynamic*.c`.
//!
//! Each module exposes a `CoreType` enum and a `detect()`. Which one is in
//! play is decided at compile time from `builtin.cpu.arch`, so only the
//! detector for the target being built ever gets analyzed. That matters:
//! these modules contain architecture-specific instructions (`cpuid`,
//! `mrs midr_el1`, `cpucfg`) that would not even assemble for a different
//! target.

const std = @import("std");
const sysinfo = @import("builtin");
const zigstd = @import("../../../zigstd.zig");
const Vtable = @import("../Vtable.zig");
const corespec = @import("_corespec.zig");
const registry = @import("../asm_registry.zig");
pub const select = @import("_select.zig");

// pub const auxv = @import("_auxv.zig");
pub const arm = @import("arm.zig");
pub const arm64 = @import("arm64.zig");
pub const generic = @import("generic.zig");
pub const loongarch = @import("loongarch.zig");
pub const mips = @import("mips.zig");
pub const power = @import("power.zig");
pub const riscv = @import("riscv.zig");
pub const s390x = @import("s390x.zig");
pub const sparc = @import("sparc.zig");
pub const wasm = @import("wasm.zig");
pub const x86 = @import("x86.zig");
pub const x86_64 = @import("x86_64.zig");

/// Panics with a standard message for a ported function that has no
/// implementation yet for the running build's `(os, arch)`. Call with
/// `@src()` from the call site (inside the dispatching `switch`) so the
/// message names the actual function that hit the unimplemented branch.
pub fn unsupported(src: zigstd.SrcInfo) noreturn {
    std.debug.panic(
        "{s}: no implementation for {s}-{s} yet ({s}:{d})",
        .{ src.fn_name, @tagName(sysinfo.os.tag), @tagName(sysinfo.cpu.arch), src.file, src.line },
    );
}

pub fn GetCpuArch(a: std.Target.Cpu.Arch) type {
    return switch (a) {
        .x86_64 => x86_64,
        .x86 => x86,
        .aarch64, .aarch64_be => arm64,
        .arm, .armeb, .thumb, .thumbeb => arm,
        .powerpc, .powerpcle, .powerpc64, .powerpc64le => power,
        .riscv32, .riscv32be, .riscv64, .riscv64be => riscv,
        .s390x => s390x,
        .mips, .mipsel, .mips64, .mips64el => mips,
        .loongarch32, .loongarch64 => loongarch,
        .sparc, .sparc64 => sparc,
        .wasm32, .wasm64 => wasm,
        else => generic,
    };
}

/// The detector module for the architecture being compiled for.
pub const arch = GetCpuArch(sysinfo.cpu.arch);

/// The core type enum for this target. Which variants exist depends on the
/// architecture -- there is no cross-architecture core type, just as the
/// reference has a separate `corename[]` table per `dynamic_*.c`.
pub const CoreType = arch.CoreType;

pub const vtableFromCoreType = arch.vtableFromCoreType;
pub const kernel = arch.kernel;

/// Detects the running CPU. Never fails: every architecture has a generic
/// fallback for the case where nothing can be determined.
pub fn detect() CoreType {
    return arch.detect();
}

/// The architecture's detector (`_detect/<arch>.zig`): `Caps`, `probe`,
/// `decide`, `runnable`.
pub const detector = corespec.Detector(sysinfo.cpu.arch);

/// What `choose` decided: the core the CPU or the user asked for, and the core
/// whose kernels serve it (null: the portable ones).
pub const Chosen = select.Chosen(CoreType);

/// The core to run, as `gotoblas_dynamic_init` picks it: `forced`
/// (`OPENBLAS_CORETYPE`) if it names a core the CPU can run, else detection,
/// each served by the closest core with a kernel object in this build.
/// `hwcap_mask` is `LD_HWCAP_MASK` (s390x only). `log.warn(level, fmt, args)`
/// reports what was not as asked. See `select.choose`.
pub fn choose(forced: ?[]const u8, hwcap_mask: ?u64, log: anytype) Chosen {
    const caps = if (comptime @hasDecl(detector, "probeMasked"))
        detector.probeMasked(hwcap_mask orelse ~@as(u64, 0))
    else
        detector.probe();
    return select.choose(detector, caps, forced, registry.LinkedCtx(CoreType){}, log);
}

/// Resolves a core type from an `OPENBLAS_CORETYPE` value, case
/// insensitively. Returns null when the name is not one this architecture
/// knows (or is the placeholder `unknown`), which the reference reports as
/// "Core not found" before falling back.
pub fn fromName(name: []const u8) ?CoreType {
    return select.byName(CoreType, name);
}

test "fromName round-trips every selectable core type for this architecture" {
    inline for (@typeInfo(CoreType).@"enum".field_names) |field_name| {
        const want: ?CoreType = if (select.selectable(@field(CoreType, field_name))) @field(CoreType, field_name) else null;
        try std.testing.expectEqual(want, fromName(field_name));
    }
}

test "fromName is case insensitive and rejects unknown names" {
    const first = comptime for (@typeInfo(CoreType).@"enum".field_names) |n| {
        if (select.selectable(@field(CoreType, n))) break n;
    } else unreachable;
    var upper: [64]u8 = undefined;
    const up = std.ascii.upperString(upper[0..first.len], first);
    try std.testing.expect(fromName(up) != null);
    try std.testing.expect(fromName("definitely-not-a-core") == null);
    try std.testing.expect(fromName("unknown") == null);
}

test "detect returns something for whatever we are built for" {
    _ = detect();
}

const NoLog = struct {
    pub fn warn(_: NoLog, _: i32, comptime _: []const u8, _: anytype) void {}
};

test "choose gives a core the CPU can run, and its server is linked or portable" {
    const r = choose(null, null, NoLog{});
    // the placeholder is a CPU nothing names: no core of its own to run
    if (select.selectable(r.core)) try std.testing.expect(detector.runnable(r.core, detector.probe()));
    if (r.served) |k| try std.testing.expect(registry.isLinked(k));
}

test "corespec.kernelCores lists the cores with kernels" {
    const tag = sysinfo.cpu.arch;
    const R = @import("../kernels/__root__.zig").ForArch(tag);
    const actual = comptime blk: {
        @setEvalBranchQuota(100_000_000);
        var cores: []const corespec.CoreType(tag) = &.{};
        if (R.supports(tag)) for (@typeInfo(R.Core).@"enum".field_names) |n| {
            if (R.hasKernels(tag, @field(R.Core, n))) cores = cores ++ .{@field(R.Core, n)};
        };
        break :blk cores;
    };
    if (!std.mem.eql(corespec.CoreType(tag), actual, corespec.kernelCores(tag))) {
        std.debug.print("kernelCores(.{t}) should be:", .{tag});
        for (actual) |c| std.debug.print(" .{t},", .{c});
        std.debug.print("\n", .{});
        return error.TestUnexpectedResult;
    }
}

// every detector's decision tests run on any host: they feed synthetic values
test {
    _ = select;
    _ = @import("_auxv.zig");
    _ = @import("_detect/x86.zig");
    _ = @import("_detect/arm.zig");
    _ = @import("_detect/arm64.zig");
    _ = @import("_detect/riscv.zig");
    _ = @import("_detect/loongarch.zig");
    _ = @import("_detect/power.zig");
    _ = @import("_detect/s390x.zig");
    _ = @import("_detect/mips.zig");
    _ = @import("_detect/sparc.zig");
    _ = @import("_detect/wasm.zig");
    _ = @import("_detect/generic.zig");
}
