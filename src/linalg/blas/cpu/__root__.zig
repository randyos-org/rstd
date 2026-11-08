//! CPU detection, one module per architecture family, as the reference splits
//! `driver/others/dynamic*.c`; only the target's detector is analyzed, since
//! detectors hold instructions (`cpuid`, `mrs midr_el1`) others cannot assemble.

const std = @import("std");
const sysinfo = @import("builtin");
const Vtable = @import("../Vtable.zig");
const coretypes = @import("coretypes.zig");
const corespec = @import("corespec.zig");
const registry = @import("../registry.zig");
const kernels = @import("../kernels/__root__.zig");
const portable = @import("../portable.zig");
const gemm_params = @import("../blas3/gemm/params.zig");
const gemv = @import("../blas2/gemv/__root__.zig");
const blas = @import("../__root__.zig");
pub const select = @import("select.zig");

// pub const auxv = @import("_auxv.zig");
pub const arm = Arch(.arm);
pub const arm64 = Arch(.arm64);
pub const generic = Arch(.generic);
pub const loongarch = Arch(.loongarch);
pub const mips = Arch(.mips);
pub const power = Arch(.power);
pub const riscv = Arch(.riscv);
pub const s390x = Arch(.s390x);
pub const sparc = Arch(.sparc);
pub const wasm = Arch(.wasm);
pub const x86 = Arch(.x86);
pub const x86_64 = Arch(.x86_64);

/// A family's module: its detector's core type, its kernel registry and its
/// vtable. Only the target's is analyzed, as only its detector assembles.
pub fn Arch(comptime family: coretypes.Family) type {
    return struct {
        // names shared with this file's own: qualified by `Self`
        const Self = @This();
        const family_detector = corespec.DetectorOf(family);

        pub const CoreType = family_detector.CoreType;
        pub const kernel = kernels.ForFamily(family);

        /// Whether this machine can run `core`'s kernels (for tests).
        pub fn canRun(core: Self.CoreType) bool {
            return select.canRun(family_detector, core);
        }

        pub fn vtableFromCoreType(core_type: Self.CoreType) Vtable {
            var vt: Vtable = .{};
            portable.install(&vt);
            // blocking first: a core's level-3 kernels set their own unrolls
            gemm_params.installCore(&vt, blockingOf(family), core_type, .fromConfig(blas.config));
            Self.kernel.install(&vt, core_type);
            // the x87 `q`/`x` gemv_n block by the core's DTB entries
            if (comptime family == .x86) gemv.installExtendedCore(&vt, core_type);
            portable.finish(&vt);
            return vt;
        }
    };
}

/// The blocking tables of `family`, of this word size and byte order where
/// those differ (power, mips).
fn blockingOf(comptime family: coretypes.Family) gemm_params.Family {
    return switch (family) {
        .x86_64 => .x86_64,
        .x86 => .x86,
        .arm64 => .arm64,
        .arm => .arm,
        .power => switch (sysinfo.target.cpu.arch) {
            .powerpc64le => .power64,
            .powerpc64 => .power64be,
            .powerpcle => .power32le,
            else => .power32,
        },
        .riscv => .riscv64,
        .s390x => .s390x,
        .mips => if (@sizeOf(usize) == 8) .mips64 else .mips32,
        .loongarch => .loongarch64,
        .sparc => .sparc,
        .wasm => .wasm,
        .generic => .generic,
    };
}

/// The module of `a`'s family.
pub fn GetCpuArch(a: std.Target.Cpu.Arch) type {
    return Arch(coretypes.familyOf(a));
}

/// The detector module for the architecture being compiled for.
pub const arch = GetCpuArch(sysinfo.target.cpu.arch);

/// The core type enum for this target. Which variants exist depends on the
/// architecture; there is no cross-architecture core type, just as the
/// reference has a separate `corename[]` table per `dynamic_*.c`.
pub const CoreType = arch.CoreType;

pub const vtableFromCoreType = arch.vtableFromCoreType;
pub const kernel = arch.kernel;

/// Detects the running CPU, under the configured `LD_HWCAP_MASK` (s390x).
/// Never fails: every architecture has a generic fallback for the case
/// where nothing can be determined.
pub fn detect() CoreType {
    return select.detect(detector, blas.config.ld_hwcap_mask);
}

/// The architecture's detector (`detect/<arch>.zig`): `Caps`, `probe`,
/// `decide`, `runnable`.
pub const detector = corespec.Detector(sysinfo.target.cpu.arch);

/// What `choose` decided: the core the CPU or the user asked for, and the core
/// whose kernels serve it (null: the portable ones).
pub const Chosen = select.Chosen(CoreType);

/// The core to run, as `gotoblas_dynamic_init` picks it: `forced`
/// (`OPENBLAS_CORETYPE`) if it names a core the CPU can run, else detection,
/// each served by the closest core with a kernel object in this build.
/// `hwcap_mask` is `LD_HWCAP_MASK` (s390x only). `log.warn(level, fmt, args)`
/// reports what was not as asked. See `select.choose`.
pub fn choose(forced: ?[]const u8, hwcap_mask: ?u64, log: anytype) Chosen {
    return select.choose(detector, select.probe(detector, hwcap_mask), forced, registry.LinkedCtx(CoreType){}, log);
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

test "choose gives a core the CPU can run, and its server is linked or portable" {
    const r = choose(null, null, select.testing.NoLog{});
    // the placeholder is a CPU nothing names: no core of its own to run
    if (select.selectable(r.core)) try std.testing.expect(detector.runnable(r.core, detector.probe()));
    if (r.served) |k| try std.testing.expect(registry.isLinked(k));
}

test "corespec.kernelCores lists the cores with kernels" {
    const tag = sysinfo.target.cpu.arch;
    const R = kernels.ForArch(tag);
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
    _ = @import("osinfo.zig");
    _ = @import("detect/x86.zig");
    _ = @import("detect/arm.zig");
    _ = @import("detect/arm64.zig");
    _ = @import("detect/riscv.zig");
    _ = @import("detect/loongarch.zig");
    _ = @import("detect/power.zig");
    _ = @import("detect/s390x.zig");
    _ = @import("detect/mips.zig");
    _ = @import("detect/sparc.zig");
    _ = @import("detect/wasm.zig");
    _ = @import("detect/generic.zig");
}
