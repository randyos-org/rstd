//! `corespec.zig`'s and `machine.zig`'s tests.

const std = @import("std");
const coretypes = @import("../coretypes.zig");
const corespec = @import("../corespec.zig");

const Spec = corespec.Spec;
const x86 = corespec.x86;
const of = corespec.of;
const kernelCores = corespec.kernelCores;
const runsOn = corespec.runsOn;
const objectQuery = corespec.objectQuery;
const findModel = corespec.findModel;
const featureSet = corespec.featureSet;
const machine = corespec.machine;
const resolve = corespec.resolve;
const excess = corespec.excess;

/// Every model and feature named for `Core` exists in `Target` (a
/// `std.Target.<arch>` namespace).
fn check(comptime Target: type, comptime Core: type) !void {
    inline for (@typeInfo(Core).@"enum".field_names) |core_name| {
        inline for (.{ comptime of(@field(Core, core_name)), comptime machine(@field(Core, core_name)) }) |s| {
            if (s.model) |m| try std.testing.expect(@hasDecl(Target.cpu, m));
            for (s.features ++ s.without) |name| try std.testing.expect(std.meta.stringToEnum(Target.Feature, name) != null);
        }
    }
}

test "every model and feature named in a spec or machine exists for its architecture" {
    try check(std.Target.x86, coretypes.X86CoreType);
    try check(std.Target.aarch64, coretypes.Arm64CoreType);
    try check(std.Target.arm, coretypes.ArmCoreType);
    try check(std.Target.powerpc, coretypes.PowerCoreType);
    try check(std.Target.riscv, coretypes.RiscvCoreType);
    try check(std.Target.mips, coretypes.MipsCoreType);
    try check(std.Target.loongarch, coretypes.LoongarchCoreType);
    try check(std.Target.s390x, coretypes.S390xCoreType);
    try check(std.Target.sparc, coretypes.SparcCoreType);
    try check(std.Target.wasm, coretypes.WasmCoreType);
}

test "cores that share flags compare equal, others do not" {
    try std.testing.expect(x86(.haswell).eql(x86(.zen)));
    try std.testing.expect(x86(.piledriver).eql(x86(.excavator)));
    try std.testing.expect(!x86(.haswell).eql(x86(.sandybridge)));
    try std.testing.expect(!x86(.skylakex).eql(x86(.cooperlake)));
}

/// How many of `arch`'s kernel objects, on `os`, need a feature a CPU they
/// serve may lack (`excess`), each printed. The targets are the baseline and
/// every kernel core's own CPU (`machine`), the weakest a build may name.
fn countExcess(comptime arch: std.Target.Cpu.Arch, os_tag: std.Target.Os.Tag) usize {
    const Model = std.Target.Cpu.Model;
    const os: std.Target.Os = .{ .tag = os_tag, .version_range = .{ .none = {} } };
    const cores = kernelCores(arch);
    var bad: usize = 0;
    for (0..cores.len + 1) |t| {
        const m: Spec = if (t == 0) .{ .model = Model.baseline(arch, os).name } else machine(cores[t - 1]);
        const model = if (m.model) |name| findModel(arch, name) else Model.generic(arch);
        var add = featureSet(arch, m.features);
        // the 32-bit cores' CPUs as 64-bit targets
        if (arch == .x86_64) add.addFeatureSet(featureSet(arch, &.{"64bit"}));
        const query: std.Target.Query = .{
            .cpu_arch = arch,
            .os_tag = os_tag,
            .cpu_model = .{ .explicit = model },
            .cpu_features_add = add,
            .cpu_features_sub = featureSet(arch, m.without),
        };
        const cpu: std.Target.Cpu = .{ .arch = arch, .model = model, .features = resolve(arch, model, add, query.cpu_features_sub) };
        for (cores) |core| {
            if (!runsOn(cpu, core)) continue;
            const q = objectQuery(query, os, cpu, of(core));
            const object_model = switch (q.cpu_model) {
                .explicit => |e| e,
                else => Model.baseline(arch, os),
            };
            const ex = excess(core, arch, resolve(arch, object_model, q.cpu_features_add, q.cpu_features_sub), cpu.features);
            if (ex.isEmpty()) continue;
            bad += 1;
            std.debug.print("{t}-{t} for {s}: {t}'s object needs", .{ arch, os_tag, model.name, core });
            for (arch.allFeaturesList(), 0..) |f, i| {
                if (ex.isEnabled(@intCast(i))) std.debug.print(" {s}", .{f.name});
            }
            std.debug.print("\n", .{});
        }
    }
    return bad;
}

test "no kernel object needs a feature the CPUs of its cores may lack" {
    const Os = std.Target.Os.Tag;
    var bad: usize = 0;
    inline for (.{
        .{ .x86, &[_]Os{ .linux, .freestanding } },
        .{ .x86_64, &[_]Os{ .linux, .macos } },
        .{ .aarch64, &[_]Os{ .linux, .macos } },
        .{ .arm, &[_]Os{ .linux, .freestanding } },
        .{ .powerpc, &[_]Os{ .linux, .freestanding } },
        .{ .powerpc64, &[_]Os{ .linux, .freestanding } },
        .{ .powerpc64le, &[_]Os{.linux} },
        .{ .mips, &[_]Os{.linux} },
        .{ .mipsel, &[_]Os{.linux} },
        .{ .mips64, &[_]Os{ .linux, .openbsd } },
        .{ .mips64el, &[_]Os{.linux} },
        .{ .riscv64, &[_]Os{.linux} },
        .{ .loongarch64, &[_]Os{.linux} },
        .{ .s390x, &[_]Os{.linux} },
        .{ .sparc, &[_]Os{ .linux, .freestanding } },
        .{ .sparc64, &[_]Os{.linux} },
        .{ .wasm32, &[_]Os{.wasi} },
    }) |a| {
        for (a[1]) |os| bad += countExcess(a[0], os);
    }
    try std.testing.expectEqual(@as(usize, 0), bad);
}
