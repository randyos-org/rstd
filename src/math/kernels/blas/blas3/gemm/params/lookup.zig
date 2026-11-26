//! Which parameter table and rules a core takes: its architecture family, its cores, the generic branch and the buffer sizes.

const std = @import("std");
const sysinfo = @import("builtin");
const Vtable = @import("../../../Vtable.zig");
const data = @import("../params_data/__root__.zig");
const coretypes = @import("../../../cpu/coretypes.zig");
const small = @import("../small.zig");
const proot = @import("../params.zig");

const Env = proot.Env;
const Table = proot.Table;
const install = proot.install;
const installRule = proot.installRule;
const ruleFor = proot.ruleFor;

// ---- lookup ----------------------------------------------------------------

/// The architecture families the data is extracted for. The 32-bit ones take
/// `Makefile.system`'s 32-bit fallbacks (x86's AVX cores Nehalem's entry, the
/// Bulldozer family and Zen Barcelona's, POWER8 to 10 POWER6's); POWER's
/// follow the byte order too (`PPC970`'s unrolls).
pub const Family = enum { x86_64, x86, arm64, arm, power64, power64be, power32, power32le, riscv64, s390x, mips64, mips32, loongarch64, sparc, wasm, generic };

/// The core enum of `family`.
pub fn CoreOf(comptime family: Family) type {
    return switch (family) {
        .x86_64, .x86 => coretypes.X86CoreType,
        .arm64 => coretypes.Arm64CoreType,
        .arm => coretypes.ArmCoreType,
        .power64, .power64be, .power32, .power32le => coretypes.PowerCoreType,
        .riscv64 => coretypes.RiscvCoreType,
        .s390x => coretypes.S390xCoreType,
        .mips64, .mips32 => coretypes.MipsCoreType,
        .loongarch64 => coretypes.LoongarchCoreType,
        .sparc => coretypes.SparcCoreType,
        .wasm => coretypes.WasmCoreType,
        .generic => coretypes.GenericCoreType,
    };
}

/// The extracted table of `core`; null where `param.h` has none for it.
pub fn tableFor(comptime family: Family, core: CoreOf(family)) ?Table {
    return switch (family) {
        inline else => |f| @field(data, @tagName(f))(core),
    };
}

/// `DTB_DEFAULT_ENTRIES` of `core` (`getarch.c`). Where `getarch.c` has no
/// entry, a build for the core takes what its `cpuid_<arch>.c` prints (riscv64:
/// 64, for every core), else the generic target's.
pub fn dtbFor(comptime family: Family, core: CoreOf(family)) i32 {
    const own = switch (family) {
        .x86_64 => data.x86_64_dtb(core),
        .x86 => data.x86_dtb(core),
        .arm64 => data.arm64_dtb(core),
        .arm => data.arm_dtb(core),
        .power64 => data.power64_dtb(core),
        .power64be => data.power64be_dtb(core),
        .power32 => data.power32_dtb(core),
        .power32le => data.power32le_dtb(core),
        .riscv64 => data.riscv64_dtb(core),
        .s390x => data.s390x_dtb(core),
        .mips64 => data.mips64_dtb(core),
        .mips32 => data.mips32_dtb(core),
        .loongarch64 => data.loongarch64_dtb(core),
        .sparc => data.sparc_dtb(core),
        .wasm => data.wasm_dtb(core),
        .generic => data.generic_dtb(core),
    };
    return own orelse switch (family) {
        .riscv64 => 64,
        else => data.generic_dtb(.generic).?,
    };
}

/// The family a target's `cpu.Arch` (`cpu/__root__.zig`) installs.
fn familyOf(arch: std.Target.Cpu.Arch, ptr_bits: u16) Family {
    return switch (arch) {
        .x86_64 => .x86_64,
        .x86 => .x86,
        .aarch64, .aarch64_be => .arm64,
        .arm, .armeb, .thumb, .thumbeb => .arm,
        .powerpc64le => .power64,
        .powerpc64 => .power64be,
        .powerpcle => .power32le,
        .powerpc => .power32,
        .riscv32, .riscv32be, .riscv64, .riscv64be => .riscv64,
        .s390x => .s390x,
        .mips, .mipsel, .mips64, .mips64el => if (ptr_bits == 64) .mips64 else .mips32,
        .loongarch32, .loongarch64 => .loongarch64,
        .sparc, .sparc64 => .sparc,
        .wasm32, .wasm64 => .wasm,
        else => .generic,
    };
}

/// The branch of `param.h`'s `GENERIC` block a static `TARGET=GENERIC`
/// build reads, and how its header variables get set: x86 by
/// `parameter.c:451`'s L2 rule; `ARCH_MIPS` (32-bit MIPS) and
/// `ARCH_LOONGARCH64` by fixed sets (`param.h:4456-4483`); every other
/// architecture keeps `parameter.c`'s initializers, `p = q = r = 128`
/// (`parameter.c:46-48`), as only x86's `blas_set_parameter` sets them for
/// `GENERIC` (`:760-954`).
pub const GenericBranch = enum { x86, mips, loongarch64, other };

/// `GenericBranch` of an architecture, as `c_check` names it from `ctest.c`:
/// mips64 is `ARCH_MIPS64` (`ctest.c:128-134`, `c_check:89-90`), loongarch32
/// has no `__loongarch64`; architectures with no OpenBLAS port read the
/// `#else` branch too.
pub fn genericBranch(arch: std.Target.Cpu.Arch) GenericBranch {
    return switch (arch) {
        .x86, .x86_64 => .x86,
        .mips, .mipsel => .mips,
        .loongarch64 => .loongarch64,
        else => .other,
    };
}

/// The branch a family's architecture reads, for a core with no table.
fn familyBranch(family: Family) GenericBranch {
    return switch (family) {
        .x86_64, .x86 => .x86,
        .mips32 => .mips,
        .loongarch64 => .loongarch64,
        .arm64, .arm, .power64, .power64be, .power32, .power32le, .riscv64, .s390x, .mips64, .sparc, .wasm, .generic => .other,
    };
}

/// Installs the generic table as a static `TARGET=GENERIC` build on `branch`
/// has it, with `buffer_size` as its `BUFFER_SIZE` (x86 derives `r` from it).
pub fn installGeneric(vt: *Vtable, comptime branch: GenericBranch, buffer_size: i64, env: Env) void {
    var table = (if (branch == .x86) data.generic(.generic) else data.generic_off_x86(branch)).?;
    table.buffer_size = buffer_size;
    installRule(vt, .generic, .generic, table, env, comptime if (branch == .x86) ruleFor(.generic, .generic) else .none);
}

/// `BUFFER_SIZE` of `family`'s architecture (`common_<arch>.h`).
fn familyBufferSize(family: Family) i64 {
    return switch (family) {
        .x86_64 => data.x86_64_buffer_size,
        .x86 => data.x86_buffer_size,
        .arm64 => data.arm64_buffer_size,
        .arm => data.arm_buffer_size,
        .power64 => data.power64_buffer_size,
        .power64be => data.power64be_buffer_size,
        .power32 => data.power32_buffer_size,
        .power32le => data.power32le_buffer_size,
        .riscv64 => data.riscv64_buffer_size,
        .s390x => data.s390x_buffer_size,
        .mips64 => data.mips64_buffer_size,
        .mips32 => data.mips32_buffer_size,
        .loongarch64 => data.loongarch64_buffer_size,
        .sparc => data.sparc_buffer_size,
        .wasm => data.wasm_buffer_size,
        .generic => data.generic_buffer_size,
    };
}

/// `BUFFER_SIZE` of the generic table on a target: its own, since OpenBLAS's
/// generic target sees its architecture's `common_<arch>.h` like any build.
/// Alpha and C-SKY have a header but no family; other architectures have no
/// OpenBLAS port and keep the generic family's x86_64 stand-in.
pub fn genericBufferSize(arch: std.Target.Cpu.Arch, ptr_bits: u16) i64 {
    return switch (arch) {
        .alpha => data.alpha_buffer_size,
        .csky => data.csky_buffer_size,
        else => familyBufferSize(familyOf(arch, ptr_bits)),
    };
}

/// The generic table's `BUFFER_SIZE` in this build.
pub const native_buffer_size = genericBufferSize(sysinfo.target.cpu.arch, sysinfo.target.ptrBitWidth());

/// The generic table's branch in this build.
pub const native_branch = genericBranch(sysinfo.target.cpu.arch);

/// Installs the parameters of `core`. A core `param.h` has no blocking for takes
/// the generic set and its rule, as OpenBLAS's generic target does on the
/// family's architecture, with the family's `BUFFER_SIZE`; the generic family
/// takes the target's branch and size.
pub fn installCore(vt: *Vtable, comptime family: Family, core: CoreOf(family), env: Env) void {
    switch (core) {
        inline else => |c| {
            defer vt.dtb_entries = dtbFor(family, c);
            defer small.installPermits(vt, small.permitFor(family, c), small.halfPermitFor(family, c));
            if (comptime family == .generic) {
                installGeneric(vt, native_branch, native_buffer_size, env);
            } else if (tableFor(family, c)) |table| {
                install(vt, family, c, table, env);
            } else {
                installGeneric(vt, comptime familyBranch(family), comptime familyBufferSize(family), env);
            }
        },
    }
}

/// `params_data`'s tables, for the tests (`_test/blas3/params.zig`).
pub const tables = data;
