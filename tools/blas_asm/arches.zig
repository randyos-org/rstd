//! The architectures the generators know: OpenBLAS's `kernel/<ob>` directory
//! and where rstd serves it. `gen.zig`'s `arch/<name>.zig`, `uk.zig`,
//! `survey.zig` and `params.zig` take their entries from here.

pub const Arch = struct {
    /// Under OpenBLAS's `kernel/`.
    ob: []const u8,
    /// Under each `blas<n>/<routine>/`; empty where rstd has no target.
    dir: []const u8 = "",
    /// The registry, `kernels/<registry>.zig`.
    registry: []const u8 = "",
    /// The core enum in `cpu/coretypes.zig`.
    cores: []const u8 = "",
};

pub const x86_64: Arch = .{ .ob = "x86_64", .dir = "x86_64", .registry = "x86_64", .cores = "X86CoreType" };
pub const x86: Arch = .{ .ob = "x86", .dir = "x86", .registry = "x86", .cores = "X86CoreType" };
pub const arm64: Arch = .{ .ob = "arm64", .dir = "arm64", .registry = "arm64", .cores = "Arm64CoreType" };
pub const arm: Arch = .{ .ob = "arm", .dir = "arm", .registry = "arm", .cores = "ArmCoreType" };
pub const power: Arch = .{ .ob = "power", .dir = "power", .registry = "power", .cores = "PowerCoreType" };
pub const mips64: Arch = .{ .ob = "mips64", .dir = "mips64", .registry = "mips", .cores = "MipsCoreType" };
/// 32-bit MIPS, served by the mips64 kernels.
pub const mips: Arch = .{ .ob = "mips", .dir = "mips64", .registry = "mips", .cores = "MipsCoreType" };
pub const riscv64: Arch = .{ .ob = "riscv64", .dir = "riscv64", .registry = "riscv", .cores = "RiscvCoreType" };
pub const loongarch64: Arch = .{ .ob = "loongarch64", .dir = "loongarch64", .registry = "loongarch", .cores = "LoongarchCoreType" };
pub const zarch: Arch = .{ .ob = "zarch", .dir = "s390x", .registry = "s390x", .cores = "S390xCoreType" };
pub const sparc: Arch = .{ .ob = "sparc", .dir = "sparc", .registry = "sparc", .cores = "SparcCoreType" };
pub const alpha: Arch = .{ .ob = "alpha" };
pub const ia64: Arch = .{ .ob = "ia64" };
pub const wasm: Arch = .{ .ob = "wasm", .dir = "wasm", .registry = "wasm", .cores = "WasmCoreType" };
pub const csky: Arch = .{ .ob = "csky" };
/// OpenBLAS's portable kernels, which every architecture without its own takes.
pub const generic: Arch = .{ .ob = "generic", .cores = "GenericCoreType" };

/// The architectures `survey.zig` lists, in its order.
pub const surveyed = [_]Arch{ x86_64, x86, arm64, arm, power, mips64, mips, riscv64, loongarch64, zarch, sparc, alpha, ia64, wasm };
