//! RISC-V CPU detection, mirroring `driver/others/dynamic_riscv64.c`.
//!
//! The reference deliberately invokes the `riscv_hwprobe` syscall directly
//! rather than through glibc, because glibc's wrapper landed later than
//! the kernel support and it wants to work on kernels that have the RVV
//! backports but not the hwprobe patches. It falls back to `AT_HWCAP`'s
//! 'V' bit when the syscall is unavailable.

const std = @import("std");
const coretypes = @import("_coretypes.zig");
const sysinfo = @import("builtin");
const auxv = @import("_auxv.zig");

pub const CoreType = coretypes.RiscvCoreType;
pub const kernel = @import("../kernels/riscv.zig");

/// `AT_HWCAP` reports base ISA extensions as letter bits: bit 0 is 'A',
/// so 'V' is bit 21.
pub const HWCAP_ISA_V: usize = 1 << ('V' - 'A');

const RISCV_HWPROBE_KEY_IMA_EXT_0: i64 = 4;
const RISCV_HWPROBE_IMA_V: u64 = 1 << 2;
const RISCV_HWPROBE_KEY_CPUPERF_0: i64 = 5;

const Pair = extern struct { key: i64, value: u64 };

/// The syscall number is stable across RISC-V Linux; the reference
/// hardcodes it for the same reason we do.
const NR_riscv_hwprobe: usize = 258;

fn hwprobe(pairs: []Pair, flags: usize) bool {
    if (sysinfo.os.tag != .linux) return false;
    const rc = std.os.linux.syscall5(
        @fromBackingInt(@intCast(NR_riscv_hwprobe)),
        @intFromPtr(pairs.ptr),
        pairs.len,
        0, // cpu set size
        0, // cpu set pointer
        flags,
    );
    return @as(isize, @bitCast(rc)) == 0;
}

pub fn supportVector() bool {
    var pairs = [_]Pair{.{ .key = RISCV_HWPROBE_KEY_IMA_EXT_0, .value = 0 }};
    if (hwprobe(&pairs, 0)) return pairs[0].value & RISCV_HWPROBE_IMA_V != 0;
    // Older kernels: fall back to the hwcap letter bit.
    return auxv.hwcap() & HWCAP_ISA_V != 0;
}

pub fn detect() CoreType {
    if (!supportVector()) return .riscv64_generic;

    // The reference distinguishes VLEN>=256 from VLEN==128 to pick between
    // the two vector kernel sets. Without a reliable VLEN query we take the
    // conservative choice, which is correct on any vector-capable part --
    // just not always the fastest one.
    return .riscv64_zvl128b;
}

const Vtable = @import("../Vtable.zig");
const axpy = @import("../blas1/axpy/__root__.zig");
pub fn vtableFromCoreType(core_type: CoreType) Vtable {
    var vt: Vtable = .{};
    axpy.install(&vt, core_type);
    return vt;
}

test "detect runs on RISC-V targets" {
    if (comptime !sysinfo.cpu.arch.isRISCV()) return error.SkipZigTest;
    _ = detect();
}
