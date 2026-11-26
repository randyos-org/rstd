//! Per-core gemm blocking parameters, extracted from OpenBLAS `v0.3.34-794-gdbf2c24c9`'s
//! `param.h` (`BUFFER_SIZE` from `common_<arch>.h`), one file per
//! architecture family, by `tools/blas_asm`'s `params.zig` (retired: issue #618);
//! maintained by hand. `.dynamic` and `.from_buffer`
//! mark the values the header leaves to `driver/others/parameter.c`
//! (the converter's `params.zig`); a field at its default in `Shape` or `Table` is
//! left out.

const x86_64_data = @import("x86_64.zig");
const x86_data = @import("x86.zig");
const arm64_data = @import("arm64.zig");
const arm_data = @import("arm.zig");
const power64_data = @import("power64.zig");
const power64be_data = @import("power64be.zig");
const power32_data = @import("power32.zig");
const power32le_data = @import("power32le.zig");
const riscv64_data = @import("riscv64.zig");
const s390x_data = @import("s390x.zig");
const mips64_data = @import("mips64.zig");
const mips32_data = @import("mips32.zig");
const loongarch64_data = @import("loongarch64.zig");
const sparc_data = @import("sparc.zig");
const wasm_data = @import("wasm.zig");
const generic_data = @import("generic.zig");

pub const x86_64 = x86_64_data.x86_64;
pub const x86_64_dtb = x86_64_data.x86_64_dtb;
pub const x86_64_buffer_size = x86_64_data.x86_64_buffer_size;
pub const x86 = x86_data.x86;
pub const x86_dtb = x86_data.x86_dtb;
pub const x86_buffer_size = x86_data.x86_buffer_size;
pub const arm64 = arm64_data.arm64;
pub const arm64_dtb = arm64_data.arm64_dtb;
pub const arm64_buffer_size = arm64_data.arm64_buffer_size;
pub const arm = arm_data.arm;
pub const arm_dtb = arm_data.arm_dtb;
pub const arm_buffer_size = arm_data.arm_buffer_size;
pub const power64 = power64_data.power64;
pub const power64_dtb = power64_data.power64_dtb;
pub const power64_buffer_size = power64_data.power64_buffer_size;
pub const power64be = power64be_data.power64be;
pub const power64be_dtb = power64be_data.power64be_dtb;
pub const power64be_buffer_size = power64be_data.power64be_buffer_size;
pub const power32 = power32_data.power32;
pub const power32_dtb = power32_data.power32_dtb;
pub const power32_buffer_size = power32_data.power32_buffer_size;
pub const power32le = power32le_data.power32le;
pub const power32le_dtb = power32le_data.power32le_dtb;
pub const power32le_buffer_size = power32le_data.power32le_buffer_size;
pub const riscv64 = riscv64_data.riscv64;
pub const riscv64_dtb = riscv64_data.riscv64_dtb;
pub const riscv64_buffer_size = riscv64_data.riscv64_buffer_size;
pub const s390x = s390x_data.s390x;
pub const s390x_dtb = s390x_data.s390x_dtb;
pub const s390x_buffer_size = s390x_data.s390x_buffer_size;
pub const mips64 = mips64_data.mips64;
pub const mips64_dtb = mips64_data.mips64_dtb;
pub const mips64_buffer_size = mips64_data.mips64_buffer_size;
pub const mips32 = mips32_data.mips32;
pub const mips32_dtb = mips32_data.mips32_dtb;
pub const mips32_buffer_size = mips32_data.mips32_buffer_size;
pub const loongarch64 = loongarch64_data.loongarch64;
pub const loongarch64_dtb = loongarch64_data.loongarch64_dtb;
pub const loongarch64_buffer_size = loongarch64_data.loongarch64_buffer_size;
pub const sparc = sparc_data.sparc;
pub const sparc_dtb = sparc_data.sparc_dtb;
pub const sparc_buffer_size = sparc_data.sparc_buffer_size;
pub const wasm = wasm_data.wasm;
pub const wasm_dtb = wasm_data.wasm_dtb;
pub const wasm_buffer_size = wasm_data.wasm_buffer_size;
pub const generic = generic_data.generic;
pub const generic_dtb = generic_data.generic_dtb;
pub const generic_buffer_size = generic_data.generic_buffer_size;
pub const alpha_buffer_size = generic_data.alpha_buffer_size;
pub const csky_buffer_size = generic_data.csky_buffer_size;
pub const generic_off_x86 = generic_data.generic_off_x86;
