//! `zgemv_n` block kernel, the assembly as written, a file per form (the source preprocessed with `CONJ`, `XCONJ`, both or neither), from OpenBLAS `kernel/x86_64/zgemv_n_microk_bulldozer-4.c` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const plain = @import("plain.zig");
const c = @import("c.zig");
const x = @import("x.zig");

pub const zgemv_kernel_4x4 = plain.zgemv_kernel_4x4;
pub const zgemv_kernel_4x2 = plain.zgemv_kernel_4x2;
pub const zgemv_kernel_4x1 = plain.zgemv_kernel_4x1;
pub const add_y_asm = plain.add_y_asm;
pub const zgemv_kernel_4x4_c = c.zgemv_kernel_4x4_c;
pub const zgemv_kernel_4x2_c = c.zgemv_kernel_4x2_c;
pub const zgemv_kernel_4x1_c = c.zgemv_kernel_4x1_c;
pub const add_y_c_asm = add_y_asm;
pub const zgemv_kernel_4x4_x = zgemv_kernel_4x4_c;
pub const zgemv_kernel_4x2_x = zgemv_kernel_4x2_c;
pub const zgemv_kernel_4x1_x = zgemv_kernel_4x1_c;
pub const add_y_x_asm = x.add_y_x_asm;
pub const zgemv_kernel_4x4_cx = zgemv_kernel_4x4;
pub const zgemv_kernel_4x2_cx = zgemv_kernel_4x2;
pub const zgemv_kernel_4x1_cx = zgemv_kernel_4x1;
pub const add_y_cx_asm = add_y_x_asm;
