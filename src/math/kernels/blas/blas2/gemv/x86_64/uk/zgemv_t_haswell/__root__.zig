//! `zgemv_t` block kernel, the assembly as written, a file per form (the source preprocessed with `CONJ`, `XCONJ`, both or neither), from OpenBLAS `kernel/x86_64/zgemv_t_microk_haswell-4.c` (`v0.3.34-794-gdbf2c24c9`); converted by `tools/blas_asm` (retired: issue #618) and maintained by hand.

const plain = @import("plain.zig");
const c = @import("c.zig");
const x = @import("x.zig");
const cx = @import("cx.zig");

pub const zgemv_kernel_4x4 = plain.zgemv_kernel_4x4;
pub const zgemv_kernel_4x2 = plain.zgemv_kernel_4x2;
pub const zgemv_kernel_4x1 = plain.zgemv_kernel_4x1;
pub const zgemv_kernel_4x4_c = c.zgemv_kernel_4x4_c;
pub const zgemv_kernel_4x2_c = c.zgemv_kernel_4x2_c;
pub const zgemv_kernel_4x1_c = c.zgemv_kernel_4x1_c;
pub const zgemv_kernel_4x4_x = x.zgemv_kernel_4x4_x;
pub const zgemv_kernel_4x2_x = x.zgemv_kernel_4x2_x;
pub const zgemv_kernel_4x1_x = x.zgemv_kernel_4x1_x;
pub const zgemv_kernel_4x4_cx = cx.zgemv_kernel_4x4_cx;
pub const zgemv_kernel_4x2_cx = cx.zgemv_kernel_4x2_cx;
pub const zgemv_kernel_4x1_cx = cx.zgemv_kernel_4x1_cx;
