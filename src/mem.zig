// 12 = log2(4096), 4 KiB page size used on every target we build for
// (x86_64, aarch64, arm, riscv64, ...)
pub const page_shift: usize = 12;
pub const page_size: usize = 1 << page_shift;
pub const page_mask: usize = page_size - 1;
