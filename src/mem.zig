const log = @import("std").log.scoped(.memory);

// 12 = log2(4096), page size shared by x86_64/aarch64
pub const page_shift: usize = 12;
pub const page_size: usize = 1 << page_shift;
pub const page_mask: usize = page_size - 1;
