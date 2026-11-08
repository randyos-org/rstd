//! The `rstd_matrix` module where the compiler lacks `@Matrix` (Zig 0.17,
//! stock 0.18); see `native.zig`. rstd's own fallback code takes over.

pub const has_matrix = false;
