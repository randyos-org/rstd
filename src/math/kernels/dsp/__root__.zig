//! `rstd.math.kernels.dsp` (`../../design.md`, #612): vector and signal
//! primitives, vDSP's; the window functions so far.

const window_ = @import("window.zig");

pub const Symmetry = window_.Symmetry;
pub const Window = window_.Window;
pub const WindowError = window_.Error;
pub const window = window_.window;

test {
    _ = window_;
}
