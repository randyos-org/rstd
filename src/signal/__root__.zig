//! Signal processing (rstd's design.md, "Planned: signal processing"), for
//! audio first: window functions; the STFT (N30) and partitioned
//! convolution (N31) come later, over `rstd.fft`.

const window_ = @import("window.zig");

pub const Symmetry = window_.Symmetry;
pub const Window = window_.Window;
pub const WindowError = window_.Error;
pub const window = window_.window;

test {
    _ = window_;
}
