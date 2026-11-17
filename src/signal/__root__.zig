//! Signal processing (rstd's design.md, "Planned: signal processing"), for
//! audio first: window functions and the STFT over `rstd.fft`; partitioned
//! convolution (N31) comes later.

const window_ = @import("window.zig");
const stft = @import("stft.zig");

pub const Symmetry = window_.Symmetry;
pub const Window = window_.Window;
pub const WindowError = window_.Error;
pub const window = window_.window;

pub const Stft = stft.Stft;
pub const StftInitError = stft.InitError;
pub const StftError = stft.Error;

test {
    _ = window_;
    _ = stft;
}
