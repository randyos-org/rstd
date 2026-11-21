//! Signal processing (rstd's design.md, "Planned: signal processing"), for
//! audio first: window functions, the STFT and uniformly partitioned
//! convolution over `rstd.math.kernels.fft`.

const window_ = @import("window.zig");
const stft = @import("stft.zig");
const convolve = @import("convolve.zig");

pub const Symmetry = window_.Symmetry;
pub const Window = window_.Window;
pub const WindowError = window_.Error;
pub const window = window_.window;

pub const Stft = stft.Stft;
pub const StftInitError = stft.InitError;
pub const StftError = stft.Error;

pub const Convolver = convolve.Convolver;
pub const ConvolverInitError = convolve.InitError;
pub const ConvolverError = convolve.Error;

test {
    _ = window_;
    _ = stft;
    _ = convolve;
}
