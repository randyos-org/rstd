//! `rstd.math.signal` (`../design.md`, #612), for audio first: window
//! functions (`kernels.dsp`'s), the STFT and uniformly partitioned
//! convolution over `kernels.fft`.

const dsp = @import("../kernels/dsp/__root__.zig");
const stft = @import("stft.zig");
const convolve = @import("convolve.zig");

pub const Symmetry = dsp.Symmetry;
pub const Window = dsp.Window;
pub const WindowError = dsp.WindowError;
pub const window = dsp.window;

pub const Stft = stft.Stft;
pub const StftInitError = stft.InitError;
pub const StftError = stft.Error;

pub const Convolver = convolve.Convolver;
pub const ConvolverInitError = convolve.InitError;
pub const ConvolverError = convolve.Error;

test {
    _ = stft;
    _ = convolve;
}
